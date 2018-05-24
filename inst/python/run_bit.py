import numpy as np
import pandas as pd
#import pyodbc
import pickle
import time
import itertools
from joblib import Parallel, delayed

from sklearn.metrics.pairwise import pairwise_distances
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import mean_squared_error as MSE
from operator import itemgetter

import operator
from sklearn import linear_model

from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score
from sklearn.tree import DecisionTreeRegressor

def match(df, covs, covs_max_list, treatment_indicator_col = 'treated', match_indicator_col = 'matched'):
    
    # this function takes a dataframe, a set of covariates to match on, 
    # the treatment indicator column and the matched indicator column.
    # it returns the array indicating whether each unit is matched (the first return value), 
    # and a list of indices for the matched units (the second return value)

    #Get Data
    arr_slice_wo_t = df[covs].values # the covariates values as a matrix
    arr_slice_w_t = df[ covs + [treatment_indicator_col] ].values # the covariate values together with the treatment indicator as a matrix

    #Compute b_u and b_u+
    lidx_wo_t = np.dot( arr_slice_wo_t, np.array([ covs_max_list[i]**(len(covs_max_list) - 1 - i) for i in range(len(covs_max_list))]) ) # matrix multiplication, get a unique number for each unit
    lidx_w_t = np.dot( arr_slice_w_t, np.array([ covs_max_list[i]**(len(covs_max_list) - i) for i in range(len(covs_max_list))] +                                               [1]
                                              ) ) # matrix multiplication, get a unique number for each unit with treatment indicator

    #Compute c_u and c_u+
    _, unqtags_wo_t, counts_wo_t = np.unique(lidx_wo_t, return_inverse=True, return_counts=True) # count how many times each number appears
    _, unqtags_w_t, counts_w_t = np.unique(lidx_w_t, return_inverse=True, return_counts=True) # count how many times each number appears (with treatment indicator)

    #A unit is matched if and only if the counts don't agree
    match_indicator = ~(counts_w_t[unqtags_w_t] == counts_wo_t[unqtags_wo_t]) 
    
    return match_indicator, lidx_wo_t[match_indicator]

# match_quality, the larger the better
def match_quality(df, holdout, covs_subset, match_indicator, ridge_reg = 0.1, tradeoff = 0.1):
    
    s = time.time() 

    #Calculate number of units unmatched (available)
    num_control = len(df[df['treated']==0]) # how many control units that are unmatched (recall matched units are removed from the data frame)
    num_treated = len(df[df['treated']==1]) # how many treated units that are unmatched (recall matched units are removed from the data frame)
    
    #Calculate number of units matched at current level
    num_control_matched = np.sum(( match_indicator ) & (df['treated']==0) ) # how many control units that are matched on this level
    num_treated_matched = np.sum(( match_indicator ) & (df['treated']==1) ) # how many treated units that are matched on this level
  
    time_BF = time.time() - s 

    #Calculate PE 
    s = time.time() #?  
    ridge_c = Ridge(alpha=ridge_reg) 
    ridge_t = Ridge(alpha=ridge_reg) 

    n_mse_t = np.mean(cross_val_score(ridge_t, holdout[holdout['treated']==1][covs_subset], 
                                holdout[holdout['treated']==1]['outcome'] , scoring = 'neg_mean_squared_error' ) )

    n_mse_c = np.mean(cross_val_score(ridge_c, holdout[holdout['treated']==0][covs_subset], 
                                holdout[holdout['treated']==0]['outcome'] , scoring = 'neg_mean_squared_error' ) )
    time_PE = time.time() - s # ? 


    # return level-wise match quality 
    return  (tradeoff * ( float(num_control_matched)/num_control + float(num_treated_matched)/num_treated ) +
             ( n_mse_t + n_mse_c ) , time_PE , time_BF)


def get_CATE_bit(df, match_indicator, index):

    d = df[ match_indicator ]
    if index is None: # when index == None, nothing is matched
        return None
    d.loc[:,'grp_id'] = index
    res = d.groupby(['grp_id', 'treated'])['outcome'].aggregate([np.size, np.mean]) # we do a groupby to get the statistics
    return res

def recover_covs(d, covs, covs_max_list, binary = True):

    ind = d.index.get_level_values(0)
    ind = [ num2vec(ind[i], covs_max_list) for i in range(len(ind)) if i%2==0]
    
    df = pd.DataFrame(ind, columns=covs ).astype(int)

    mean_list = list(d['mean'])
    size_list = list(d['size'])
        
    effect_list = [mean_list[2*i+1] - mean_list[2*i] for i in range(len(mean_list)//2) ]
    df.loc[:,'effect'] = effect_list
    df.loc[:,'size'] = [size_list[2*i+1] + size_list[2*i] for i in range(len(size_list)//2) ]
    
    return df

def cleanup_result(res_all):
    res = []
    for i in range(len(res_all)):
        r = res_all[i]
        if not r[1] is None:
            res.append(recover_covs( r[1], r[0][0], r[0][1] ) )
    return res

def num2vec(num, covs_max_list):
    res = []
    for i in reversed(range(len(covs_max_list))):
        res.append(num//pow(covs_max_list[i],i))
        num = num % pow(covs_max_list[i],i)
    return res

def run_bit(df, holdout, covs, covs_max_list, num_treated, num_control, tradeoff_param = 0.1):

    #Change column names into panda index (object)
    col = list(range(len(covs)))
    col.extend(["outcome","treated","matched"])
    col = pd.Index(col)
    df.columns = col
    holdout.columns = col

    #Change row names into panda index (int64)
    row = list(range(num_control))
    row.extend(list(range(num_treated)))
    row = pd.Index(row)
    df.index = row
    holdout.index = row

    constant_list = ['outcome', 'treated']
    
    covs_dropped = []
    cur_covs = covs[:]
    cur_covs_max_list = covs_max_list[:]

    timings = [0]*5 # first entry - match (matrix multiplication and value counting and comparison), 
                    # second entry - regression (compute PE),
                    # third entry - compute BF
                    # fourth entry - keep track of CATE,
                    # fifth entry - update dataframe (remove matched units)
    
    level = 0
    s = time.time()
    match_indicator, index = match(df, cur_covs, covs_max_list) # match without dropping anything
    timings[0] = timings[0] + time.time() - s
    
    s = time.time()
    res = get_CATE_bit(df, match_indicator, index) # get the CATEs without dropping anything
    timings[3] = timings[3] + time.time() - s
    
    matching_res = [[( cur_covs, cur_covs_max_list, None, match_indicator, index), res]] # result on first level, None says nothing is dropped
    
    s = time.time()
    df = df[~match_indicator][ cur_covs + constant_list ] # remove matched units
    timings[4] = timings[4] + time.time() - s
    
    level_scores = []
    
    while len(cur_covs)>1:
        
        print(cur_covs)
        
        best_score = np.inf
        level += 1
        matching_result_tmp = []
        
        if (np.sum(df['treated'] == 0) == 0 ) | (np.sum(df['treated'] == 1) == 0 ): # the early stopping condition
            print('no more matches')
            break
        
        for i in range(len(cur_covs)):
            
            cur_covs_no_c = cur_covs[:i] + cur_covs[i+1:]
            
            cur_covs_max_list_no_c = cur_covs_max_list[:i] + cur_covs_max_list[i+1:]
            
            s = time.time()
            match_indicator, index = match(df, cur_covs_no_c, cur_covs_max_list_no_c)
            timings[0] = timings[0] + time.time() - s 
            
            score, time_PE, time_BF = match_quality(df, holdout, cur_covs_no_c, match_indicator, tradeoff=tradeoff_param)
            timings[1] = timings[1] + time_PE 
            timings[2] = timings[2] + time_BF 
                                    
            matching_result_tmp.append( (cur_covs_no_c, cur_covs_max_list_no_c, score, match_indicator, index) )
        
        best_res = max(matching_result_tmp, key=itemgetter(2)) # use the one with largest MQ as the one to drop
        
        level_scores.append(max( [t[2] for t in matching_result_tmp] ))
        
        del matching_result_tmp
        
        new_matching_res = get_CATE_bit(df, best_res[-2], best_res[-1])
        
        cur_covs = best_res[0] 
        cur_covs_max_list = best_res[1]
        matching_res.append([best_res, new_matching_res])
        
        s = time.time()
        df = df[~ best_res[-2] ]
        timings[4] = timings[4] + time.time() - s
        
    return (cleanup_result(matching_res), level_scores)









