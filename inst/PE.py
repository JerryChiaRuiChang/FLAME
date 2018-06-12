import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import mean_squared_error as MSE
from sklearn import linear_model
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score
from sklearn.tree import DecisionTreeRegressor

def predictive_error(holdout, covs, covs_subset, ridge_reg = 0.1):

	#Change column names into panda index (object)
    col = list(range(len(covs)))
    col.extend(["outcome","treated"])
    col = pd.Index(col)
    holdout.columns = col

    ridge_c = Ridge(alpha=ridge_reg) 
    ridge_t = Ridge(alpha=ridge_reg) 


    holdout_treated = holdout[holdout['treated']==1][covs_subset]
    holdout_control = holdout[holdout['treated']==0][covs_subset]


    mse_t = np.mean(cross_val_score(ridge_t, holdout_treated, 
                                holdout[holdout['treated']==1]['outcome'] , scoring = 'neg_mean_squared_error' ) )


    mse_c = np.mean(cross_val_score(ridge_c, holdout_control, 
                                holdout[holdout['treated']==0]['outcome'] , scoring = 'neg_mean_squared_error' ) )

     
    return ((mse_t + mse_c))
