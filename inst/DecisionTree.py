import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error as MSE
from sklearn import tree
from sklearn.tree import DecisionTreeRegressor
from sklearn.model_selection import cross_val_score

def predictive_error(holdout, num_covs, covs_subset, tree_depth):

	# Change column names into panda index (object)
    col = list(range(num_covs))
    col.extend(["outcome","treated"])
    col = pd.Index(col)
    holdout.columns = col

    # Ridge Regression Model
    tree_c = DecisionTreeRegressor(max_depth = tree_depth, random_state = 0)
    tree_t = DecisionTreeRegressor(max_depth = tree_depth, random_state = 0)


    holdout_treated = holdout[holdout['treated']==1][covs_subset]
    holdout_control = holdout[holdout['treated']==0][covs_subset]


    mse_t = np.mean(cross_val_score(tree_t, holdout_treated,
                                holdout[holdout['treated']==1]['outcome'] , scoring = 'neg_mean_squared_error' ) )


    mse_c = np.mean(cross_val_score(tree_c, holdout_control,
                                holdout[holdout['treated']==0]['outcome'] , scoring = 'neg_mean_squared_error' ) )


    return ((mse_t + mse_c))
