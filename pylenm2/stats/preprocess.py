import numpy as np
import scipy.stats as stats

from pylenm2.stats import metrics
from pylenm2.stats import gp as stats_gp
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

preprocess_logger = logger_config.setup_logging(
    module_name=__name__,
    # level=logging.INFO,
    level=logging.DEBUG,
    logfile_dir=c.LOGFILE_DIR,
)



# Helper function for plot_correlation
# Sorts analytes in a specific order: 'TRITIUM', 'URANIUM-238','IODINE-129','SPECIFIC CONDUCTANCE', 'PH', 'DEPTH_TO_WATER'
def _custom_analyte_sort(analytes):
    my_order = 'TURISPDABCEFGHJKLMNOQVWXYZ-_abcdefghijklmnopqrstuvwxyz135790 2468'
    return sorted(
        analytes, 
        key=lambda word: [
            my_order.index(c) for c in word
        ],
    )


def get_MCL(analyte_name):
    """Returns the Maximum Concentration Limit value for the specified analyte. Example: 'TRITIUM' returns 1.3

    Args:
        analyte_name (str): name of the analyte to be processed

    Returns:
        float: MLC value
    """
    mcl_dictionary = {
        'TRITIUM': 1.3, 
        'URANIUM-238': 1.31, 
        'NITRATE-NITRITE AS NITROGEN': 1,
        'TECHNETIUM-99': 2.95, 
        'IODINE-129': 0, 
        'STRONTIUM-90': 0.9,
    }
    return mcl_dictionary[analyte_name]


def remove_outliers(
        data, 
        z_threshold=4,
        nan_policy="propagate",     # same as the default of stats.zscore()
    ):
    """Removes outliers from a dataframe based on the z_scores and returns the new dataframe.
    NOTE: The new logic is not same as the previous logic. 
        As per the new logic, it sets NA for the values that are greater than 
        z_threshold, whereas, in the previous logic, all the rows containing 
        even a single column of z > threshold are dropped (and then NaNs are added internally for whatever rows are dropped.)
    TODO: Confirm if this change is okay or if we need to revert back to the previous implementation. @Zexuan
    
    TODO: REWRITE!!! 
        REMOVE `nan_policy`.
        TEST AGAIN AFTER REWRITING TO MAKE SURE THAT THE OUTPUTS MATCH.

    Args:
        data (pd.DataFrame): data for the outliers to removed from
        z_threshold (int, optional): z_score threshold to eliminate. Values above this threshold are elimited. Defaults to 4. 
            NOTE: Get the threshold and logic confirmed by the @Zexuan and @Haruko.
        nan_policy (str, optional): specifies how to handle `nan` values. Passed in the stats.zscore() function.
            Options are one of:
                'propagate': returns nan.
                'raise': throws an error.
                'omit': performs the calculations ignoring nan values.
            Default is 'omit'.

    Returns:
        pd.DataFrame: data with outliers removed
    """

    # Old implementation
    z = np.abs(stats.zscore(data, nan_policy=nan_policy))
    row_loc = np.unique(np.where(z > z_threshold)[0])
    
    # Setting values outside threshold to `nan` values.
    # data = data.drop(data.index[row_loc])
    data = data.drop(data.index[row_loc]).reindex(data.index)     # NOTE: Reindexing is necessary to make sure that the size mismatch is handled by adding `NaN` values for the dropped rows.

    # # New implementation
    # z = np.abs(stats.zscore(data, nan_policy=nan_policy))
    # data[z > z_threshold] = np.nan

    return data


# Helper fucntion for get_Best_Wells
def _get_Best_Well(
        # self, 
        X, 
        y, 
        xx, 
        ref, 
        selected, 
        leftover, 
        ft=['Elevation'], 
        regression='linear', 
        verbose=True, 
        smooth=True, 
        model=None,
    ):

    num_selected=len(selected)
    errors = []
    if(model==None):
        if(len(selected)<5):
            model, pred = stats_gp.fit_gp(X, y, xx)
        else:
            model = None
    else:
        model=model
    if(verbose):  
        print("# of wells to choose from: ", len(leftover))
    if(num_selected==0):
        if(verbose): 
            print("Selecting first well")
        for ix in leftover:
            y_pred, r_map, residuals, lr_trend = stats_gp.interpolate_topo(X=X.iloc[ix:ix+1,:], y=y[ix:ix+1], xx=xx, ft=ft, regression=regression, model=model, smooth=smooth)
            y_err = stats_gp.mse(ref, y_pred)
            errors.append((ix, y_err))
    
    if(num_selected > 0):
        for ix in leftover:
            joined = selected + [ix]
            y_pred, r_map, residuals, lr_trend = stats_gp.interpolate_topo(X=X.iloc[joined,:], y=y[joined], xx=xx, ft=ft, regression=regression, model=model, smooth=smooth)
            y_err = metrics.mse(ref, y_pred)
            errors.append((ix, y_err))
        
    err_ix = [x[0] for x in errors]
    err_vals = [x[1] for x in errors]
    min_val = min(err_vals)
    min_ix = err_ix[err_vals.index(min(err_vals))]
    if(verbose):
        print("Selected well: {} with a MSE error of {}\n".format(min_ix, min_val))
    return min_ix, min_val


def get_Best_Wells(
        # self, 
        X, 
        y, 
        xx, 
        ref, 
        initial, 
        max_wells, 
        ft=['Elevation'], 
        regression='linear', 
        verbose=True, 
        smooth=True, 
        model=None,
    ):
    """Greedy optimization function to select a subset of wells as to minimizes the MSE from a reference map

    Args:
        X (numpy.array): array of dimension (number of wells, 2) where each element is a pair of UTM coordinates.
        y (numpy.array): array of size (number of wells) where each value corresponds to a concentration value at a well.
        xx (numpy.array): prediction locations
        ref (numpy.array): reference field to optimize for (aka best/true map)
        initial (list): indices of wells as the starting wells for optimization
        max_wells (int): number of wells to optimize for
        ft (list, optional): feature names to train on. Defaults to ['Elevation'].
        regression (str, optional): choice between 'linear' for linear regression, 'rf' for random forest regression, 'ridge' for ridge regression, or 'lasso' for lasso regression.. Defaults to 'linear'.
        verbose (bool, optional): v. Defaults to True.
        smooth (bool, optional): flag to toggle WhiteKernel on and off. Defaults to True.
        model (GaussianProcessRegressor, optional): model to fit. Defaults to None.

    Returns:
        list: index of best wells in order from best to worst
    """
    tot_err = []
    selected = initial
    leftover = list(range(0, X.shape[0])) # all indexes from 0 to number of well
    
    # Remove the initial set of wells from pool of well indices to choose from
    for i in initial:
        leftover.remove(i)

    for i in range(max_wells-len(selected)):
        if(i==0): # select first well will min error
            well_ix, err = _get_Best_Well(X=X,y=y, xx=xx, ref=ref, selected=selected, leftover=leftover, ft=ft, regression=regression, verbose=verbose, smooth=smooth, model=model)
            selected.append(well_ix)
            leftover.remove(well_ix)
            tot_err.append(err)
        else:
            well_ix, err = _get_Best_Well(X=X,y=y, xx=xx, ref=ref, selected=selected, leftover=leftover, ft=ft, regression=regression, verbose=verbose, smooth=smooth, model=model)
            selected.append(well_ix)
            leftover.remove(well_ix)
            tot_err.append(err)
    print(selected)
    return selected, tot_err