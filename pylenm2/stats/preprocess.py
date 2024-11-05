import numpy as np
import pandas as pd
from tqdm import tqdm
import scipy.stats as stats
from statsmodels.nonparametric.smoothers_lowess import lowess

from pylenm2.stats import metrics
from pylenm2.stats import gp as stats_gp
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

preprocess_logger = logger_config.setup_logging(
    module_name=__name__,
    level=logging.INFO,
    # level=logging.DEBUG,
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
    try:
        data = data.drop(data.index[row_loc]).reindex(data.index)     # NOTE: Reindexing is necessary to make sure that the size mismatch is handled by adding `NaN` values for the dropped rows.
    
    except ValueError as ve:
        preprocess_logger.error(ve)
        data = data.drop(data.index[row_loc])   # NOTE: Using this as a hack currently.
    
    except Exception as e:
        raise(e)

    # # New implementation
    # z = np.abs(stats.zscore(data, nan_policy=nan_policy))
    # data[z > z_threshold] = np.nan

    return data


## by K. Whiteaker 2024-08, kwhit@alum.mit.edu
def remove_outliers_lowess(
        # self, 
        data, 
        lowess_frac=0.1, 
        std_thresh=2.2, 
        return_difference=False,
    ):
    """Identifies outliers in time series data as deviations of std_thresh standard deviations from a LOWESS fit, then sets these outliers to np.nan

    Args:
        data (pd.Series, pd.DataFrame): data for the outliers to removed from. 
            MUST BE indexed by datetime (ex. collection date). 
            Can convert pylenm concentration dataframe (for one station and 
            one analyte) into a valid input via 
            dataframe.set_index('COLLECTION_DATE').RESULT. 
            Expects a single column of data if pd.DataFrame.
        lowess_frac (float, optional): fraction of total data points considered 
            in local fits of LOWESS smoother. A smaller value captures more 
            local behaviour, and may be required for large datasets. 
            Defaults to 0.1.
        std_thresh (int, optional): number of standard deviations in 
            (observation - fit) outside of which is considered an outlier. 
            Defaults to 2.2.
        return_difference (bool, optional): if True, return a pd.Series 
            containing the difference between data and lowess fit

    Returns:
        if return_difference:
            pd.Series: input data with outliers set to np.nan
            pd.Series: input data - lowess fit
        else:
            pd.Series: input data with outliers set to np.nan
    """

    # create a copy so the data isn't modified inplace
    working_data = data.copy()
    
    # lowess() reads datetime as nanoseconds and has issues with large numbers & low frac, so need to scale down x-axis during fitting
    scaleDown = 1e17
    x_data = pd.to_datetime(working_data.index)
    x_readable = x_data.astype(int).to_numpy()/scaleDown
    data_lowess = lowess(
        working_data, 
        x_readable, 
        frac=lowess_frac, 
        return_sorted=False,
    )

    # identify outlier locations and mark them as nan
    difference = working_data - data_lowess  # this is a pd.Series
    thresh = std_thresh * np.std(difference)
    difference_ignoreNaN = np.ma.array(difference, mask=np.isnan(difference)) # Use a mask to mark & ignore the NaNs as per https://stackoverflow.com/questions/37749900/how-to-disregard-the-nan-data-point-in-numpy-array-and-generate-the-normalized-d
    outliers_iloc = np.where(np.abs(difference_ignoreNaN)>thresh)[0]
    outliers = working_data.iloc[outliers_iloc]
    working_data.iloc[outliers_iloc] = np.nan

    if return_difference:
        return working_data, difference
    else:
        return working_data


# Helper fucntion for get_Best_Stations
def _get_Best_Station(
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

    num_selected = len(selected)
    errors = []
    
    if(model==None):
        if(len(selected)<5):
            model, pred = stats_gp.fit_gp(X, y, xx)
        else:
            model = None
    else:
        model=model
    
    if(verbose):  
        print("# of stations to choose from: ", len(leftover))
    
    if num_selected == 0:
        if(verbose): 
            print("Selecting first station")
        for ix in leftover:
            y_pred, r_map, residuals, lr_trend = stats_gp.interpolate_topo(
                X=X.iloc[ix:ix+1,:], 
                y=y[ix:ix+1], 
                xx=xx, 
                ft=ft, 
                regression=regression, 
                model=model, 
                smooth=smooth,
            )
            # y_err = stats_gp.mse(ref, y_pred)
            y_err = metrics.mse(ref, y_pred)
            errors.append((ix, y_err))
    
    if num_selected > 0:
        for ix in leftover:
            joined = selected + [ix]
            y_pred, r_map, residuals, lr_trend = stats_gp.interpolate_topo(
                X=X.iloc[joined,:], 
                y=y[joined], 
                xx=xx, 
                ft=ft, 
                regression=regression, 
                model=model, 
                smooth=smooth,
            )
            y_err = metrics.mse(ref, y_pred)
            errors.append((ix, y_err))
        
    err_ix = [x[0] for x in errors]
    err_vals = [x[1] for x in errors]
    min_val = min(err_vals)
    min_ix = err_ix[err_vals.index(min(err_vals))]
    
    preprocess_logger.info(f"Selected station: {min_ix} with a MSE error of {min_val}.")
    
    if(verbose):
        # print("Selected station: {} with a MSE error of {}\n".format(min_ix, min_val))
        print(f"Selected station: {min_ix} with a MSE error of {min_val}.")
    
    return min_ix, min_val


def get_Best_Stations(
        # self, 
        X, 
        y, 
        xx, 
        ref, 
        initial, 
        max_stations, 
        ft=['Elevation'], 
        regression='linear', 
        verbose=True, 
        smooth=True, 
        model=None,
    ):
    """Greedy optimization function to select a subset of stations as to minimizes the MSE from a reference map

    Args:
        X (numpy.array): array of dimension (number of stations, 2) where each element is a pair of UTM coordinates.
        y (numpy.array): array of size (number of stations) where each value corresponds to a concentration value at a station.
        xx (numpy.array): prediction locations
        ref (numpy.array): reference field to optimize for (aka best/true map)
        initial (list): indices of stations as the starting stations for optimization
        max_stations (int): number of stations to optimize for
        ft (list, optional): feature names to train on. Defaults to ['Elevation'].
        regression (str, optional): choice between 'linear' for linear regression, 'rf' for random forest regression, 'ridge' for ridge regression, or 'lasso' for lasso regression.. Defaults to 'linear'.
        verbose (bool, optional): v. Defaults to True.
        smooth (bool, optional): flag to toggle WhiteKernel on and off. Defaults to True.
        model (GaussianProcessRegressor, optional): model to fit. Defaults to None.

    Returns:
        list: index of best stations in order from best to worst
    """

    BAD_RETURN = None, None

    if not isinstance(initial, list):
        preprocess_logger.error(f"`initial` parameter must be a list. Found `{type(initial)}`.")
        return BAD_RETURN

    tot_err = []
    selected = initial
    leftover = list(range(0, X.shape[0])) # all indexes from 0 to number of station
    
    # Remove the initial set of stations from pool of station indices to choose from
    for i in initial:
        leftover.remove(i)

    station_itr_count = max_stations - len(selected)
    for i in tqdm(range(station_itr_count), desc="Station", total=station_itr_count):

        # For some reason both the if...else conditions contain the same piece of code!!!
        # if i == 0: # select first station will min error
        #     station_ix, err = _get_Best_Station(
        #         X=X, 
        #         y=y, 
        #         xx=xx, 
        #         ref=ref, 
        #         selected=selected, 
        #         leftover=leftover, 
        #         ft=ft, 
        #         regression=regression, 
        #         verbose=verbose, 
        #         smooth=smooth, 
        #         model=model,
        #     )
        #     selected.append(station_ix)
        #     leftover.remove(station_ix)
        #     tot_err.append(err)
        
        # else:
        station_ix, err = _get_Best_Station(
            X=X, 
            y=y, 
            xx=xx, 
            ref=ref, 
            selected=selected, 
            leftover=leftover, 
            ft=ft, 
            regression=regression, 
            verbose=verbose, 
            smooth=smooth, 
            model=model,
        )
        selected.append(station_ix)
        leftover.remove(station_ix)
        tot_err.append(err)
    
    preprocess_logger.info(f"{selected = }")
    
    return selected, tot_err
