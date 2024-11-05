import numpy as np

from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, WhiteKernel
from sklearn.linear_model import (
    LinearRegression, 
    Ridge, 
    Lasso, 
    RidgeCV, 
    LassoCV,
)

from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

gp_logger = logger_config.setup_logging(
    module_name=__name__,
    level=logging.INFO,
    # level=logging.DEBUG,
    logfile_dir=c.LOGFILE_DIR,
)


def get_Best_GP(
        # self, 
        X, 
        y, 
        smooth=True, 
        seed=42,
    ):
    """Returns the best Gaussian Process model for a given X and y.
    TODO: ASK: Ask if the bound values for `Matern` and `WhiteKernel` are chosen particularly or arbitrarily. @Zexuan
    
    Args:
        X (numpy.array): array of dimension (number of stations, 2) where each element is a pair of UTM coordinates.
        y (numpy.array): array of size (number of stations) where each value corresponds to a concentration value at a station.
        smooth (bool, optional): flag to toggle WhiteKernel on and off. Defaults to True.
        seed (int, optional): random state setting. Defaults to 42.

    Returns:
        GaussianProcessRegressor: best GP model
    """

    gp = GaussianProcessRegressor(normalize_y=True, random_state=seed, optimizer='fmin_l_bfgs_b')
    
    # Kernel models
    if(smooth):
        k1 = Matern(length_scale=400, nu=1.5, length_scale_bounds=(100.0, 5000.0)) + WhiteKernel() #noise_level_bounds=(1e-10, 1e5))
        k2 = Matern(length_scale=800, nu=1.5, length_scale_bounds=(100.0, 5000.0)) + WhiteKernel() #noise_level_bounds=(1e-10, 1e5))
        k3 = Matern(length_scale=1200, nu=1.5, length_scale_bounds=(100.0, 5000.0)) + WhiteKernel() #noise_level_bounds=(1e-10, 1e5))
        k4 = Matern(length_scale=400, nu=1.5, length_scale_bounds=(100.0, 5000.0)) + WhiteKernel() #noise_level_bounds=(1e-10, 1e5))
        k5 = Matern(length_scale=800, nu=1.5, length_scale_bounds=(100.0, 5000.0)) + WhiteKernel() #noise_level_bounds=(1e-10, 1e5))
        k6 = Matern(length_scale=1200, nu=2.5, length_scale_bounds=(100.0, 5000.0)) + WhiteKernel() #noise_level_bounds=(1e-10, 1e5))
        k7 = Matern(length_scale=400, nu=2.5, length_scale_bounds=(100.0, 5000.0)) + WhiteKernel() #noise_level_bounds=(1e-10, 1e5))
        k8 = Matern(length_scale=800, nu=2.5, length_scale_bounds=(100.0, 5000.0)) + WhiteKernel() #noise_level_bounds=(1e-10, 1e5))
        k9 = Matern(length_scale=1200, nu=2.5, length_scale_bounds=(100.0, 5000.0)) + WhiteKernel() #noise_level_bounds=(1e-10, 1e5))
        k10 = Matern(length_scale=400, nu=np.inf, length_scale_bounds=(100.0, 5000.0)) + WhiteKernel() #noise_level_bounds=(1e-10, 1e5))
        k11 = Matern(length_scale=800, nu=np.inf, length_scale_bounds=(100.0, 5000.0)) + WhiteKernel() #noise_level_bounds=(1e-10, 1e5))
        k12 = Matern(length_scale=1200, nu=np.inf, length_scale_bounds=(100.0, 5000.0)) + WhiteKernel() #noise_level_bounds=(1e-10, 1e5))
    else:
        k1 = Matern(length_scale=400, nu=1.5, length_scale_bounds=(100.0, 5000.0))
        k2 = Matern(length_scale=800, nu=1.5, length_scale_bounds=(100.0, 5000.0))
        k3 = Matern(length_scale=1200, nu=1.5, length_scale_bounds=(100.0, 5000.0))
        k4 = Matern(length_scale=400, nu=1.5, length_scale_bounds=(100.0, 5000.0))
        k5 = Matern(length_scale=800, nu=1.5, length_scale_bounds=(100.0, 5000.0))
        k6 = Matern(length_scale=1200, nu=2.5, length_scale_bounds=(100.0, 5000.0))
        k7 = Matern(length_scale=400, nu=2.5, length_scale_bounds=(100.0, 5000.0))
        k8 = Matern(length_scale=800, nu=2.5, length_scale_bounds=(100.0, 5000.0))
        k9 = Matern(length_scale=1200, nu=2.5, length_scale_bounds=(100.0, 5000.0))
        k10 = Matern(length_scale=400, nu=np.inf, length_scale_bounds=(100.0, 5000.0))
        k11 = Matern(length_scale=800, nu=np.inf, length_scale_bounds=(100.0, 5000.0))
        k12 = Matern(length_scale=1200, nu=np.inf, length_scale_bounds=(100.0, 5000.0))
    
    parameters = {'kernel': [k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12]}
    
    model = GridSearchCV(gp, parameters)
    model.fit(X, y)
    
    return model


def fit_gp(
        # self, 
        X, 
        y, 
        xx, 
        model=None, 
        smooth=True,
    ):
    """Fits Gaussian Process for X and y and returns both the GP model and the predicted values

    Args:
        X (numpy.array): array of dimension (number of stations, 2) where each element is a pair of UTM coordinates.
        y (numpy.array): array of size (number of stations) where each value corresponds to a concentration value at a station.
        xx (numpy.array): prediction locations
        model (GaussianProcessRegressor, optional): model to fit. Defaults to None.
        smooth (bool, optional): flag to toggle WhiteKernel on and off. Defaults to True.

    Returns:
        GaussianProcessRegressor, numpy.array: GP model, prediction of xx
    """
    
    if(model==None):
        gp = get_Best_GP(X, y, smooth) # selects best kernel params to fit
    else:
        gp = model
    
    gp.fit(X, y)
    y_pred = gp.predict(xx)
    
    return gp, y_pred


def interpolate_topo(
        # self, 
        X, 
        y, 
        xx, 
        ft=['Elevation'], 
        model=None, 
        smooth=True, 
        regression='linear', 
        seed = 42,
    ):
    """Spatially interpolate the water table as a function of topographic metrics using Gaussian Process. Uses regression to generate trendline adds the values to the GP map.

    Args:
        X (numpy.array): training values. Must include "Easting" and "Northing" columns.
        y (numpy.array): array of size (number of stations) where each value corresponds to a concentration value at a station.
        xx (numpy.array): prediction locations
        ft (list, optional): eature names to train on. Defaults to ['Elevation'].
        model (GaussianProcessRegressor, optional): model to fit. Defaults to None.
        smooth (bool, optional): flag to toggle WhiteKernel on and off. Defaults to True.
        regression (str, optional): choice between 'linear' for linear regression, 'rf' for random forest regression, 'ridge' for ridge regression, or 'lasso' for lasso regression.. Defaults to 'linear'.
        seed (int, optional): random state setting. Defaults to 42.

    Returns:
        numpy.array: predicton of locations xx
    """

    # `alpha_values` are only used for RidgeCV and LassoCV models
    alpha_Values = [1e-5, 5e-5, 0.0001, 0.0005, 0.005, 0.05, 0.1, 0.3, 1, 3, 5, 10, 15, 30, 80]
    
    if(regression.lower()=='linear'):
        reg = LinearRegression()
    
    if(regression.lower()=='rf'):
        reg = RandomForestRegressor(n_estimators=200, random_state=seed)
    
    if(regression.lower()=='ridge'):
        # reg = make_pipeline(PolynomialFeatures(3), Ridge())
        reg = RidgeCV(alphas=alpha_Values)
    
    if(regression.lower()=='lasso'):
        # reg = make_pipeline(PolynomialFeatures(3), Lasso())
        reg = LassoCV(alphas=alpha_Values)
    
    if(all(elem in list(xx.columns) for elem in ft)):
        reg.fit(X[ft], y)
        y_est = reg.predict(X[ft])
        residuals = y - y_est
        if(model==None):
            model = get_Best_GP(X[['Easting','Northing']], residuals, smooth=smooth, seed=seed)
        else:
            model = model
        reg_trend = reg.predict(xx[ft])
    else:
        reg.fit(X[['Easting','Northing']], y)
        y_est = reg.predict(X[['Easting','Northing']])
        residuals = y - y_est
        if(model==None):
            model = get_Best_GP(X[['Easting','Northing']], residuals, smooth=smooth, seed=seed)
        else:
            model = model
        reg_trend = reg.predict(xx)
    
    r_map = model.predict(xx[['Easting','Northing']])
    y_map = reg_trend + r_map
    
    return y_map, r_map, residuals, reg_trend