import os
import scipy
import random
import itertools
import numpy as np
import pandas as pd
from tqdm import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import date2num, num2date

from supersmoother import SuperSmoother
from statsmodels.nonparametric.smoothers_lowess import lowess

from pylenm2.data import filters
from pylenm2.stats import preprocess
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

plots_logger = logger_config.setup_logging(
    module_name=__name__,
    level=logging.INFO,
    # level=logging.DEBUG,
    logfile_dir=c.LOGFILE_DIR,
)


# def _plotUpperHalf(*args, **kwargs):
def _plotUpperHalf(x_data, y_data, *args, **kwargs):   # TODO: confirm that this change does not break any other function.
    """<Function docstring>
    TODO: Complete the function docstring.

    Args:
        x_data (pd.Series): X-axis data.
        y_data (pd.Series): Y-axis data.
        args (Tuple): Tuple of positional arguments.
        kwargs (Dict): Dictionary of keyword arguments.
    """

    # corr_r = args[0].corr(args[1], 'pearson')
    corr_r = x_data.corr(y_data, 'pearson')
    corr_text = f"{corr_r:2.2f}"
    
    ax = plt.gca()
    ax.set_axis_off()
    
    marker_size = abs(corr_r) * 10000
    ax.scatter([.5], [.5], marker_size, [corr_r], alpha=0.6, cmap="coolwarm",
                vmin=-1, vmax=1, transform=ax.transAxes)
    
    font_size = abs(corr_r) * 40 + 5
    ax.annotate(corr_text, [.5, .48,],  xycoords="axes fraction", # [.5, .48,]
                ha='center', va='center', fontsize=font_size, fontweight='bold')


def plot_data(
        data_pylenm_dm, 
        station_name, 
        analyte_name, 
        log_transform=True, 
        alpha=0,
        plot_inline=True, 
        year_interval=2, 
        x_label='Years', 
        y_label='',
        save_plot=False, 
        save_dir='plot_data', 
        filter=False, 
        col=None, 
        equals=[],
    ):
    """Plot concentrations over time of a specified station and analyte with a smoothed curve on interpolated data points.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        station_name (str): name of the station to be processed
        analyte_name (str): name of the analyte to be processed
        log_transform (bool, optional): choose whether or not the data should 
            be transformed to log base 10 values. Defaults to `True`.
        alpha (int, optional): value between 0 and 10 for line smoothing. 
            Defaults to 0.
        plot_inline (bool, optional): `True` to show plot inline else `False`. 
            Defaults to `True`.
        year_interval (int, optional): plot by how many years to appear in the 
            axis e.g.(1 = every year, 5 = every 5 years, ...). Defaults to 2.
        x_label (str, optional): x axis label. Defaults to 'Years'.
        y_label (str, optional): y axis label. Defaults to ''.
        save_plot (bool, optional): `True` to save the plot to the `save_dir` 
            else `False`. Defaults to `False`.
        save_dir (str, optional): name of the directory you want to save the 
            plot to. `save_plot` must be `True` for this to be useful. 
            Defaults to 'plot_data'.
        filter (bool, optional): flag to indicate filtering. 
            Defaults to `False`.
        col (str, optional): column to filter. Example: col='STATION_ID'. 
            Defaults to None.
        equals (list, optional): values to filter col by. 
            Examples: equals=['FAI001A', 'FAI001B']. Defaults to [].

    Returns:
        None
    """

    # Gets appropriate data (station_name and analyte_name)
    query_raw = filters.query_data(
        data_pylenm_dm=data_pylenm_dm, 
        station_name=station_name, 
        analyte_name=analyte_name,
    )       # TODO: Handle better for empty values.

    if query_raw is None:
        plots_logger.debug(f"No results found for {station_name} and {analyte_name}")
        return None
    
    query = filters.simplify_data(data=query_raw)
    
    # Check if the query returned any results
    # if(type(query)==int and query == 0):
    if query is None:
        plots_logger.debug(f"No results found for {station_name} and {analyte_name}")
        return None
    
    else:
        # Filter the data by the specified column and values (if necessary)
        if(filter):
            filter_res = filters.filter_by_column(
                data=data_pylenm_dm.construction_data, 
                col=col, 
                equals=equals,
            )
            
            # if('ERROR:' in str(filter_res)):
            if filter_res is None:
                plots_logger.debug("Ran into ERROR when calling filter_by_column()!")
                return filter_res

            # Get the intersection of the query stations and filter stations
            query_stations = list(query.STATION_ID.unique())
            filter_stations = list(filter_res.index.unique())
            intersect_stations = list(set(query_stations) & set(filter_stations))
            if(len(intersect_stations)<=0):
                plots_logger.debug('ERROR: No results for this query with the specifed filter parameters.')
                # return 'ERROR: No results for this query with the specifed filter parameters.'
                return None
            query = query[query['STATION_ID'].isin(intersect_stations)]
        
        # Extract the date and result values from the query data
        x_data = query.COLLECTION_DATE
        x_data = pd.to_datetime(x_data)
        y_data = query.RESULT
        if(log_transform):
            y_data = np.log10(y_data)

        # Prepare data for Supersmoother
        x_RR = x_data.astype(int).to_numpy()
        nu = x_data.shape[0]

        # Create and fit supersmoother model
        model = SuperSmoother(alpha=alpha)
        try:
            model.fit(x_RR, y_data)

            # Get smoothed predictions, residuals and outliers from the model.
            y_pred = model.predict(x_RR)
        
        except Exception as e:
            plots_logger.debug(e)
            return None

        r = model.cv_residuals()
        out = abs(r) > (2.2 * np.std(r))
        out_x = x_data[out]
            
        out_y = y_data[out]

        # Plot data
        plt.figure(figsize=(8,8))
        ax = plt.axes()
        years = mdates.YearLocator(year_interval)  # every year
        months = mdates.MonthLocator()  # every month
        yearsFmt = mdates.DateFormatter('%Y')

        for label in ax.get_xticklabels():
            label.set_rotation(30)
            label.set_horizontalalignment('center')

        ax = plt.gca()
        ax.xaxis.set_major_locator(years)
        ax.xaxis.set_major_formatter(yearsFmt)
        ax.autoscale_view()

        unit = query.RESULT_UNITS.values[0]

        ax.set_title(
            # str(station_name) + ' - ' + analyte_name, 
            f"{station_name} - {analyte_name}", fontweight='bold',
        )
        ttl = ax.title
        ttl.set_position([.5, 1.05])
        
        if(y_label==''):    
            if(log_transform):
                ax.set_ylabel('log-Concentration (' + unit + ')')
            else:
                ax.set_ylabel('Concentration (' + unit + ')')
        
        else:
            ax.set_ylabel(y_label)
        
        ax.set_xlabel(x_label)
        small_fontSize = 15
        large_fontSize = 20
        
        plt.rc('axes', titlesize=large_fontSize)
        plt.rc('axes', labelsize=large_fontSize)
        plt.rc('legend', fontsize=small_fontSize)
        plt.rc('xtick', labelsize=small_fontSize)
        plt.rc('ytick', labelsize=small_fontSize) 
        
        ax.plot(x_data, y_data, ls='', marker='o', ms=5, color='black', alpha=1)
        ax.plot(x_data, y_pred, ls='-', marker='', ms=5, lw=2, color='blue', alpha=0.5, label="Super Smoother")
        ax.plot(out_x , out_y, ls='', marker='o', ms=5, color='red', alpha=1, label="Outliers")
        ax.legend(bbox_to_anchor=(1.04, 1), loc='upper left', borderaxespad=0.)
        
        props = dict(boxstyle='round', facecolor='grey', alpha=0.15)       
        ax.text(1.05, 0.85, 'Samples: {}'.format(nu), transform=ax.transAxes, 
                fontsize=small_fontSize,
                fontweight='bold',
                verticalalignment='top', 
                bbox=props)
        
        # Save the plot to `save_dir` (if necessary)
        if save_plot:
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            plt.savefig(
                f"{save_dir}/{str(station_name)}-{analyte_name}.png", 
                bbox_inches="tight",
            )
        
        # Show plot inline (used in notebooks)
        if(plot_inline):
            plt.show()
        
        # Clear all the plots from memory
        plt.clf()
        plt.cla()
        plt.close()

        # Returning a non-None value to mark a successful execution
        return 200


def plot_all_data(
        # self, 
        data_pylenm_dm,
        log_transform=True, 
        alpha=0, 
        year_interval=2, 
        plot_inline=True, 
        save_dir='plot_data',
    ):
    """Plot concentrations over time for every station and analyte with a smoothed curve on interpolated data points.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        log_transform (bool, optional): choose whether or not the data should be transformed to log base 10 values. Defaults to True.
        alpha (int, optional): alue between 0 and 10 for line smoothing. Defaults to 0.
        plot_inline (bool, optional): choose whether or not to show plot inline. Defaults to True.
        year_interval (int, optional): plot by how many years to appear in the axis e.g.(1 = every year, 5 = every 5 years, ...). Defaults to 2.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_data'.
    """
    
    analytes = ['TRITIUM','URANIUM-238','IODINE-129','SPECIFIC CONDUCTANCE', 'PH', 'DEPTH_TO_WATER']
    
    stations = np.array(data_pylenm_dm.data.STATION_ID.values)
    stations = np.unique(stations)

    stations_analytes_iter = list(itertools.product(stations, analytes))
    
    success = 0
    errors = 0
    
    # for station in stations:
    #     for analyte in analytes:
    for station, analyte in tqdm(stations_analytes_iter, desc="Station Analyte Pair", total=(len(analytes)*len(stations))):
        plot = plot_data(
            data_pylenm_dm=data_pylenm_dm,
            station_name=station, 
            analyte_name=analyte, 
            log_transform=log_transform, 
            alpha=alpha, 
            year_interval=year_interval, 
            plot_inline=plot_inline, 
            save_plot=True,
            save_dir=save_dir,
        )
        
        if plot is None:
            errors = errors + 1
        else:
            success = success + 1

    plots_logger.info(f"Success: {success}")
    plots_logger.info(f"Errors: {errors}")

        
def plot_MCL(
        # self, 
        data_pylenm_dm,
        station_name, 
        analyte_name, 
        year_interval=5, 
        save_dir='plot_MCL',
    ):
    """Plots the linear regression line of data given the analyte_name and station_name. The plot includes the prediction where the line of best fit intersects with the Maximum Concentration Limit (MCL).

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        station_name (str): ame of the station to be processed
        analyte_name (str): name of the analyte to be processed
        year_interval (int, optional): lot by how many years to appear in the axis e.g.(1 = every year, 5 = every 5 years, ...). Defaults to 5.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_MCL'.
    """
    
    # data = self.data
    data = data_pylenm_dm.data
    
    # finds the intersection point of 2 lines given the slopes and y-intercepts
    def line_intersect(m1, b1, m2, b2):
        if m1 == m2:
            print ('The lines are parallel')
            return None
        x = (b2 - b1) / (m1 - m2)
        y = m1 * x + b1
        return x,y

    # Gets appropriate data (station_name and analyte_name)
    # query = self.query_data(station_name, analyte_name)
    query = filters.query_data(
        data_pylenm_dm=data_pylenm_dm,
        station_name=station_name, 
        analyte_name=analyte_name,
    )

    # if(type(query)==int and query == 0):
    if isinstance(query, int) and query==0:
        plots_logger.info(
            f"No results found for {station_name} and {analyte_name}"
        )
        return None
    
    else:
        test = query.groupby(['COLLECTION_DATE'])['RESULT'].mean()
        test.index = pd.to_datetime(test.index)

        x = date2num(test.index)
        # y = np.log10(test.RESULT)
        y = np.log10(test)
        ylabel = 'log-Concentration (' + data_pylenm_dm.get_unit(analyte_name) + ')'
        y = y.rename(ylabel)

        p, cov = np.polyfit(x, y, 1, cov=True)  # parameters and covariance from of the fit of 1-D polynom.

        m_unc = np.sqrt(cov[0][0])
        b_unc = np.sqrt(cov[1][1])

        f = np.poly1d(p)

        try:
            # MCL = self.get_MCL(analyte_name)
            MCL = preprocess.get_MCL(analyte_name)
            m1, b1 = f # line of best fit
            m2, b2 = 0, MCL # MCL constant

            intersection = line_intersect(m1, b1, m2, b2)

            ## Get confidence interval intersection points with MCL
            data = list(zip(x,y))
            n = len(data)
            list_slopes = []
            list_intercepts = []
            random.seed(50)
            
            for _ in range(80):
                sampled_data = [ random.choice(data) for _ in range(n) ]
                x_s, y_s = zip(*sampled_data)
                x_s = np.array(x_s)
                y_s = np.array(y_s)

                m_s, b_s, r, p, err = scipy.stats.linregress(x_s,y_s)
                ymodel = (m_s * x_s) + b_s
                
                list_slopes.append(m_s)
                list_intercepts.append(b_s)

            max_index = list_slopes.index(max(list_slopes))
            min_index = list_slopes.index(min(list_slopes))
            
            intersection_left = line_intersect(list_slopes[min_index], list_intercepts[min_index], m2, b2)
            intersection_right = line_intersect(list_slopes[max_index], list_intercepts[max_index], m2, b2)
            ##

            fig, ax = plt.subplots(figsize=(10, 6))

            ax.set_title(station_name + ' - ' + analyte_name, fontweight='bold')
            ttl = ax.title
            ttl.set_position([.5, 1.05])
            years = mdates.YearLocator(year_interval)  # every year
            months = mdates.MonthLocator()  # every month
            yearsFmt = mdates.DateFormatter('%Y') 

            for label in ax.get_xticklabels():
                label.set_rotation(30)
                label.set_horizontalalignment('center')

            ax.xaxis.set_major_locator(years)
            ax = plt.gca()
            ax.xaxis.set_major_locator(years)
            ax.xaxis.set_major_formatter(yearsFmt)
            ax.autoscale_view()
            ax.grid(True, alpha=0.4)
            small_fontSize = 15
            large_fontSize = 20
            
            plt.rc('axes', titlesize=large_fontSize)
            plt.rc('axes', labelsize=large_fontSize)
            plt.rc('legend', fontsize=small_fontSize)
            plt.rc('xtick', labelsize=small_fontSize)
            plt.rc('ytick', labelsize=small_fontSize)

            ax.set_xlabel('Years')
            ax.set_ylabel('log-Concentration (' + data_pylenm_dm.get_unit(analyte_name) + ')')

            if(intersection[0] < min(x)):
                temp = intersection_left
                intersection_left = intersection_right
                intersection_right = temp
                ax.set_ylim([0, max(y)+1])
                ax.set_xlim([intersection_left[0]-1000, max(x)+1000])
            
            elif(intersection[0] < max(x) and intersection[0] > min(x)):
                ax.set_ylim([0, max(y)+1])
                ax.set_xlim(min(x)-1000, max(x)+1000)

            else:
                ax.set_ylim([0, max(y)+1])
                ax.set_xlim([min(x)-1000, intersection_right[0]+1000])

            ax = sns.regplot(
                x=x, y=y, logx=True, truncate=False, 
                seed=42, n_boot=1000, ci=95,
            ) # Line of best fit
            ax.plot(x, y, ls='', marker='o', ms=5, color='black', alpha=1) # Data
            ax.axhline(y=MCL, color='r', linestyle='--') # MCL
            ax.plot(intersection[0], intersection[1], color='blue', marker='o', ms=10)
            ax.plot(intersection_left[0], intersection_left[1], color='green', marker='o', ms=5)
            ax.plot(intersection_right[0], intersection_right[1], color='green', marker='o', ms=5)

            predict = num2date(intersection[0]).date()
            l_predict = num2date(intersection_left[0]).date()
            u_predict = num2date(intersection_right[0]).date()
            
            ax.annotate(
                predict, 
                (intersection[0], intersection[1]), 
                xytext=(intersection[0], intersection[1]+1), 
                bbox=dict(boxstyle="round", alpha=0.1), 
                ha='center', 
                arrowprops=dict(arrowstyle="->", color='blue'), 
                fontsize=small_fontSize, 
                fontweight='bold',
            )
            
            props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
            ax.text(
                1.1, 
                0.5, 
                # 'Lower confidence: {}\nPrediction: {}\nUpper confidence: {}'.format(l_predict, predict, u_predict), 
                f"Lower confidence: {l_predict}\nPrediction: {predict}\nUpper confidence: {u_predict}", 
                transform=ax.transAxes, 
                fontsize=small_fontSize, 
                fontweight='bold', 
                verticalalignment='bottom', 
                bbox=props,
            )

            # Save figure
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            
            plt.savefig(
                # save_dir + '/' + station_name + '-' + analyte_name +'.png', 
                f"{save_dir}/{station_name}-{analyte_name}.png", 
                bbox_inches="tight",
            )

        except Exception as e:
            plots_logger.error('ERROR: Something went wrong')
            plots_logger.error(e)
            return None



## by K. Whiteaker 2024-08, kwhit@alum.mit.edu
# Helper function for plot_data_weekAvg() and plot_data_lowess()
# copied from self.plot_data(), with a few additional plot items added
def _plot_data_xOutliers_fit(
        # self, 
        station_name, 
        analyte_name, 
        concentration_data_xOutliers, 
        concentration_data_xOutliers_fit, 
        outliers, 
        units,
        difference, 
        fitColor, 
        fitName, 
        rm_outliers, 
        std_thresh, 
        lowess_frac, 
        show_difference, 
        x_label, 
        y_label, 
        year_interval, 
        log_transform, 
        y_zoom, 
        return_data, 
        save, 
        save_dir, 
        plot_inline, 
        save_as_pdf,
    ):
    
    # define figure shape and size
    plt.figure(figsize=(8,8))
    ax = plt.axes()
    
    years = mdates.YearLocator(year_interval)  # every year
    months = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y')
    
    for label in ax.get_xticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment('center')
    
    ax = plt.gca()
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(yearsFmt)
    ax.autoscale_view()
    
    if y_label is None:
        y_label = analyte_name + f' [{units}]'
        if(log_transform):
            y_label = 'log_10 of ' + analyte_name + f' [{units}]'
    
    ax.set_ylabel(y_label)
    ax.set_xlabel(x_label)
    
    small_fontSize = 15
    large_fontSize = 20
    plt.rc('axes', titlesize=large_fontSize)
    plt.rc('axes', labelsize=large_fontSize)
    plt.rc('legend', fontsize=small_fontSize)
    plt.rc('xtick', labelsize=small_fontSize)
    plt.rc('ytick', labelsize=small_fontSize)
    
    props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
    ax.text(
        1.05, 
        0.85, 
        f"Samples: {concentration_data_xOutliers.size}", 
        transform=ax.transAxes, 
        fontsize=small_fontSize,
        fontweight='bold',
        verticalalignment='top', 
        bbox=props,
    )
    
    # ax.set_title(str(station_name) + ' - ' + analyte_name, fontweight='bold')
    ax.set_title(f"{station_name} - {analyte_name}", fontweight='bold')
    ttl = ax.title
    ttl.set_position([.5, 1.05])

    # plot items
    ax.plot(
        concentration_data_xOutliers.index, 
        concentration_data_xOutliers, 
        ls='', 
        marker='o', 
        ms=2, 
        color='black', 
        alpha=1,
    )
    
    if rm_outliers:
        ax.plot(
            outliers.index, 
            outliers, 
            ls='', 
            marker='x', 
            ms=5, 
            color='red', 
            alpha=0.75, 
            label=f"Outliers ({outliers.size})",
        )
    
    # need to use pd.Series.dropna() when plotting so that matplotlib doesn't insert any gaps
    ax.plot(
        concentration_data_xOutliers_fit.dropna().index, 
        concentration_data_xOutliers_fit.dropna(), 
        ls='-', 
        marker='', 
        ms=2, 
        lw=1, 
        color=fitColor, 
        alpha=1, 
        label=fitName,
    )
    
    if y_zoom:
        plt.ylim(min(concentration_data_xOutliers), max(concentration_data_xOutliers))
    ax.legend(bbox_to_anchor=(1.04, 1), loc='upper left', borderaxespad=0.)
    
    fileType = '.pdf' if save_as_pdf else '.png'
    
    if save:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(save_dir + '/' + str(station_name) + '_' + analyte_name + '_' + fitName.replace(" ", "") + fileType, bbox_inches="tight")
    
    if(plot_inline):
        plt.show()
    
    plt.clf()
    plt.cla()
    plt.close()
    
    # plot difference between fit and observation, with outliers marked
    if rm_outliers and show_difference:
        outliers_loc = outliers.index
        
        plt.scatter(difference.index, difference, s=0.5, label='Difference')
        
        plt.scatter(difference[outliers_loc].index, difference[outliers_loc], s=10, marker='x', label=f'Outliers outside {std_thresh} SD')
        
        plt.title("(obs - fit), fraction "+ str(lowess_frac) + ", station "+ str(station_name))
        
        plt.ylabel('Observation - LOWESS Fit [ft]')
        plt.legend()
        
        if save:
            # plt.savefig(save_dir + '/' + str(station_name) + '_' + analyte_name + '_difference' + fileType, bbox_inches="tight")
            plt.savefig(
                f"{save_dir}/{station_name}_{analyte_name}_difference{fileType}", 
                bbox_inches="tight",
            )
        
        if(plot_inline):
            plt.show()
        
        plt.clf()
        plt.cla()
        plt.close()


## by K. Whiteaker 2024-08, kwhit@alum.mit.edu
def plot_data_rollAvg(
        # self, 
        data_pylenm_dm,
        station_name, 
        analyte_name, 
        window='1W', 
        rm_outliers=True, 
        std_thresh=2.2, 
        lowess_frac=0.1, 
        show_difference=False, 
        x_label=None, 
        y_label=None, 
        year_interval=2, 
        log_transform=False, 
        y_zoom=False, 
        return_data=False, 
        save=False, 
        save_dir='plot_data_lowess',
        plot_inline=True, 
        save_as_pdf=False,
    ):
    """Plot time series data for a specified station and analyte, alongside a curve calculated via a rolling time average (default 1 week) centred around each point

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        station_name (str): name of the station to be processed
        analyte_name (str): name of the analyte to be processed
        window (str): the time length of the rolling window used for averaging, ex. '1W', '2D'. Defaults to '1W'. Months are not supported; valid inputs are listed here: https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html
        rm_outliers (bool, optional): if True, use lowess function to remove outliers. Defaults to True
        std_thresh (int, optional): number of standard deviations in (observation - LOWESS fit) outside of which is considered an outlier. Defaults to 2.2
        lowess_frac (float, optional): fraction of total data points considered in local fits of LOWESS smoother. A smaller value captures more local behaviour, and may be required for large datasets. Defaults to 0.1
        show_difference (bool, optional): if True, plots (observation - LOWESS fit) at each measurement time
        x_label (str, optional): x axis label. Defaults to None
        y_label (str, optional): y axis label. Defaults to None. If left as None, will be set to analyte_name + f' [{units}]'
        year_interval (int, optional): plot by how many years to appear in the axis e.g.(1 = every year, 5 = every 5 years, ...). Defaults to 2
        log_transform (bool, optional): choose whether or not the data should be transformed to log base 10 values. Defaults to False
        y_zoom (bool, optional): if True, plot y axes will zoom in to [minimum y value, maximum y value] after removing outliers. Defaults to False
        return_data (bool, optional): if True, return the data used to make the plots. Defaults to False
        save (bool, optional): if True, save plots to file in save_dir. Defaults to False
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_data_lowess'
        plot_inline (bool, optional): choose whether or not to show plot inline. Defaults to True
        save_as_pdf (bool, optional): if True, saves figure as vectorized pdf instead of png. Defaults to False

    Returns:
        if return_data:
            if rm_outliers:
                pd.Series: concentration data with outliers set to np.nan, indexed by date
                pd.Series: rolling average of post-outlier-removal concentration data, indexed by date
                pd.Series: outliers, indexed by date
            else:
                pd.Series: concentration data, indexed by date
                pd.Series: rolling average of concentration data, indexed by date
        else:
            None
    """

    # import data and filter for the chosen analyte
    query = filters.query_data(
        data_pylenm_dm=data_pylenm_dm,
        station_name=station_name, 
        analyte_name=analyte_name,
    )
    query = filters.simplify_data(data=query)
    
    # if(type(query)==int and query == 0):
    if query is None:
        # return 'No results found for {} and {}'.format(station_name, analyte_name)
        plots_logger.debug(f"No results found for {station_name} and {analyte_name}.")
        return None
    
    # reshape data so it can be passed into lowess() and remove_outliers_lowess()
    concentration_data = query.set_index('COLLECTION_DATE').RESULT
    time_data = concentration_data.index
    units = query.RESULT_UNITS.values[0]  # assumes same units for every entry
    
    if(log_transform):
        concentration_data = concentration_data.apply(np.log10)

    # if rm_outliers toggled, call remove_outliers_lowess() to set outliers to np.nan via lowess fit
    if rm_outliers:
        concentration_data_xOutliers, difference = preprocess.remove_outliers_lowess(
            concentration_data, 
            lowess_frac=lowess_frac,
            std_thresh=std_thresh, 
            return_difference=True,
        )
    else:
        concentration_data_xOutliers = concentration_data
        difference = pd.Series  # will not be used, just need something to pass into _plot_data_xOutliers_fit
    
    outliers_loc = concentration_data_xOutliers.compare(concentration_data).index  # dates at which outliers occurred
    outliers = concentration_data[outliers_loc]  # outlier values indexed by their date of occurrence
    num_outliers = outliers.size
    analyte_max_xOutliers, analyte_min_xOutliers = max(concentration_data_xOutliers), min(concentration_data_xOutliers)
            
    # now that outliers have been dropped, calculate rolling average of y values (result of average is assigned to the centre of the window)
    avg_window = pd.Timedelta(window)  # could also do pd.Timedelta(7, "d")
    rolling_ser = concentration_data_xOutliers
    roll = rolling_ser.rolling(
        window=avg_window, 
        min_periods=1, 
        center=True,
    )
    concentration_data_xOutliers_avg = roll.mean()
    
    # plot rolling average
    _plot_data_xOutliers_fit(
        station_name=station_name, 
        analyte_name=analyte_name, 
        concentration_data_xOutliers=concentration_data_xOutliers, 
        concentration_data_xOutliers_fit=concentration_data_xOutliers_avg, 
        outliers=outliers, 
        units=units, 
        difference=difference, 
        fitColor='green', 
        fitName=f"{window} Average", 
        rm_outliers=rm_outliers, 
        std_thresh=std_thresh, 
        lowess_frac=lowess_frac, 
        show_difference=show_difference, 
        x_label=x_label, 
        y_label=y_label, 
        year_interval=year_interval, 
        log_transform=log_transform, 
        y_zoom=y_zoom, 
        return_data=return_data, 
        save=save, 
        save_dir=save_dir, 
        plot_inline=plot_inline, 
        save_as_pdf=save_as_pdf,
    )
    
    # return data
    if return_data:
        if rm_outliers:
            return concentration_data_xOutliers, concentration_data_xOutliers_avg, outliers
        else:
            return concentration_data_xOutliers, concentration_data_xOutliers_avg
        

## by K. Whiteaker 2024-08, kwhit@alum.mit.edu
def plot_data_lowess(
        # self, 
        data_pylenm_dm, 
        station_name, 
        analyte_name, 
        rm_outliers=True, 
        std_thresh=2.2, 
        lowess_frac=0.1, 
        show_difference=False, 
        x_label=None, 
        y_label=None, 
        year_interval=2, 
        log_transform=False, 
        y_zoom=False, 
        return_data=False, 
        save=False, 
        save_dir='plot_data_lowess', 
        plot_inline=True, 
        save_as_pdf=False,
    ):
    """Plot time series data for a specified station and analyte, alongside a 
    smoothed curve on interpolated data points (LOWESS).
    
    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object 
            containing the concentration and construction data.
        station_name (str): name of the station to be processed
        analyte_name (str): name of the analyte to be processed
        rm_outliers (bool, optional): if True, use lowess function to 
            remove outliers. Defaults to True
        std_thresh (int, optional): number of standard deviations in 
            (observation - fit) outside of which is considered an outlier. 
            Defaults to 2.2
        lowess_frac (float, optional): fraction of total data points 
            considered in local fits of LOWESS smoother. A smaller value 
            captures more local behaviour, and may be required for large 
            datasets. Defaults to 0.1
        show_difference (bool, optional): if True, plots (observation - 
            LOWESS fit) at each measurement time
        x_label (str, optional): x axis label. Defaults to None
        y_label (str, optional): y axis label. Defaults to None. If left as 
            None, will be set to analyte_name + f' [{units}]'
        year_interval (int, optional): plot by how many years to appear in 
            the axis e.g.(1 = every year, 5 = every 5 years, ...). 
            Defaults to 2.
        log_transform (bool, optional): choose whether or not the data 
            should be transformed to log base 10 values. Defaults to False
        y_zoom (bool, optional): if True, plot y axes will zoom in to 
            [minimum y value, maximum y value] after removing outliers. 
            Defaults to False.
        return_data (bool, optional): if True, return the data used to make 
            the plots. Defaults to False.
        save (bool, optional): if True, save plots to file in save_dir. 
            Defaults to False.
        save_dir (str, optional): name of the directory you want to save 
            the plot to. Defaults to 'plot_data_lowess'.
        plot_inline (bool, optional): choose whether or not to show plot 
            inline. Defaults to True.
        save_as_pdf (bool, optional): if True, saves figure as vectorized 
            pdf instead of png. Defaults to False.

    Returns:
        if return_data:
            if rm_outliers:
                pd.Series: concentration data with outliers set to np.nan, indexed by date
                pd.Series: LOWESS fit of post-outlier-removal concentration data, indexed by date
                pd.Series: outliers, indexed by date
            else:
                pd.Series: concentration data, indexed by date
                pd.Series: LOWESS fit of concentration data, indexed by date
        else:
            None
    """
    
    # import data and filter for the chosen analyte
    query = filters.query_data(
        data_pylenm_dm=data_pylenm_dm, 
        station_name=station_name, 
        analyte_name=analyte_name,
    )
    query = filters.simplify_data(data=query)
    
    # if(type(query)==int and query == 0):
    if query is None:
        # return 'No results found for {} and {}'.format(station_name, analyte_name)
        plots_logger.debug(f"No results found for {station_name} and {analyte_name}.")
        return None
    
    # reshape data so it can be passed into lowess() and remove_outliers_lowess()
    concentration_data = query.set_index('COLLECTION_DATE').RESULT
    time_data = concentration_data.index
    units = query.RESULT_UNITS.values[0]  # assumes same units for every entry
    
    if(log_transform):
        concentration_data = concentration_data.apply(np.log10)

    # if rm_outliers toggled, call remove_outliers_lowess() to set outliers to np.nan via lowess fit
    if rm_outliers:
        concentration_data_xOutliers, difference = preprocess.remove_outliers_lowess(
            concentration_data, 
            lowess_frac=lowess_frac, 
            std_thresh=std_thresh, 
            return_difference=True,
        )
    else:
        concentration_data_xOutliers = concentration_data
        difference = pd.Series  # will not be used, just need something to pass into __plot_data_xOutliers_fit
    
    outliers_loc = concentration_data_xOutliers.compare(concentration_data).index  # dates at which outliers occurred
    outliers = concentration_data[outliers_loc]  # outlier values indexed by their date of occurrence
    num_outliers = outliers.size
    analyte_max_xOutliers, analyte_min_xOutliers = max(concentration_data_xOutliers), min(concentration_data_xOutliers)
            
    # now that outliers have been dropped, calculate a new lowess fit for plotting
    # lowess() reads datetime as nanoseconds and has issues with large numbers & low frac, so need to scale down x-axis during fitting
    scaleDown = 1e17
    x_readable = time_data.to_numpy().astype(int)/scaleDown
    y_data_lowess = lowess(
        concentration_data_xOutliers, 
        x_readable, 
        frac=lowess_frac, 
        return_sorted=False,
    )
    concentration_data_xOutliers_lowess = pd.Series(
        data=y_data_lowess, 
        index=concentration_data_xOutliers.index, 
        dtype=float, 
        name=concentration_data_xOutliers.name,
    )

    # plot lowess fit
    _plot_data_xOutliers_fit(
        station_name=station_name, 
        analyte_name=analyte_name, 
        concentration_data_xOutliers=concentration_data_xOutliers, 
        concentration_data_xOutliers_fit=concentration_data_xOutliers_lowess, 
        outliers=outliers, 
        units=units, 
        difference=difference, 
        fitColor='violet', 
        fitName='LOWESS Fit', 
        rm_outliers=rm_outliers, 
        std_thresh=std_thresh,
        lowess_frac=lowess_frac, 
        show_difference=show_difference, 
        x_label=x_label, 
        y_label=y_label,
        year_interval=year_interval, 
        log_transform=log_transform, 
        y_zoom=y_zoom, 
        return_data=return_data,
        save=save, 
        save_dir=save_dir, 
        plot_inline=plot_inline, 
        save_as_pdf=save_as_pdf,
    )
    
    if return_data:
        if rm_outliers:
            return concentration_data_xOutliers, concentration_data_xOutliers_lowess, outliers
        else:
            return concentration_data_xOutliers, concentration_data_xOutliers_lowess
