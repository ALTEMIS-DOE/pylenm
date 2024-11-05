import os
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

import pylenm2
from pylenm2.data import filters
from pylenm2.data import fetchers
from pylenm2.stats import preprocess
from pylenm2.data import transformation
from pylenm2.visualization import plots
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

plot_corr_logger = logger_config.setup_logging(
    module_name=__name__,
    level=logging.INFO,
    # level=logging.DEBUG,
    logfile_dir=c.LOGFILE_DIR,
)


def _format_ticks_and_labels(axes, fontsize):
    """Utility function to format all the ticks and labels along the `x` and `y` axes.

    Args:
        axes (matplotlib.Axes): List of Axes to format.
        fontsize (int): Fontsize for all the labels and ticks.
    """
    # for ax in g.axes.flat:
    for ax in axes.flat:
        ax.tick_params("y", labelrotation=0, labelsize=fontsize)
        
        # ax.set_xticklabels(
        #     ax.get_xticklabels(), rotation=45, fontsize=fontsize,
        # )   # NOTE: DISCOURAGED version of `xticks` - https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xticklabels.html
        # ax.set_xticks(
        #     ax.get_xticks(),
        #     [f"{num:.2f}" for num in ax.get_xticks()],    # TODO: Check with Zexuan if this is supposed to be normalized.
        #     rotation=45, 
        #     fontsize=fontsize,
        # )     # NOTE: Non-discouraged version of `xticks`
        ax.tick_params("x", labelrotation=45, labelsize=fontsize)   # NOTE: Cleaner version of `xticks`. TODO: Check if this one works or if we need to revert back to the previous ticklabels. @Zexuan
        ax.set_xlabel(
            ax.get_xlabel(), fontsize=fontsize, fontweight='bold',
        )
        ax.set_ylabel(
            ax.get_ylabel(), fontsize=fontsize, fontweight='bold',
        )


def plot_corr_by_station(
        data_pylenm_dm, 
        station_name, 
        analytes, 
        remove_outliers=True, 
        z_threshold=4, 
        interpolate=False, 
        frequency='2W', 
        save_dir='plot_correlation', 
        log_transform=False, 
        fontsize=20, 
        # returnData=False, 
        remove=[], 
        no_log=None,
    ):
    """Plots the correlations with the physical plots as station as the correlations of the important analytes over time for a specified station.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        station_name (str): name of the station to be processed
        analytes (list): list of analyte names to use
        remove_outliers (bool, optional): choose whether or to remove the outliers. Defaults to True.
        z_threshold (int, optional): z_score threshold to eliminate outliers. Defaults to 4.
        interpolate (bool, optional): choose whether or to interpolate the data. Defaults to False.
        frequency (str, optional): {‘D’, ‘W’, ‘M’, ‘Y’} frequency to interpolate. Note: See https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html for valid frequency inputs. (e.g. ‘W’ = every week, ‘D ’= every day, ‘2W’ = every 2 weeks). Defaults to '2W'.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_correlation'.
        log_transform (bool, optional): flag for log base 10 transformation. Defaults to False.
        fontsize (int, optional): font size. Defaults to 20.
        returnData (bool, optional): flag to return data used to perfrom correlation analysis. Defaults to False.
        remove (list, optional): stations to remove. Defaults to [].
        no_log (list, optional): list of column names to not apply log transformation to. Defaults to None.

    Returns:
        None
    """
    
    # Prepare data
    data = data_pylenm_dm.data
    query = data[data.STATION_ID == station_name]
    a = list(np.unique(query.ANALYTE_NAME.values))# get all analytes from dataset
    
    for value in analytes:
        # if((value in a)==False):
        if value not in a:
            plot_corr_logger.error(f"ERROR: No analyte named {value} in data.")
            # return 'ERROR: No analyte named "{}" in data.'.format(value)
            return None

    analytes = sorted(analytes)
    
    query = query.loc[query.ANALYTE_NAME.isin(analytes)]
    x = query[['COLLECTION_DATE', 'ANALYTE_NAME']]
    unique = ~x.duplicated()
    query = query[unique]
    piv = query.reset_index().pivot(index='COLLECTION_DATE',columns='ANALYTE_NAME', values='RESULT')
    piv = piv[analytes]
    piv.index = pd.to_datetime(piv.index)
    totalSamples = piv.shape[0]
    piv = piv.dropna()
    
    # Interpolate (if specified) and define file extension and titles
    if(interpolate):
        piv = transformation.interpolate_station_data(
            data_pylenm_dm=data_pylenm_dm,
            station_name=station_name, 
            analytes=analytes, 
            frequency=frequency,
        )
        # file_extension = '_interpolated_' + frequency
        file_extension = f"_interpolated_{frequency}"
        # title = station_name + '_correlation - interpolated every ' + frequency
        title = f"{station_name}_correlation - interpolated every {frequency}"
    else:
        file_extension = "_correlation"
        # title = station_name + '_correlation'
        title = f"{station_name}_correlation"
    samples = piv.shape[0]
    
    # Plot if there are enough samples
    if(samples < 5):
        if(interpolate):
            # return 'ERROR: {} does not have enough samples to plot.\n Try a different interpolation frequency'.format(station_name)
            plot_corr_logger.error(f"ERROR: {station_name} does not have enough samples to plot.\n Try a different interpolation frequency")
            return None
        
        # return 'ERROR: {} does not have enough samples to plot.'.format(station_name)
        plot_corr_logger.error(f"ERROR: {station_name} does not have enough samples to plot.")
        return None

    else:
        # scaler = StandardScaler()
        # pivScaled = scaler.fit_transform(piv)
        # pivScaled = pd.DataFrame(pivScaled, columns=piv.columns)
        # pivScaled.index = piv.index
        # piv = pivScaled
        
        # Perform log transformation (if specified)
        if(log_transform):
            piv[piv <= 0] = 0.00000001
            temp = piv.copy()
            piv = np.log10(piv)
            if(no_log !=None):
                for col in no_log:
                    piv[col] = temp[col]

        # Remove outliers
        if(remove_outliers):
            piv = preprocess.remove_outliers(piv, z_threshold=z_threshold)

        # Prepare data to plot
        samples = piv.shape[0]

        idx = piv.index.date
        dates = [dates.strftime('%Y-%m-%d') for dates in idx]
        remaining = [i for i in dates if i not in remove]
        piv = piv.loc[remaining]
        
        # Create the seaborn PairGrid plot
        sns.set_style("white", {"axes.facecolor": "0.95"})
        g = sns.PairGrid(piv, aspect=1.2, diag_sharey=False, despine=False)
        g.figure.suptitle(title, fontweight='bold', y=1.08, fontsize=25)
        g.map_lower(
            sns.regplot, 
            lowess=True, 
            ci=False, 
            line_kws={'color': 'red', 'lw': 3}, 
            scatter_kws={'color': 'black', 's': 20},
        )
        # g.map_diag(
        #     sns.distplot, 
        #     kde_kws={'color': 'black', 'lw': 3}, 
        #     hist_kws={'histtype': 'bar', 'lw': 2, 'edgecolor': 'k', 'facecolor':'grey'},
        # )   # NOTE: Replacing with `sns.histplot` because `sns.distplot` is deprecated.
        g.map_diag(
            sns.histplot, 
            stat="density",
            color="black",
            element="bars",
            edgecolor="dimgray",
            facecolor="lightgray",
            lw=2,
            kde=True, 
            kde_kws=dict(cut=3),
            line_kws=dict(color="black", lw=3),
        )

        g.map_upper(plots._plotUpperHalf)

        # for ax in g.axes.flat:
        #     ax.tick_params("y", labelrotation=0, labelsize=fontsize)
            
        #     # ax.set_xticklabels(
        #     #     ax.get_xticklabels(), rotation=45, fontsize=fontsize,
        #     # )   # NOTE: DISCOURAGED version of `xticks` - https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xticklabels.html
        #     # ax.set_xticks(
        #     #     ax.get_xticks(),
        #     #     [f"{num:.2f}" for num in ax.get_xticks()],    # TODO: Check with Zexuan if this is supposed to be normalized.
        #     #     rotation=45, 
        #     #     fontsize=fontsize,
        #     # )     # NOTE: Non-discouraged version of `xticks`
        #     ax.tick_params("x", labelrotation=45, labelsize=fontsize)   # NOTE: Cleaner version of `xticks`
        #     ax.set_xlabel(
        #         ax.get_xlabel(), fontsize=fontsize, fontweight='bold',
        #     )
        #     ax.set_ylabel(
        #         ax.get_ylabel(), fontsize=fontsize, fontweight='bold',
        #     )
        _format_ticks_and_labels(axes=g.axes, fontsize=fontsize)
            
        # g.figure.subplots_adjust(wspace=0.3, hspace=0.3)
        g.figure.subplots_adjust(
            wspace=0.3, hspace=0.3, 
            left=0.1, right=0.75, 
            bottom=0.1, top=0.9,
        )

        # Add the meta-data information box
        ax = plt.gca()
        props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
        ax.text(
            1.3, 6.2, 
            # 'Start date:  {}\nEnd date:    {}\n\nOriginal samples:     {}\nSamples used:     {}'.format(piv.index[0].date(), piv.index[-1].date(), totalSamples, samples), 
            f"{'Start date:':<12} {str(piv.index[0].date())}\n{'End date:':<12} {str(piv.index[-1].date())}\n\n{'Original samples:':<20} {totalSamples}\n{'Samples used:':<20} {samples}",
            transform=ax.transAxes, 
            fontsize=fontsize, 
            fontweight='bold', 
            verticalalignment='bottom', 
            bbox=props,
        )
        
        # Add titles to the diagonal axes/subplots
        for ax, col in zip(np.diag(g.axes), piv.columns):
            ax.set_title(col, y=0.82, fontsize=15)
        
        # Save figure (if specified)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        g.figure.savefig(
            # save_dir + '/' + station_name + file_extension + '.png', 
            f"{save_dir}/{station_name}{file_extension}.png", 
            bbox_inches="tight",
        )
        
        # if(returnData):
        #     return piv
    
        return piv
    

def plot_all_corr_by_station(
        # self, 
        data_pylenm_dm,
        analytes, 
        remove_outliers=True, 
        z_threshold=4, 
        interpolate=False, 
        frequency='2W', 
        save_dir='plot_correlation', 
        log_transform=False, 
        fontsize=20,
    ):
    """Plots the correlations with the physical plots as station as the important analytes over time for each station in the dataset.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        analytes (list): list of analyte names to use
        remove_outliers (bool, optional): choose whether or to remove the outliers. Defaults to True.
        z_threshold (int, optional): z_score threshold to eliminate outliers. Defaults to 4.
        interpolate (bool, optional): choose whether or to interpolate the data. Defaults to False.
        frequency (str, optional): {‘D’, ‘W’, ‘M’, ‘Y’} frequency to interpolate. Note: See https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html for valid frequency inputs. (e.g. ‘W’ = every week, ‘D ’= every day, ‘2W’ = every 2 weeks). Defaults to '2W'.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_correlation'.
        log_transform (bool, optional): flag for log base 10 transformation. Defaults to False.
        fontsize (int, optional): font size. Defaults to 20.
    """
    
    # data = self.data
    data = data_pylenm_dm.data
    
    stations = np.array(data.STATION_ID.values)
    stations = np.unique(stations)
    
    for station in stations:
        plot_corr_by_station(
            data_pylenm_dm=data_pylenm_dm,
            station_name=station, 
            analytes=analytes,remove_outliers=remove_outliers, 
            z_threshold=z_threshold, 
            interpolate=interpolate, 
            frequency=frequency, 
            save_dir=save_dir, 
            log_transform=log_transform, 
            fontsize=fontsize,
        )

 
def plot_corr_by_date_range(
        # self, 
        data_pylenm_dm, 
        date, 
        analytes, 
        lag=0, 
        min_samples=10, 
        save_dir='plot_corr_by_date', 
        log_transform=False, 
        fontsize=20, 
        # returnData=False, 
        no_log=None,
    ):
    """Plots the correlations with the physical plots as station as the correlations of the important analytes for ALL the stations on a specified date or range of dates if a lag greater than 0 is specifed.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        date (str): date to be analyzed
        analytes (_type_): list of analyte names to use
        lag (int, optional): number of days to look ahead and behind the specified date (+/-). Defaults to 0.
        min_samples (int, optional): minimum number of samples the result should contain in order to execute.. Defaults to 10.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_corr_by_date'.
        log_transform (bool, optional): flag for log base 10 transformation. Defaults to False.
        fontsize (int, optional): font size. Defaults to 20.
        returnData (bool, optional): flag to return data used to perfrom correlation analysis. Defaults to False.
        no_log (list, optional): list of column names to not apply log transformation to. Defaults to None.
    """
    
    if (lag==0):
        
        # Prepare data
        data = data_pylenm_dm.data
        data = filters.simplify_data(data=data)
        query = data[data.COLLECTION_DATE == date]
        a = list(np.unique(query.ANALYTE_NAME.values))# get all analytes from dataset
        
        for value in analytes:
            # if((value in a)==False):
            if value not in a:
                plot_corr_logger.error(f"ERROR: No analyte named {value} in data.")
                # return 'ERROR: No analyte named "{}" in data.'.format(value)
                # return f"ERROR: No analyte named {value} in data."
                return None
        
        analytes = sorted(analytes)
        query = query.loc[query.ANALYTE_NAME.isin(analytes)]
        
        if(query.shape[0] == 0):
            plot_corr_logger.error(f"ERROR: {date} has no data for all of the analytes.")
            # return 'ERROR: {} has no data for all of the analytes.'.format(date)
            # return f"ERROR: {date} has no data for all of the analytes."
            return None
        
        samples = query[['COLLECTION_DATE', 'STATION_ID', 'ANALYTE_NAME']].duplicated().value_counts()[0]
        
        # Transform data if there are enough samples
        if(samples < min_samples):
            plot_corr_logger.error(f"ERROR: {date} does not have at least {min_samples} samples.")
            # return 'ERROR: {} does not have at least {} samples.'.format(date, min_samples)
            # return f"ERROR: {date} does not have at least {min_samples} samples."
            return None
        else:
            piv = query.reset_index().pivot_table(
                index='STATION_ID', 
                columns='ANALYTE_NAME', 
                values='RESULT', 
                aggfunc=np.mean,
            )
            # return piv
    
    else:
        # If the data has already been calculated with the lag specified, retrieve it
        # if(self.jointData_is_set(lag=lag)==True): 
        #     data = self.__jointData[0]
        if data_pylenm_dm.is_set_jointData(lag=lag): 
            data = data_pylenm_dm.jointData[0]
        
        # Otherwise, calculate it
        else:
            # data = self.getJointData(analytes, lag=lag)
            # self.__set_jointData(data=data, lag=lag)
            data = fetchers.getJointData(data_pylenm_dm=data_pylenm_dm, analytes=analytes, lag=lag)
            data_pylenm_dm.set_jointData(data=data, lag=lag)
        
        # Get new range based on the lag and create the pivot table to be able to do the correlation
        # dateStart, dateEnd = self.__getLagDate(date, lagDays=lag)
        dateStart, dateEnd = fetchers._getLagDate(date, lagDays=lag)
        dateRange_key = str(dateStart.date()) + " - " + str(dateEnd.date())
        piv = pd.DataFrame(data.loc[dateRange_key]).unstack().T
        piv.index = piv.index.droplevel()
        piv = pd.DataFrame(piv).dropna(axis=0, how='all')
        num_NaNs = int(piv.isnull().sum().sum())
        samples = (piv.shape[0] * piv.shape[1]) - num_NaNs
    
        # Convert columns to float and raise errors if the conversion fails
        for col in piv.columns:
            piv[col] = piv[col].astype('float64', errors='raise')
    
        if lag > 0:
            date = dateRange_key
        # return piv
    
    title = date + '_correlation'
    # scaler = StandardScaler()
    # pivScaled = scaler.fit_transform(piv)
    # pivScaled = pd.DataFrame(pivScaled, columns=piv.columns)
    # pivScaled.index = piv.index
    # piv = pivScaled

    # Perform log transform (if specified)
    if(log_transform):
        piv[piv <= 0] = 0.00000001
        temp = piv.copy()
        piv = np.log10(piv)
    
        if(no_log !=None):
            for col in no_log:
                piv[col] = temp[col]

    # Create seaborn PairGrid plot
    sns.set_style("white", {"axes.facecolor": "0.95"})
    g = sns.PairGrid(piv, aspect=1.2, diag_sharey=False, despine=False)
    g.figure.suptitle(title, fontweight='bold', y=1.08, fontsize=25)
    g.map_lower(
        sns.regplot, 
        lowess=True, 
        ci=False, 
        line_kws={'color': 'red', 'lw': 3}, 
        scatter_kws={'color': 'black', 's': 20},
    )
    # g.map_diag(
    #     sns.distplot, 
    #     kde_kws={'color': 'black', 'lw': 3}, 
    #     hist_kws={'histtype': 'bar', 'lw': 2, 'edgecolor': 'k', 'facecolor':'grey'},
    # )   # NOTE: Replacing with `sns.histplot` because `sns.distplot` is deprecated.
    g.map_diag(
        sns.histplot, 
        stat="density",
        color="black",
        element="bars",
        edgecolor="dimgray",
        facecolor="lightgray",
        lw=2,
        kde=True, 
        kde_kws=dict(cut=3),
        line_kws=dict(color="black", lw=3),
    )
    # g.map_upper(self.__plotUpperHalf)
    g.map_upper(plots._plotUpperHalf)
    
    # for ax in g.axes.flat:
    #     ax.tick_params("y", labelrotation=0, labelsize=fontsize)
    #     # ax.set_xticklabels(
    #     #     ax.get_xticklabels(), rotation=45, fontsize=fontsize,
    #     # )   # DISCOURAGED - https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xticklabels.html
    #     ax.set_xticks(
    #         ax.get_xticks(),
    #         [f"{num:.2f}" for num in ax.get_xticks()],    # TODO: Check with Zexuan if this is supposed to be normalized.
    #         rotation=45, 
    #         fontsize=fontsize,
    #     )
    #     ax.set_xlabel(
    #         ax.get_xlabel(), fontsize=fontsize, fontweight='bold',
    #     ) #HERE
    #     ax.set_ylabel(
    #         ax.get_ylabel(), fontsize=fontsize, fontweight='bold',
    #     )
    _format_ticks_and_labels(axes=g.axes, fontsize=fontsize)

    # g.figure.subplots_adjust(wspace=0.3, hspace=0.3)
    g.figure.subplots_adjust(
        wspace=0.3, hspace=0.3, 
        left=0.1, right=0.75, 
        bottom=0.1, top=0.9,
    )
    ax = plt.gca()

    # Add the meta-data information box
    props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
    ax.text(
        1.3, 6.2, 
        f"Date: {date}\n\nStations: {piv.shape[0]}\nSamples used: {samples}",
        # 'Date:  {}\n\nStations:     {}\nSamples used:     {}'.format(date, piv.shape[0] ,samples), 
        transform=ax.transAxes, 
        fontsize=fontsize, 
        fontweight='bold', 
        verticalalignment='bottom', 
        bbox=props,
    )
    
    # Add titles to the diagonal axes/subplots
    for ax, col in zip(np.diag(g.axes), piv.columns):
        ax.set_title(col, y=0.82, fontsize=15)
    
    # Save figure (if specified)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    g.figure.savefig(
        # save_dir + '/' + date + '.png',
        f"{save_dir}/{date}.png", 
        bbox_inches="tight",
    )
    
    # if(returnData):
    #     return piv

    return piv


def plot_corr_by_year(
        # self, 
        data_pylenm_dm,
        year, 
        analytes, 
        remove_outliers=True, 
        z_threshold=4, 
        min_samples=10, 
        save_dir='plot_corr_by_year', 
        log_transform=False, 
        fontsize=20, 
        returnData=False, 
        no_log=None,
    ):
    """Plots the correlations with the physical plots as station as the correlations of the important analytes for ALL the stations in specified year.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        year (int): year to be analyzed
        analytes (list): list of analyte names to use
        remove_outliers (bool, optional): choose whether or to remove the outliers.. Defaults to True.
        z_threshold (int, optional): z_score threshold to eliminate outliers. Defaults to 4.
        min_samples (int, optional): minimum number of samples the result should contain in order to execute.. Defaults to 10.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_correlation'.
        log_transform (bool, optional): flag for log base 10 transformation. Defaults to False.
        fontsize (int, optional): font size. Defaults to 20.
        returnData (bool, optional): flag to return data used to perfrom correlation analysis. Defaults to False.
        no_log (list, optional): list of column names to not apply log transformation to. Defaults to None.
    """
    
    # Prepare data
    # data = self.data
    data = data_pylenm_dm.data
    query = data
    # query = self.simplify_data(data=query)
    query = filters.simplify_data(data=query)
    query.COLLECTION_DATE = pd.to_datetime(query.COLLECTION_DATE)
    query = query[query.COLLECTION_DATE.dt.year == year]
    a = list(np.unique(query.ANALYTE_NAME.values))# get all analytes from dataset
    
    for value in analytes:
        if((value in a)==False):
            plot_corr_logger.error(f"ERROR: No analyte named {value} in data.")
            # return 'ERROR: No analyte named "{}" in data.'.format(value)
            return None
    
    analytes = sorted(analytes)
    
    # Performing checks for enough data for the plots
    query = query.loc[query.ANALYTE_NAME.isin(analytes)]
    if(query.shape[0] == 0):
        plot_corr_logger.error(f"ERROR: {year} has no data for the 6 analytes.")
        # return 'ERROR: {} has no data for the 6 analytes.'.format(year)
        return None

    samples = query[['COLLECTION_DATE', 'STATION_ID', 'ANALYTE_NAME']].duplicated().value_counts()[0]
    if(samples < min_samples):
        plot_corr_logger.error(f"ERROR: {year} does not have at least {min_samples} samples.")
        # return 'ERROR: {} does not have at least {} samples.'.format(date, min_samples)
        return None

    else:
        # Create pivot table if there is enough data
        piv = query.reset_index().pivot_table(
            index='STATION_ID', 
            columns='ANALYTE_NAME', 
            values='RESULT', 
            aggfunc=np.mean,
        )
        # return piv
        
        # Remove outliers
        if(remove_outliers):
            # piv = self.remove_outliers(piv, z_threshold=z_threshold)
            piv = preprocess.remove_outliers(piv, z_threshold=z_threshold)
        
        samples = piv.shape[0] * piv.shape[1]

        title = str(year) + '_correlation'
        # scaler = StandardScaler()
        # pivScaled = scaler.fit_transform(piv)
        # pivScaled = pd.DataFrame(pivScaled, columns=piv.columns)
        # pivScaled.index = piv.index
        # piv = pivScaled

        # Perform log transformation (if specified)
        if(log_transform):
            piv[piv <= 0] = 0.00000001
            temp = piv.copy()
            piv = np.log10(piv)
            
            if(no_log !=None):
                for col in no_log:
                    piv[col] = temp[col]

        # Create the seaborn PairGrid plot
        sns.set_style("white", {"axes.facecolor": "0.95"})
        g = sns.PairGrid(piv, aspect=1.2, diag_sharey=False, despine=False)
        g.figure.suptitle(title, fontweight='bold', y=1.08, fontsize=25)
        g.map_lower(
            sns.regplot, 
            lowess=True, 
            ci=False, 
            line_kws={'color': 'red', 'lw': 3}, 
            scatter_kws={'color': 'black', 's': 20},
        )
        # g.map_diag(
        #     sns.distplot, 
        #     kde_kws={'color': 'black', 'lw': 3}, 
        #     hist_kws={'histtype': 'bar', 'lw': 2, 'edgecolor': 'k', 'facecolor':'grey'},
        # )   # NOTE: Replacing with `sns.histplot` because `sns.distplot` is deprecated.
        g.map_diag(
            sns.histplot, 
            stat="density",
            color="black",
            element="bars",
            edgecolor="dimgray",
            facecolor="lightgray",
            lw=2,
            kde=True, 
            kde_kws=dict(cut=3),
            line_kws=dict(color="black", lw=3),
        )
        # g.map_upper(self.__plotUpperHalf)
        g.map_upper(plots._plotUpperHalf)
        
        # for ax in g.axes.flat:
        #     ax.tick_params("y", labelrotation=0, labelsize=fontsize)
        #     # ax.set_xticklabels(
        #     #     ax.get_xticklabels(), rotation=45, fontsize=fontsize,
        #     # )   # DISCOURAGED - https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xticklabels.html
        #     ax.set_xticks(
        #         ax.get_xticks(),
        #         [f"{num:.2f}" for num in ax.get_xticks()],    # TODO: Check with Zexuan if this is supposed to be normalized.
        #         rotation=45, 
        #         fontsize=fontsize,
        #     )
        #     ax.set_xlabel(
        #         ax.get_xlabel(), fontsize=fontsize, fontweight='bold',
        #     )
        #     ax.set_ylabel(
        #         ax.get_ylabel(), fontsize=fontsize, fontweight='bold',
        #     )
        _format_ticks_and_labels(axes=g.axes, fontsize=fontsize)

        # g.figure.subplots_adjust(wspace=0.3, hspace=0.3)
        g.figure.subplots_adjust(
            wspace=0.3, hspace=0.3, 
            left=0.1, right=0.75, 
            bottom=0.1, top=0.9,
        )

        # Add the meta-data information box
        ax = plt.gca()
        props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
        ax.text(
            1.3, 
            # 3, 
            6.2, 
            # 'Date:  {}\n\nSamples used:     {}'.format(year, samples), 
            f"Year: {year}\n\nSamples used: {samples}",
            transform=ax.transAxes, 
            fontsize=fontsize, 
            fontweight='bold', 
            verticalalignment='bottom', 
            bbox=props,
        )
        
        # Add titles to the diagonal axes/subplots
        for ax, col in zip(np.diag(g.axes), piv.columns):
            ax.set_title(col, y=0.82, fontsize=15)
        
        # Save figure (if specified)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        g.figure.savefig(save_dir + '/' + str(year) + '.png', bbox_inches="tight")
        
        # if(returnData):
        #     return piv

        return piv