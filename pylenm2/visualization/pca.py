import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from pylenm2.data import filters, fetchers, transformation
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

plot_pca_logger = logger_config.setup_logging(
    module_name=__name__,
    level=logging.INFO,
    # level=logging.DEBUG,
    logfile_dir=c.LOGFILE_DIR,
)


def plot_PCA_by_date(
        # self, 
        data_pylenm_dm,
        date, 
        analytes, 
        lag=0, 
        n_clusters=4, 
        return_clusters=False, 
        min_samples=3, 
        show_labels=True, 
        save_dir='plot_PCA_by_date', 
        filter=False, 
        col=None, 
        equals=[],
    ):
    """Generates a PCA biplot (PCA score plot + loading plot) of the data given a date in the dataset. The data is also clustered into n_clusters.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        date (str): date to be analyzed
        analytes (str): list of analyte names to use
        lag (int, optional): number of days to look ahead and behind the specified date (+/-). Defaults to 0.
        n_clusters (int, optional): number of clusters to split the data into.. Defaults to 4.
        return_clusters (bool, optional): Flag to return the cluster data to be used for spatial plotting.. Defaults to False.
        min_samples (int, optional): minimum number of samples the result should contain in order to execute.. Defaults to 3.
        show_labels (bool, optional): choose whether or not to show the name of the stations.. Defaults to True.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_PCA_by_date'.
        filter (bool, optional): flag to indicate filtering. Defaults to False.
        col (str, optional): column to filter. Example: col='STATION_ID'. Defaults to None.
        equals (list, optional): values to filter col by. Examples: equals=['FAI001A', 'FAI001B']. Defaults to [].
    """
    
    if(lag==0):
        
        # Prepare data
        # data = self.data
        data = data_pylenm_dm.data
        # data = self.simplify_data(data=data)
        data = filters.simplify_data(data=data)
        query = data[data.COLLECTION_DATE == date]
        
        # Filter data (if specified)
        if(filter):
            # filter_res = self.filter_by_column(data=self.construction_data, col=col, equals=equals)
            filter_res = filters.filter_by_column(
                data=data_pylenm_dm.construction_data, 
                col=col, 
                equals=equals,
            )
            
            # if('ERROR:' in str(filter_res)):
            #     return filter_res
            if filter_res is None:
                return None

            query_stations = list(query.STATION_ID.unique())
            filter_stations = list(filter_res.index.unique())
            intersect_stations = list(set(query_stations) & set(filter_stations))
            if len(intersect_stations) <= 0:
                plot_pca_logger.error("ERROR: No results for this query with the specifed filter parameters.")
                # return 'ERROR: No results for this query with the specifed filter parameters.'
                return None
            
            query = query[query['STATION_ID'].isin(intersect_stations)]
        
        a = list(np.unique(query.ANALYTE_NAME.values))  # get all analytes from dataset
        
        for value in analytes:
            if((value in a)==False):
                plot_pca_logger.error(f"ERROR: No analyte named {value} in data.")
                # return 'ERROR: No analyte named "{}" in data.'.format(value)
                return None
        
        analytes = sorted(analytes)
        query = query.loc[query.ANALYTE_NAME.isin(analytes)]

        if query.shape[0] == 0:
            plot_pca_logger.error(f"ERROR: {date} has no data for the {len(analytes)} analytes.")
            # return 'ERROR: {} has no data for the {len(analytes)} analytes.'.format(date)
            return None
        
        samples = query[['COLLECTION_DATE', 'STATION_ID', 'ANALYTE_NAME']].duplicated().value_counts()[0]
        
        # Check if there are enough samples
        if samples < min_samples:
            plot_pca_logger.error(f"ERROR: {date} does not have at least {min_samples} samples.")
            # return 'ERROR: {} does not have at least {} samples.'.format(date, min_samples)
            return None

        # if(len(np.unique(query.ANALYTE_NAME.values)) < 6):
        #     return 'ERROR: {} has less than the {len(analytes)} analytes we want to analyze.'.format(date)
        else:
            # analytes = preprocess._custom_analyte_sort(np.unique(query.ANALYTE_NAME.values))
            analytes = sorted(analytes)
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
            data = fetchers.getJointData(
                data_pylenm_dm=data_pylenm_dm, 
                analytes=analytes, 
                lag=lag,
            )
            # self.__set_jointData(data=data, lag=lag)
            data_pylenm_dm.set_jointData(data=data, lag=lag)
        
        # get new range based on the lag and create the pivor table to be able to do the correlation
        # dateStart, dateEnd = self.__getLagDate(date, lagDays=lag)
        dateStart, dateEnd = fetchers._getLagDate(date, lagDays=lag)
        # dateRange_key = str(dateStart.date()) + " - " + str(dateEnd.date())
        dateRange_key = f"{str(dateStart.date())} - {str(dateEnd.date())}"
        piv = pd.DataFrame(data.loc[dateRange_key]).unstack().T
        piv.index = piv.index.droplevel()
        piv = pd.DataFrame(piv).dropna(axis=0, how='all')
        num_NaNs = int(piv.isnull().sum().sum())
        samples = (piv.shape[0]*piv.shape[1])-num_NaNs
        for col in piv.columns:
            piv[col] = piv[col].astype('float64', errors='raise')
        if lag > 0:
            date = dateRange_key
        # return piv

    main_data = piv.dropna()
    
    # Fit PCA and KMeans models
    scaler = StandardScaler()
    X = scaler.fit_transform(main_data)
    
    pca = PCA(n_components=2)
    x_new = pca.fit_transform(X)
    
    pca_points = pd.DataFrame(x_new, columns=["x1", "x2"])
    k_Means = KMeans(n_init='auto', n_clusters=n_clusters, random_state=42)
    model = k_Means.fit(pca_points[['x1', 'x2']])
    predict = model.predict(pca_points[['x1', 'x2']])
    
    # attach predicted cluster to original points
    pca_points['predicted'] = model.labels_
    
    # Create a dataframe for cluster_centers (centroids)
    centroids = pd.DataFrame(model.cluster_centers_, columns=["x1", "x2"])
    colors = [
        'red', 'blue', 'orange', 'purple', 'green', 'beige', 
        'pink', 'black', 'cadetblue', 'lightgreen',
    ]
    pca_points['color'] = pca_points['predicted'].map(lambda p: colors[p])

    fig, ax = plt.subplots(figsize=(10,10))
    # ax = plt.axes()   # NOTE: Causes redundant overlapping x and y ticks.

    small_fontSize = 15
    large_fontSize = 20
    plt.rc('axes', titlesize=large_fontSize)
    plt.rc('axes', labelsize=large_fontSize)
    plt.rc('legend', fontsize=small_fontSize)
    plt.rc('xtick', labelsize=small_fontSize)
    plt.rc('ytick', labelsize=small_fontSize)

    def myplot(score, coeff, labels=None, c='r', centroids=None):
        xs = score.iloc[:,0]
        ys = score.iloc[:,1]
        n = coeff.shape[0]
        scalex = 1.0/(xs.max() - xs.min())
        scaley = 1.0/(ys.max() - ys.min())
        scatt_X = xs * scalex
        scatt_Y = ys * scaley
        scatter = plt.scatter(
            scatt_X, 
            scatt_Y, 
            alpha=0.8, 
            label='Stations', 
            c=c,
        )
        centers = plt.scatter(
            centroids.iloc[:,0] * scalex, 
            centroids.iloc[:,1] * scaley, 
            c=colors[0:n_clusters], 
            marker='X', 
            s=550,
        )

        for i in range(n):
            arrow = plt.arrow(
                0, 0, coeff[i,0], coeff[i,1], color='r', alpha=0.9, 
                head_width=0.05, head_length=0.05, label='Loadings',
            )
            if labels is None:
                plt.text(
                    coeff[i,0]* 1.15, coeff[i,1] * 1.15, f"Var{i+1}", 
                    color='g', ha='center', va='center',
                )
            else:
                plt.text(
                    coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], 
                    color='g', ha='center', va='bottom',
                )

        if(show_labels):
            for x_pos, y_pos, label in zip(scatt_X, scatt_Y, main_data.index):
                ax.annotate(
                    label, # The label for this point
                    xy=(x_pos, y_pos), # Position of the corresponding point
                    xytext=(7, 0),     # Offset text by 7 points to the right
                    textcoords='offset points', # tell it to use offset points
                    ha='left',         # Horizontally aligned to the left
                    va='center',       # Vertical alignment is centered
                    color='black', 
                    alpha=0.8,
                )
        plt.legend(
            [scatter, centers, arrow], 
            ['Stations', 'Station centroids','Loadings'],
        )

    samples = x_new.shape[0] * piv.shape[1]

    # Add the meta-data information box
    props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
    ax.text(
        1.1, 
        0.5, 
        # 'Date:  {}\n\nSamples:          {}\nStations:               {}'.format(date, samples, x_new.shape[0]), 
        f"Date: {date}\n\nSamples: {samples}\nStations: {x_new.shape[0]}",
        transform=ax.transAxes, 
        fontsize=20, 
        fontweight='bold', 
        verticalalignment='bottom', 
        bbox=props,
    )

    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.xlabel("PC{}".format(1))
    plt.ylabel("PC{}".format(2))
    # ax.set_title('PCA Biplot - ' + date, fontweight='bold')
    ax.set_title(f"PCA Biplot - {date}", fontweight='bold')
    plt.grid(alpha=0.5)

    #Call the function. Use only the 2 Principal Components.
    myplot(
        pca_points, 
        np.transpose(pca.components_[0:2, :]), 
        labels=piv.columns, 
        c=pca_points['color'], 
        centroids=centroids,
    )
    plt.show()
    
    # Save figure (if specified)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    fig.savefig(
        f"{save_dir}/PCA Biplot - {date}.png", bbox_inches="tight",
        # save_dir + '/' + 'PCA Biplot - '+ date +'.png', bbox_inches="tight"
    )
    
    # Return clusters (if specified)
    if(return_clusters):
        stations = list(main_data.index)
        color_stations = list(pca_points.color)
        
        def merge(list1, list2): 
            merged_list = [(list1[i], list2[i]) for i in range(0, len(list1))] 
            return merged_list
        
        color_df = pd.DataFrame(
            merge(stations, color_stations), 
            columns=['STATION_ID', 'color'],
        )
        
        # if(self.get_Construction_Data==None):
        if data_pylenm_dm.get_construction_data is None:
            plot_pca_logger.error("You need to set the GPS data first for using the get_construction_data() function.")
            # print('You need to set the GPS data first using the getConstructionData function.')
            return None
        
        else:
            gps_color = pd.merge(
                data_pylenm_dm.get_construction_data(), 
                color_df, 
                on=['STATION_ID'],
            )
            return gps_color


def plot_PCA_by_year(
        # self, 
        data_pylenm_dm,
        year, 
        analytes, 
        n_clusters=4, 
        return_clusters=False, 
        min_samples=10, 
        show_labels=True, 
        save_dir='plot_PCA_by_year', 
        filter=False, 
        col=None, 
        equals=[],
    ):
    """Gernates a PCA biplot (PCA score plot + loading plot) of the data given a year in the dataset. The data is also clustered into n_clusters.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        year (int): year to be analyzed
        analytes (str): list of analyte names to use
        n_clusters (int, optional): number of clusters to split the data into.. Defaults to 4.
        return_clusters (bool, optional): Flag to return the cluster data to be used for spatial plotting.. Defaults to False.
        min_samples (int, optional): minimum number of samples the result should contain in order to execute.. Defaults to 3.
        show_labels (bool, optional): choose whether or not to show the name of the stations.. Defaults to True.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_PCA_by_date'.
        filter (bool, optional): flag to indicate filtering. Defaults to False.
        col (str, optional): column to filter. Example: col='STATION_ID'. Defaults to None.
        equals (list, optional): values to filter col by. Examples: equals=['FAI001A', 'FAI001B']. Defaults to [].
    """

    # Prepare data
    # data = self.data
    data = data_pylenm_dm.data
    # query = self.simplify_data(data=data)
    query = filters.simplify_data(data=data)
    query.COLLECTION_DATE = pd.to_datetime(query.COLLECTION_DATE)
    query = query[query.COLLECTION_DATE.dt.year == year]
    
    # Filter data (if specified)
    if(filter):
        # filter_res = self.filter_by_column(data=self.construction_data, col=col, equals=equals)
        filter_res = filters.filter_by_column(
            data=data_pylenm_dm.construction_data, 
            col=col, 
            equals=equals,
        )
        
        # if('ERROR:' in str(filter_res)):
        #     return filter_res
        if filter_res is None:
            return None

        query_stations = list(query.STATION_ID.unique())
        filter_stations = list(filter_res.index.unique())
        intersect_stations = list(set(query_stations) & set(filter_stations))

        if len(intersect_stations) <= 0:
            plot_pca_logger.error("ERROR: No results for this query with the specifed filter parameters.")
            # return 'ERROR: No results for this query with the specifed filter parameters.'
            return None
        
        query = query[query['STATION_ID'].isin(intersect_stations)]
    
    a = list(np.unique(query.ANALYTE_NAME.values))  # get all analytes from dataset
    
    for value in analytes:
        if((value in a)==False):
            plot_pca_logger.error(f"ERROR: No analyte named {value} in data.")
            # return 'ERROR: No analyte named "{}" in data.'.format(value)
            return None
    
    analytes = sorted(analytes)
    query = query.loc[query.ANALYTE_NAME.isin(analytes)]

    if(query.shape[0] == 0):
        plot_pca_logger.error(f"ERROR: {year} has no data for the {len(analytes)} analytes.")
        # return 'ERROR: {} has no data for the {len(analytes)} analytes.'.format(year)
        return None
    
    samples = query[['COLLECTION_DATE', 'STATION_ID', 'ANALYTE_NAME']].duplicated().value_counts()[0]

    if samples < min_samples:
        plot_pca_logger.error(f"ERROR: {year} does not have at least {min_samples} samples.")
        # return 'ERROR: {} does not have at least {} samples.'.format(year, min_samples)
        return None

    # if(len(np.unique(query.ANALYTE_NAME.values)) < 6):
    #     return 'ERROR: {} has less than the {len(analytes)} analytes we want to analyze.'.format(year)
    else:
        # analytes = preprocess._custom_analyte_sort(np.unique(query.ANALYTE_NAME.values))
        analytes = sorted(analytes)
        piv = query.reset_index().pivot_table(
            index='STATION_ID', 
            columns='ANALYTE_NAME', 
            values='RESULT', 
            aggfunc=np.mean,
        )
        
        main_data = piv.dropna()
        # # FILTERING CODE
        # if(filter):
        #     res_stations = self.filter_stations(filter_station_by)
        #     main_data = main_data.loc[main_data.index.isin(res_stations)]
        
        # Fit PCA and KMeans models
        scaler = StandardScaler()
        X = scaler.fit_transform(main_data)
        
        pca = PCA(n_components=2)
        x_new = pca.fit_transform(X)

        pca_points = pd.DataFrame(x_new, columns=["x1", "x2"])
        k_Means = KMeans(n_init='auto', n_clusters=n_clusters, random_state=42)
        model = k_Means.fit(pca_points[['x1', 'x2']])
        predict = model.predict(pca_points[['x1', 'x2']])
        
        # attach predicted cluster to original points
        pca_points['predicted'] = model.labels_
        
        # Create a dataframe for cluster_centers (centroids)
        centroids = pd.DataFrame(model.cluster_centers_, columns=["x1", "x2"])
        colors = [
            'red', 'blue', 'orange', 'purple', 'green', 'beige', 
            'pink', 'black', 'cadetblue', 'lightgreen',
        ]
        pca_points['color'] = pca_points['predicted'].map(lambda p: colors[p])
        
        fig, ax = plt.subplots(figsize=(15,15))
        # ax = plt.axes()   # NOTE: Causes redundant overlapping x and y ticks.

        small_fontSize = 15
        large_fontSize = 20
        plt.rc('axes', titlesize=large_fontSize)
        plt.rc('axes', labelsize=large_fontSize)
        plt.rc('legend', fontsize=small_fontSize)
        plt.rc('xtick', labelsize=small_fontSize)
        plt.rc('ytick', labelsize=small_fontSize) 

        def myplot(score,coeff,labels=None,c='r', centroids=None):
            xs = score[:,0]
            ys = score[:,1]
            n = coeff.shape[0]
            scalex = 1.0/(xs.max() - xs.min())
            scaley = 1.0/(ys.max() - ys.min())
            scatt_X = xs * scalex
            scatt_Y = ys * scaley
            scatter = plt.scatter(
                scatt_X, 
                scatt_Y, 
                alpha=0.8, 
                label='Stations', 
                c=c,
            )
            centers = plt.scatter(
                centroids.iloc[:,0] * scalex, 
                centroids.iloc[:,1] * scaley, 
                c = colors[0:n_clusters], 
                marker='X', 
                s=550,
            )

            for i in range(n):
                arrow = plt.arrow(
                    0, 0, coeff[i,0], coeff[i,1], color='r', alpha=0.9, 
                    head_width=0.05, head_length=0.05,
                )
                if labels is None:
                    plt.text(
                        coeff[i,0] * 1.15, coeff[i,1] * 1.15, f"Var{i+1}", 
                        color='g', ha='center', va='center',
                    )
                else:
                    plt.text(
                        coeff[i,0] * 1.15, coeff[i,1] * 1.15, labels[i], 
                        color='g', ha='center', va='bottom',
                    )

            if(show_labels):
                for x_pos, y_pos, label in zip(scatt_X, scatt_Y, main_data.index):
                    ax.annotate(
                        label, # The label for this point
                        xy=(x_pos, y_pos), # Position of the corresponding point
                        xytext=(7, 0),     # Offset text by 7 points to the right
                        textcoords='offset points', # tell it to use offset points
                        ha='left',         # Horizontally aligned to the left
                        va='center', 
                        color='black', 
                        alpha=0.8,
                    )       # Vertical alignment is centered
            plt.legend(
                [scatter, centers, arrow], 
                ['Stations', 'Station centroids','Loadings'],
            )

        samples = x_new.shape[0]*piv.shape[1]    
        
        # Add meta-data information box
        props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
        ax.text(
            1.1, 
            0.5, 
            # 'Date:  {}\n\nSamples:          {}\nStations:               {}'.format(year,samples, x_new.shape[0]), 
            f"Year: {year}\n\nSamples: {samples}\nStations: {x_new.shape[0]}",
            transform=ax.transAxes, 
            fontsize=20, 
            fontweight='bold', 
            verticalalignment='bottom', 
            bbox=props,
        )

        plt.xlim(-1,1)
        plt.ylim(-1,1)
        plt.xlabel("PC{}".format(1))
        plt.ylabel("PC{}".format(2))
        # ax.set_title('PCA Biplot - ' + str(year), fontweight='bold')
        ax.set_title(f"PCA Biplot - {str(year)}", fontweight='bold')
        plt.grid(alpha=0.5)

        #Call the function. Use only the 2 PCs.
        myplot(
            x_new[:,0:2], 
            np.transpose(pca.components_[0:2, :]), 
            labels=piv.columns, 
            c=pca_points['color'], 
            centroids=centroids,
        )
        plt.show()

        # Save figure (if specified)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        fig.savefig(
            f"{save_dir}/PCA Biplot - {year}.png", bbox_inches="tight",
            # save_dir + '/' + 'PCA Biplot - '+ str(year) +'.png', bbox_inches="tight",
        )
        
        # Return clusters (if specified)
        if(return_clusters):
            stations = list(main_data.index)
            color_stations = list(pca_points.color)
            
            def merge(list1, list2): 
                merged_list = [
                    (list1[i], list2[i]) 
                    for i in range(0, len(list1))
                ] 
                return merged_list
            
            color_df = pd.DataFrame(
                merge(stations, color_stations), 
                columns=['STATION_ID', 'color'],
            )
            
            # if(self.get_Construction_Data==None):
            if data_pylenm_dm.get_construction_data is None:
                plot_pca_logger.error('You need to set the GPS data first for using the setConstructionData function.')
                # print('You need to set the GPS data first for using the setConstructionData function.')
                return None
            
            else:
                gps_color = pd.merge(
                    data_pylenm_dm.get_construction_data(), 
                    color_df, 
                    on=['STATION_ID'],
                )
                return gps_color


def plot_PCA_by_station(
        # self, 
        data_pylenm_dm, 
        station_name, 
        analytes, 
        interpolate=False, 
        frequency='2W', 
        min_samples=10, 
        show_labels=True, 
        save_dir='plot_PCA_by_station',
    ):
    """Gernates a PCA biplot (PCA score plot + loading plot) of the data given a station_name in the dataset. Only uses the 6 important analytes.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        station_name (str): name of the station to be processed
        analytes (str): list of analyte names to use
        interpolate (bool, optional): choose to interpolate the data. Defaults to False.
        frequency (str, optional): {‘D’, ‘W’, ‘M’, ‘Y’} frequency to interpolate. See https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html for valid frequency inputs. (e.g. ‘W’ = every week, ‘D ’= every day, ‘2W’ = every 2 weeks). Defaults to '2W'.
        min_samples (int, optional): minimum number of samples the result should contain in order to execute.. Defaults to 3.
        show_labels (bool, optional): choose whether or not to show the name of the stations.. Defaults to True.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_PCA_by_date'.
    """

    # data = self.data
    data = data_pylenm_dm.data
    query = data[data.STATION_ID == station_name]
    a = list(np.unique(query.ANALYTE_NAME.values))# get all analytes from dataset
    
    for value in analytes:
        if((value in a)==False):
            plot_pca_logger.error(f"ERROR: No analyte named {value} in data.")
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
    totalSamples = piv.stack().shape[0]
    piv = piv.dropna()
    
    # Interpolate (if specified)
    if(interpolate):
        # piv = self.interpolate_station_data(station_name, analytes, frequency=frequency)
        piv = transformation.interpolate_station_data(
            station_name, 
            analytes, 
            frequency=frequency,
        )
        # title = 'PCA Biplot - ' + station_name + ' - interpolated every ' + frequency
        title = f"PCA Biplot - {station_name} - interpolated every {frequency}."
    else:
        # title = 'PCA Biplot - ' + station_name
        title = f"PCA Biplot - {station_name}"

    if query.shape[0] == 0:
        plot_pca_logger.error(f"ERROR: Station {station_name} has no data for the {len(analytes)} analytes.")
        # return 'ERROR: {} has no data for the {len(analytes)} analytes.'.format(date)
        return None
    
    samples = query[['COLLECTION_DATE', 'STATION_ID', 'ANALYTE_NAME']].duplicated().value_counts()[0]
    if samples < min_samples:
        plot_pca_logger.error(f"ERROR: Station {station_name} does not have at least {min_samples} samples.")
        # return 'ERROR: {} does not have at least {} samples.'.format(date, min_samples)
        return None
    
    # if(len(np.unique(query.ANALYTE_NAME.values)) < 6):
    #     return 'ERROR: {} has less than the {len(analytes)} analytes we want to analyze.'.format(station_name)
    else:
        
        # Scale the values
        scaler = StandardScaler()
        X = scaler.fit_transform(piv.dropna())
        
        # Fit PCA model
        pca = PCA(n_components=2)
        x_new = pca.fit_transform(X)

        # Generate Plots
        fig, ax = plt.subplots(figsize=(15,15))
        # ax = plt.axes()   # NOTE: Causes redundance overlapping x and y ticks.

        small_fontSize = 15
        large_fontSize = 20
        plt.rc('axes', titlesize=large_fontSize)
        plt.rc('axes', labelsize=large_fontSize)
        plt.rc('legend', fontsize=small_fontSize)
        plt.rc('xtick', labelsize=small_fontSize)
        plt.rc('ytick', labelsize=small_fontSize) 

        def myplot(score, coeff, labels=None):
            xs = score[:,0]
            ys = score[:,1]
            n = coeff.shape[0]
            scalex = 1.0/(xs.max() - xs.min())
            scaley = 1.0/(ys.max() - ys.min())
            scatt_X = xs * scalex
            scatt_Y = ys * scaley
            scatter = plt.scatter(
                scatt_X, 
                scatt_Y, 
                alpha=0.8, 
                label='Date samples',
            )

            for i in range(n):
                arrow = plt.arrow(
                    0, 0, coeff[i,0], coeff[i,1], color = 'r', alpha = 0.9, 
                    head_width=0.05, head_length=0.05, label='Loadings',
                )
                if labels is None:
                    plt.text(
                        coeff[i,0] * 1.15, coeff[i,1] * 1.15, "Var"+str(i+1), 
                        color='g', ha='center', va='center',
                    )
                else:
                    plt.text(
                        coeff[i,0] * 1.15, coeff[i,1] * 1.15, labels[i], 
                        color='g', ha='center', va='bottom',
                    )

            if(show_labels):
                for x_pos, y_pos, label in zip(scatt_X, scatt_Y, piv.dropna().index.date):
                    ax.annotate(
                        label, # The label for this point
                        xy=(x_pos, y_pos), # Position of the corresponding point
                        xytext=(7, 0),     # Offset text by 7 points to the right
                        textcoords='offset points', # tell it to use offset points
                        ha='left',         # Horizontally aligned to the left
                        va='center',       # Vertical alignment is centered
                        color='black', 
                        alpha=0.8,
                    )
            plt.legend(
                [scatter, arrow], 
                ['Date samples', 'Loadings'],
            )

        samples = x_new.shape[0] * piv.shape[1]
        
        # Add the meta-data information box
        props = dict(boxstyle='round', facecolor='grey', alpha=0.15)      
        ax.text(
            1.1, 
            0.5, 
            # 'Start date:  {}\nEnd date:    {}\n\nOriginal samples:     {}\nSamples used:     {}\nDate samples:               {}'.format(piv.index[0].date(), piv.index[-1].date(), totalSamples, samples, x_new.shape[0]), 
            f"Start date: {piv.index[0].date()}\nEnd date: {piv.index[-1].date()}\n\nOriginal samples: {totalSamples}\nSamples used: {samples}\nDate samples: {x_new.shape[0]}",
            transform=ax.transAxes, 
            fontsize=20, 
            fontweight='bold', 
            verticalalignment='bottom', 
            bbox=props,
        )

        plt.xlim(-1,1)
        plt.ylim(-1,1)
        plt.xlabel("PC{}".format(1))
        plt.ylabel("PC{}".format(2))
        ax.set_title(title, fontweight='bold')
        plt.grid(alpha=0.5)

        #Call the function. Use only the 2 Principal Components.
        myplot(
            x_new[:,0:2], 
            np.transpose(pca.components_[0:2, :]), 
            labels=piv.columns,
        )
        plt.show()

        # Save figure (if specified)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        fig.savefig(
            f"{save_dir}/{title}.png", bbox_inches="tight",
            # save_dir + '/' + title +'.png', bbox_inches="tight"
        )