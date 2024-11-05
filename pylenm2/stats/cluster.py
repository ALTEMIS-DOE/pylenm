import pandas as pd
import matplotlib.dates as mdates
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt
from matplotlib.dates import date2num

from pylenm2.data import filters
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

cluster_logger = logger_config.setup_logging(
    module_name=__name__,
    level=logging.INFO,
    # level=logging.DEBUG,
    logfile_dir=c.LOGFILE_DIR,
)


# NOT in use!
def __cluster_data_OLD(self, data, n_clusters=4, log_transform=False, filter=False, filter_station_by=['D'], return_clusters=False):
    if(filter):
        res_stations = self.filter_stations(filter_station_by)
        data = data.T
        data = data.loc[data.index.isin(res_stations)]
        data = data.T
    if(log_transform):
        data = np.log10(data)
        data = data.dropna(axis=1)
    temp = data.T
    k_Means = KMeans(n_clusters=n_clusters, random_state=42)
    km = k_Means.fit(temp)
    predict = km.predict(temp)
    temp['predicted'] = km.labels_
    colors = ['red', 'blue', 'orange', 'purple', 'green', 'beige', 'pink', 'black', 'cadetblue', 'lightgreen']
    temp['color'] = temp['predicted'].map(lambda p: colors[p])

    fig, ax = plt.subplots(figsize=(20,10))
    ax = plt.axes()

    color = temp['color']
    for x in range(temp.shape[0]):
        curr = data.iloc[:,x]
        ax.plot(curr, label=curr.name, color=color[x])
        ax.legend()
    
    if(return_clusters):
        color_df = pd.DataFrame(temp['color'])
        color_df['STATION_ID'] = color_df.index
        if(self.get_Construction_Data==None):
            print('You need to set the GPS data first using the setConstructionData function.')
            return None
        else:
            gps_color = pd.merge(self.get_Construction_Data(), color_df, on=['STATION_ID'])
            return gps_color


def cluster_data(
        # self, 
        data_pylenm_dm,
        data, 
        analyte_name=["ANALYTE_NAME"], 
        n_clusters=4, 
        filter=False, 
        col=None, 
        equals=[], 
        year_interval=5, 
        y_label = 'Concentration', 
        return_clusters=False,
    ):
    """Clusters time series concentration data using kmeans algorithm and plots it.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        data (pd.DataFrame): data to be used in clustering.
        analyte_name (list, optional): analytes to use to cluster. Defaults to ["ANALYTE_NAME"].
        n_clusters (int, optional): number of clusters for kmeans. Defaults to 4.
        filter (bool, optional): flag to indicate filtering. Defaults to False.
        col (str, optional): column to filter. Example: col='STATION_ID'. Defaults to None.
        equals (list, optional): values to filter col by. Examples: equals=['FAI001A', 'FAI001B']. Defaults to [].
        year_interval (int, optional): plot x_label interval in years. Defaults to 5.
        y_label (str, optional): y axis label. Defaults to 'Concentration'.
        return_clusters (bool, optional): flag to return cluster assignemnt. Defaults to False.
    """
    
    data = data.copy()

    # Filter data (if specified)
    if(filter):
        # filter_res = self.filter_by_column(data=self.get_Construction_Data(), col=col, equals=equals)
        filter_res = filters.filter_by_column(
            data=data_pylenm_dm.get_construction_data(), 
            col=col, 
            equals=equals,
        )

        # if('ERROR:' in str(filter_res)):
        #     return filter_res
        if filter_res is None:
            return None
        
        query_stations = list(data.columns)
        filter_stations = list(filter_res.index.unique())
        intersect_stations = list(set(query_stations) & set(filter_stations))
        
        cluster_logger.info(f"{intersect_stations = }")
        cluster_logger.info(f"{len(intersect_stations) = }")
        
        if len(intersect_stations) <= 0:
            # return 'ERROR: No results for this query with the specifed filter parameters.'
            cluster_logger.error("ERROR: No results for this query with the specifed filter parameters.")
            return None

        data = data[intersect_stations]
    
    data.index = date2num(data.index)
    temp = data.T

    # KMeans clustering
    k_Means = KMeans(n_init="auto", n_clusters=n_clusters, random_state=43)
    km = k_Means.fit(temp)
    predict = km.predict(temp)
    
    temp['predicted'] = km.labels_
    colors = [
        'red', 'blue', 'orange', 'purple', 'green', 'pink', 
        'black', 'cadetblue', 'lightgreen','beige',
    ]
    temp['color'] = temp['predicted'].map(lambda p: colors[p])
    
    # Generate plots
    fig, ax = plt.subplots(figsize=(10,10), dpi=100)
    ax.minorticks_off()
    # ax = plt.axes()

    color = temp['color']
    for x in range(temp.shape[0]):
        curr = data.iloc[:,x]
        ax.plot(curr, label=curr.name, color=color[x])
        
    years = mdates.YearLocator(year_interval)  # every year
    months = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y') 
    
    ax = plt.gca()
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(yearsFmt)
    ax.autoscale_view()
    
    ax.set_xlabel("Years", fontsize=20)
    plt.xticks(fontsize=20)
    
    ax.set_ylabel(y_label, fontsize=20)
    plt.yticks(fontsize=20)
    
    # ax.set_title("{}: {} clusters".format(analyte_name, n_clusters), fontsize=20)
    ax.set_title(f"{analyte_name}: {n_clusters} clusters", fontsize=20)
    
    # Return clusters (if specified)
    if(return_clusters):
        color_df = pd.DataFrame(temp['color'])
        color_df['STATION_ID'] = color_df.index
        
        if data_pylenm_dm.get_construction_data() is None:
            cluster_logger.error("You need to set the GPS data first using the setConstructionData function.")
            return None
        else:
            gps_color = pd.merge(
                data_pylenm_dm.get_construction_data(), 
                color_df, 
                on=['STATION_ID'],
            )
            return gps_color
    
    # if(return_data):
    #     return color, temp