import os
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
from matplotlib import pyplot as plt
from sklearn.preprocessing import StandardScaler

from pylenm2.stats import preprocess
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

heatmap_logger = logger_config.setup_logging(
    module_name=__name__,
    level=logging.INFO,
    # level=logging.DEBUG,
    logfile_dir=c.LOGFILE_DIR,
)


def plot_correlation_heatmap(
        # self, 
        data_pylenm_dm, 
        station_name, 
        show_symmetry=True, 
        color=True, 
        save_dir='plot_correlation_heatmap',
    ):
    """ Plots a heatmap of the correlations of the important analytes over time for a specified station.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        station_name (str): name of the station to be processed
        show_symmetry (bool, optional): choose whether or not the heatmap should show the same information twice over the diagonal. Defaults to True.
        color (bool, optional): choose whether or not the plot should be in color or in greyscale. Defaults to True.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_correlation_heatmap'.

    Returns:
        None
    """
    
    data = data_pylenm_dm.data
    query = data[data.STATION_ID == station_name]
    
    a = list(np.unique(query.ANALYTE_NAME.values))
    b = ['TRITIUM','IODINE-129','SPECIFIC CONDUCTANCE', 'PH','URANIUM-238', 'DEPTH_TO_WATER']
    
    analytes = preprocess._custom_analyte_sort(list(set(a) and set(b)))
    query = query.loc[query.ANALYTE_NAME.isin(analytes)]
    analytes = preprocess._custom_analyte_sort(np.unique(query.ANALYTE_NAME.values))
    
    x = query[['COLLECTION_DATE', 'ANALYTE_NAME']]
    unique = ~x.duplicated()
    query = query[unique]
    
    piv = query.reset_index().pivot(
        index='COLLECTION_DATE', 
        columns='ANALYTE_NAME', 
        values='RESULT',
    )
    piv = piv[analytes]
    totalSamples = piv.shape[0]
    piv = piv.dropna()
    samples = piv.shape[0]
    
    if(samples < 5):
        heatmap_logger.debug(f"ERROR: {station_name} does not have enough samples to plot.")
        return None
    
    else:
        scaler = StandardScaler()
        pivScaled = scaler.fit_transform(piv)
        pivScaled = pd.DataFrame(pivScaled, columns=piv.columns)
        pivScaled.index = piv.index
        
        corr_matrix = pivScaled.corr()
        
        if(show_symmetry):
            mask = None
        else:
            mask = np.triu(corr_matrix)
        if(color):
            cmap = 'RdBu'
        else:
            cmap = 'binary'
        
        fig, ax = plt.subplots(figsize=(8,6))
        ax.set_title(station_name + '_correlation', fontweight='bold')
        ttl = ax.title
        ttl.set_position([.5, 1.05])
        
        props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
        ax.text(
            1.3, 
            1.05, 
            # 'Start date: {}\nEnd date: {}\n\nSamples: {} of {}'.format(piv.index[0], piv.index[-1], samples, totalSamples), 
            f"Start date: {piv.index[0]}\nEnd date: {piv.index[-1]}\n\nSamples: {samples} of {totalSamples}", 
            transform=ax.transAxes, 
            fontsize=15, 
            fontweight='bold', 
            verticalalignment='bottom', 
            bbox=props,
        )
        
        ax = sns.heatmap(
            corr_matrix,
            ax=ax,
            mask=mask,
            vmin=-1, vmax=1,
            xticklabels=corr_matrix.columns,
            yticklabels=corr_matrix.columns,
            cmap=cmap,
            annot=True,
            linewidths=1,
            cbar_kws={'orientation': 'vertical'},
        )

        # Save plot
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        fig.savefig(
            # save_dir + '/' + station_name + '_correlation.png', 
            f"{save_dir}/{station_name}_correlation.png", 
            bbox_inches="tight",
        )

        # Returning a random non-None value to mark successful execution
        return 200


def plot_all_correlation_heatmap(
        # self, 
        data_pylenm_dm, 
        show_symmetry=True, 
        color=True, 
        save_dir='plot_correlation_heatmap',
    ):
    """Plots a heatmap of the correlations of the important analytes over time for each station in the dataset.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        show_symmetry (bool, optional): choose whether or not the heatmap should show the same information twice over the diagonal. Defaults to True.
        color (bool, optional): choose whether or not the plot should be in color or in greyscale. Defaults to True.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_correlation_heatmap'.
    """
    
    data = data_pylenm_dm.data
    
    stations = np.array(data.STATION_ID.values)
    stations = np.unique(stations)

    successful_plots = 0
    failed_plots = 0
    
    for station in tqdm(stations, desc="Stations"):
        plot_status = plot_correlation_heatmap(
            data_pylenm_dm=data_pylenm_dm,
            station_name=station,
            show_symmetry=show_symmetry,
            color=color,
            save_dir=save_dir,
        )

        if plot_status is None:
            failed_plots += 1
        else:
            successful_plots += 1
        
    heatmap_logger.info(f"Successful plots: {successful_plots}")
    heatmap_logger.info(f"Failed plots: {failed_plots}")