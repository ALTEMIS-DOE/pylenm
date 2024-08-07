def plot_corr_by_well(self, well_name, analytes, remove_outliers=True, z_threshold=4, interpolate=False, frequency='2W', save_dir='plot_correlation', log_transform=False, fontsize=20, returnData=False, remove=[], no_log=None):
    """Plots the correlations with the physical plots as well as the correlations of the important analytes over time for a specified well.

    Args:
        well_name (str): name of the well to be processed
        analytes (list): list of analyte names to use
        remove_outliers (bool, optional): choose whether or to remove the outliers. Defaults to True.
        z_threshold (int, optional): z_score threshold to eliminate outliers. Defaults to 4.
        interpolate (bool, optional): choose whether or to interpolate the data. Defaults to False.
        frequency (str, optional): {‘D’, ‘W’, ‘M’, ‘Y’} frequency to interpolate. Note: See https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html for valid frequency inputs. (e.g. ‘W’ = every week, ‘D ’= every day, ‘2W’ = every 2 weeks). Defaults to '2W'.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_correlation'.
        log_transform (bool, optional): flag for log base 10 transformation. Defaults to False.
        fontsize (int, optional): font size. Defaults to 20.
        returnData (bool, optional): flag to return data used to perfrom correlation analysis. Defaults to False.
        remove (list, optional): wells to remove. Defaults to [].
        no_log (list, optional): list of column names to not apply log transformation to. Defaults to None.

    Returns:
        None
    """
    data = self.data
    query = data[data.STATION_ID == well_name]
    a = list(np.unique(query.ANALYTE_NAME.values))# get all analytes from dataset
    for value in analytes:
        if((value in a)==False):
            return 'ERROR: No analyte named "{}" in data.'.format(value)
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
    if(interpolate):
        piv = self.interpolate_well_data(well_name, analytes, frequency=frequency)
        file_extension = '_interpolated_' + frequency
        title = well_name + '_correlation - interpolated every ' + frequency
    else:
        file_extension = '_correlation'
        title = well_name + '_correlation'
    samples = piv.shape[0]
    if(samples < 5):
        if(interpolate):
            return 'ERROR: {} does not have enough samples to plot.\n Try a different interpolation frequency'.format(well_name)
        return 'ERROR: {} does not have enough samples to plot.'.format(well_name)
    else:
        # scaler = StandardScaler()
        # pivScaled = scaler.fit_transform(piv)
        # pivScaled = pd.DataFrame(pivScaled, columns=piv.columns)
        # pivScaled.index = piv.index
        # piv = pivScaled
        
        if(log_transform):
            piv[piv <= 0] = 0.00000001
            temp = piv.copy()
            piv = np.log10(piv)
            if(no_log !=None):
                for col in no_log:
                    piv[col] = temp[col]

        # Remove outliers
        if(remove_outliers):
            piv = self.remove_outliers(piv, z_threshold=z_threshold)
        samples = piv.shape[0]

        idx = piv.index.date
        dates = [dates.strftime('%Y-%m-%d') for dates in idx]
        remaining = [i for i in dates if i not in remove]
        piv = piv.loc[remaining]
        
        sns.set_style("white", {"axes.facecolor": "0.95"})
        g = sns.PairGrid(piv, aspect=1.2, diag_sharey=False, despine=False)
        g.fig.suptitle(title, fontweight='bold', y=1.08, fontsize=25)
        g.map_lower(sns.regplot, lowess=True, ci=False, line_kws={'color': 'red', 'lw': 3},
                                                        scatter_kws={'color': 'black', 's': 20})
        g.map_diag(sns.distplot, kde_kws={'color': 'black', 'lw': 3}, hist_kws={'histtype': 'bar', 'lw': 2, 'edgecolor': 'k', 'facecolor':'grey'})
        g.map_upper(self.__plotUpperHalf)
        for ax in g.axes.flat:
            ax.tick_params("y", labelrotation=0, labelsize=fontsize)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize=fontsize)
            ax.set_xlabel(ax.get_xlabel(), fontsize=fontsize, fontweight='bold') #HERE
            ax.set_ylabel(ax.get_ylabel(), fontsize=fontsize,fontweight='bold')
            
        g.fig.subplots_adjust(wspace=0.3, hspace=0.3)
        ax = plt.gca()

        props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
        ax.text(1.3, 6.2, 'Start date:  {}\nEnd date:    {}\n\nOriginal samples:     {}\nSamples used:     {}'.format(piv.index[0].date(), piv.index[-1].date(), totalSamples, samples), transform=ax.transAxes, fontsize=20, fontweight='bold', verticalalignment='bottom', bbox=props)
        # Add titles to the diagonal axes/subplots
        for ax, col in zip(np.diag(g.axes), piv.columns):
            ax.set_title(col, y=0.82, fontsize=15)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        g.fig.savefig(save_dir + '/' + well_name + file_extension + '.png', bbox_inches="tight")
        if(returnData):
            return piv
        

def plot_all_corr_by_well(self, analytes, remove_outliers=True, z_threshold=4, interpolate=False, frequency='2W', save_dir='plot_correlation', log_transform=False, fontsize=20):
    """Plots the correlations with the physical plots as well as the important analytes over time for each well in the dataset.

    Args:
        analytes (list): list of analyte names to use
        remove_outliers (bool, optional): choose whether or to remove the outliers. Defaults to True.
        z_threshold (int, optional): z_score threshold to eliminate outliers. Defaults to 4.
        interpolate (bool, optional): choose whether or to interpolate the data. Defaults to False.
        frequency (str, optional): {‘D’, ‘W’, ‘M’, ‘Y’} frequency to interpolate. Note: See https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html for valid frequency inputs. (e.g. ‘W’ = every week, ‘D ’= every day, ‘2W’ = every 2 weeks). Defaults to '2W'.
        save_dir (str, optional): name of the directory you want to save the plot to. Defaults to 'plot_correlation'.
        log_transform (bool, optional): flag for log base 10 transformation. Defaults to False.
        fontsize (int, optional): font size. Defaults to 20.
    """
    data = self.data
    wells = np.array(data.STATION_ID.values)
    wells = np.unique(wells)
    for well in wells:
        self.plot_corr_by_well(well_name=well, analytes=analytes,remove_outliers=remove_outliers, z_threshold=z_threshold, interpolate=interpolate, frequency=frequency, save_dir=save_dir, log_transform=log_transform, fontsize=fontsize)
    
def plot_corr_by_date_range(self, date, analytes, lag=0, min_samples=10, save_dir='plot_corr_by_date', log_transform=False, fontsize=20, returnData=False, no_log=None):
    """Plots the correlations with the physical plots as well as the correlations of the important analytes for ALL the wells on a specified date or range of dates if a lag greater than 0 is specifed.

    Args:
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
    if(lag==0):
        data = self.data
        data = self.simplify_data(data=data)
        query = data[data.COLLECTION_DATE == date]
        a = list(np.unique(query.ANALYTE_NAME.values))# get all analytes from dataset
        for value in analytes:
            if((value in a)==False):
                return 'ERROR: No analyte named "{}" in data.'.format(value)
        analytes = sorted(analytes)
        query = query.loc[query.ANALYTE_NAME.isin(analytes)]
        if(query.shape[0] == 0):
            return 'ERROR: {} has no data for all of the analytes.'.format(date)
        samples = query[['COLLECTION_DATE', 'STATION_ID', 'ANALYTE_NAME']].duplicated().value_counts()[0]
        if(samples < min_samples):
            return 'ERROR: {} does not have at least {} samples.'.format(date, min_samples)
        else:
            piv = query.reset_index().pivot_table(index = 'STATION_ID', columns='ANALYTE_NAME', values='RESULT',aggfunc=np.mean)
            # return piv
    else:
        # If the data has already been calculated with the lag specified, retrieve it
        if(self.jointData_is_set(lag=lag)==True): 
            data = self.__jointData[0]
        # Otherwise, calculate it
        else:
            data = self.getJointData(analytes, lag=lag)
            self.__set_jointData(data=data, lag=lag)
        # get new range based on the lag and create the pivor table to be able to do the correlation
        dateStart, dateEnd = self.__getLagDate(date, lagDays=lag)
        dateRange_key = str(dateStart.date()) + " - " + str(dateEnd.date())
        piv = pd.DataFrame(data.loc[dateRange_key]).unstack().T
        piv.index = piv.index.droplevel()
        piv = pd.DataFrame(piv).dropna(axis=0, how='all')
        num_NaNs = int(piv.isnull().sum().sum())
        samples = (piv.shape[0]*piv.shape[1])-num_NaNs
        for col in piv.columns:
            piv[col] = piv[col].astype('float64', errors = 'raise')
        if(lag>0):
            date = dateRange_key
        # return piv
    title = date + '_correlation'
    # scaler = StandardScaler()
    # pivScaled = scaler.fit_transform(piv)
    # pivScaled = pd.DataFrame(pivScaled, columns=piv.columns)
    # pivScaled.index = piv.index
    # piv = pivScaled

    if(log_transform):
        piv[piv <= 0] = 0.00000001
        temp = piv.copy()
        piv = np.log10(piv)
        if(no_log !=None):
            for col in no_log:
                piv[col] = temp[col]

    sns.set_style("white", {"axes.facecolor": "0.95"})
    g = sns.PairGrid(piv, aspect=1.2, diag_sharey=False, despine=False)
    g.fig.suptitle(title, fontweight='bold', y=1.08, fontsize=25)
    g.map_lower(sns.regplot, lowess=True, ci=False, line_kws={'color': 'red', 'lw': 3},
                                                    scatter_kws={'color': 'black', 's': 20})
    g.map_diag(sns.distplot, kde_kws={'color': 'black', 'lw': 3}, hist_kws={'histtype': 'bar', 'lw': 2, 'edgecolor': 'k', 'facecolor':'grey'})
    g.map_upper(self.__plotUpperHalf)
    for ax in g.axes.flat:
            ax.tick_params("y", labelrotation=0, labelsize=fontsize)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize=fontsize)
            ax.set_xlabel(ax.get_xlabel(), fontsize=fontsize, fontweight='bold') #HERE
            ax.set_ylabel(ax.get_ylabel(), fontsize=fontsize,fontweight='bold')
    g.fig.subplots_adjust(wspace=0.3, hspace=0.3)
    ax = plt.gca()

    props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
    ax.text(1.3, 3, 'Date:  {}\n\nWells:     {}\nSamples used:     {}'.format(date, piv.shape[0] ,samples), transform=ax.transAxes, fontsize=20, fontweight='bold', verticalalignment='bottom', bbox=props)
    # Add titles to the diagonal axes/subplots
    for ax, col in zip(np.diag(g.axes), piv.columns):
        ax.set_title(col, y=0.82, fontsize=15)
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    g.fig.savefig(save_dir + '/' + date + '.png', bbox_inches="tight")
    if(returnData):
        return piv


def plot_corr_by_year(self, year, analytes, remove_outliers=True, z_threshold=4, min_samples=10, save_dir='plot_corr_by_year', log_transform=False, fontsize=20, returnData=False, no_log=None):
    """Plots the correlations with the physical plots as well as the correlations of the important analytes for ALL the wells in specified year.

    Args:
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
    data = self.data
    query = data
    query = self.simplify_data(data=query)
    query.COLLECTION_DATE = pd.to_datetime(query.COLLECTION_DATE)
    query = query[query.COLLECTION_DATE.dt.year == year]
    a = list(np.unique(query.ANALYTE_NAME.values))# get all analytes from dataset
    for value in analytes:
        if((value in a)==False):
            return 'ERROR: No analyte named "{}" in data.'.format(value)
    analytes = sorted(analytes)
    query = query.loc[query.ANALYTE_NAME.isin(analytes)]
    if(query.shape[0] == 0):
        return 'ERROR: {} has no data for the 6 analytes.'.format(year)
    samples = query[['COLLECTION_DATE', 'STATION_ID', 'ANALYTE_NAME']].duplicated().value_counts()[0]
    if(samples < min_samples):
        return 'ERROR: {} does not have at least {} samples.'.format(date, min_samples)
    else:
        piv = query.reset_index().pivot_table(index = 'STATION_ID', columns='ANALYTE_NAME', values='RESULT',aggfunc=np.mean)
        # return piv
        # Remove outliers
        if(remove_outliers):
            piv = self.remove_outliers(piv, z_threshold=z_threshold)
        samples = piv.shape[0] * piv.shape[1]

        title = str(year) + '_correlation'
        # scaler = StandardScaler()
        # pivScaled = scaler.fit_transform(piv)
        # pivScaled = pd.DataFrame(pivScaled, columns=piv.columns)
        # pivScaled.index = piv.index
        # piv = pivScaled

        if(log_transform):
            piv[piv <= 0] = 0.00000001
            temp = piv.copy()
            piv = np.log10(piv)
            if(no_log !=None):
                for col in no_log:
                    piv[col] = temp[col]

        sns.set_style("white", {"axes.facecolor": "0.95"})
        g = sns.PairGrid(piv, aspect=1.2, diag_sharey=False, despine=False)
        g.fig.suptitle(title, fontweight='bold', y=1.08, fontsize=25)
        g.map_lower(sns.regplot, lowess=True, ci=False, line_kws={'color': 'red', 'lw': 3},
                                                        scatter_kws={'color': 'black', 's': 20})
        g.map_diag(sns.distplot, kde_kws={'color': 'black', 'lw': 3}, hist_kws={'histtype': 'bar', 'lw': 2, 'edgecolor': 'k', 'facecolor':'grey'})
        g.map_upper(self.__plotUpperHalf)
        for ax in g.axes.flat:
            ax.tick_params("y", labelrotation=0, labelsize=fontsize)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, fontsize=fontsize)
            ax.set_xlabel(ax.get_xlabel(), fontsize=fontsize, fontweight='bold') #HERE
            ax.set_ylabel(ax.get_ylabel(), fontsize=fontsize,fontweight='bold')
        g.fig.subplots_adjust(wspace=0.3, hspace=0.3)
        ax = plt.gca()

        props = dict(boxstyle='round', facecolor='grey', alpha=0.15)
        ax.text(1.3, 3, 'Date:  {}\n\nSamples used:     {}'.format(year, samples), transform=ax.transAxes, fontsize=20, fontweight='bold', verticalalignment='bottom', bbox=props)
        # Add titles to the diagonal axes/subplots
        for ax, col in zip(np.diag(g.axes), piv.columns):
            ax.set_title(col, y=0.82, fontsize=15)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        g.fig.savefig(save_dir + '/' + str(year) + '.png', bbox_inches="tight")
        if(returnData):
            return piv