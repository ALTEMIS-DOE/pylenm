def interpolate_well_data(self, well_name, analytes, frequency='2W'):
    """Resamples the data based on the frequency specified and interpolates the values of the analytes.

    Args:
        well_name (str): name of the well to be processed.
        analytes (list): list of analyte names to use
        frequency (str, optional): {‘D’, ‘W’, ‘M’, ‘Y’} frequency to interpolate. See https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html for valid frequency inputs. (e.g. ‘W’ = every week, ‘D ’= every day, ‘2W’ = every 2 weeks). Defaults to '2W'.

    Returns:
        pd.DataFrame
    """
    data = self.data
    inter_series = {}
    query = data[data.STATION_ID == well_name]
    for analyte in analytes:
        series = query[query.ANALYTE_NAME == analyte]
        series = (series[['COLLECTION_DATE', 'RESULT']])
        series.COLLECTION_DATE = pd.to_datetime(series.COLLECTION_DATE)
        series.index = series.COLLECTION_DATE
        original_dates = series.index
        series = series.drop('COLLECTION_DATE', axis=1)
        series = series.rename({'RESULT': analyte}, axis=1)
        upsampled = series.resample(frequency).mean()
        interpolated = upsampled.interpolate(method='linear', order=2)
        inter_series[analyte] = interpolated
    join = inter_series[analytes[0]]
    join = join.drop(analytes[0], axis=1)
    for analyte in analytes:
        join = join.join(inter_series[analyte])
    join = join.dropna()
    return join


def interpolate_wells_by_analyte(self, analyte, frequency='2W', rm_outliers=True, z_threshold=3):
    """Resamples analyte data based on the frequency specified and interpolates the values in between. NaN values are replaced with the average value per well.

    Args:
        analyte (_type_): analyte name for interpolation of all present wells.
        frequency (str, optional): {‘D’, ‘W’, ‘M’, ‘Y’} frequency to interpolate. See https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html for valid frequency inputs. (e.g. ‘W’ = every week, ‘D ’= every day, ‘2W’ = every 2 weeks). Defaults to '2W'.
        rm_outliers (bool, optional): flag to remove outliers in the data. Defaults to True.
        z_threshold (int, optional): z_score threshold to eliminate outliers. Defaults to 3.

    Returns:
        pd.DataFrame: interpolated data
    """
    data = self.data
    df_t, dates = self.__transform_time_series( 
                                                analytes=[analyte], 
                                                resample=frequency, 
                                                rm_outliers=True, 
                                                z_threshold=z_threshold)
    res_interp = self.__get_individual_analyte_df(data=df_t, dates=dates, analyte=analyte)
    res_interp = res_interp.dropna(axis=1, how='all')
    return res_interp


# IN THE WORKS
def __transform_time_series(self, analytes=[], resample='2W', rm_outliers=False, z_threshold=4):
    data = self.data
    def transform_time_series_by_analyte(data, analyte_name):
        wells_analyte = np.unique(data[data.ANALYTE_NAME == analyte_name].STATION_ID)
        condensed = data[data.ANALYTE_NAME == analyte_name].groupby(['STATION_ID','COLLECTION_DATE']).mean()
        analyte_df_resample = pd.DataFrame(index=wells_analyte, columns=t)
        analyte_df_resample.sort_index(inplace=True)
        for well in wells_analyte:
            for date in condensed.loc[well].index:
                analyte_df_resample.at[well, date] = condensed.loc[well,date].RESULT
        analyte_df_resample = analyte_df_resample.astype('float').T
        analyte_df_resample = analyte_df_resample.interpolate(method='linear')
        return analyte_df_resample

    all_dates = np.unique(data.COLLECTION_DATE)
    # Create array of equally spaced dates
    start = pd.Timestamp(all_dates.min())
    end = pd.Timestamp(all_dates.max())
    delta = end - start
    t = np.linspace(start.value, end.value, delta.days)
    t = pd.to_datetime(t)
    t = pd.Series(t)
    t = t.apply(lambda x: x.replace(minute=0, hour=0, second=0, microsecond=0, nanosecond=0))

    cutoff_dates = []
    # Save each analyte data
    analyte_data = []
    for analyte in analytes:
        ana_data = transform_time_series_by_analyte(data, analyte)
        if(rm_outliers):
            col_num = ana_data.shape[1]
            for col in range(col_num):
                ana_data.iloc[:,col] = self.remove_outliers(ana_data.iloc[:,col].dropna(), z_threshold=z_threshold)
            ana_data = ana_data.interpolate(method='linear')
        ana_data.index = pd.to_datetime(ana_data.index)
        # Resample
        ana_data_resample = ana_data.resample(resample).mean()
        # Save data
        analyte_data.append(ana_data_resample)
        # Determine cuttoff point for number of NaNs in dataset
        passes_limit = []
        for date in ana_data_resample.index:
            limit = 0.7 * ana_data_resample.shape[1]
            curr = ana_data_resample.isna().loc[date,:].value_counts()
            if('False' in str(curr)):
                curr_total = ana_data_resample.isna().loc[date,:].value_counts()[0]
                if curr_total > limit:
                    passes_limit.append(date)
        passes_limit = pd.to_datetime(passes_limit)
        cutoff_dates.append(passes_limit.min())
    start_index = pd.Series(cutoff_dates).max()

    # Get list of shared wells amongst all the listed analytes
    combined_well_list = []
    for x in range(len(analytes)):
        combined_well_list = combined_well_list + list(analyte_data[x].columns)
    combined_count = pd.Series(combined_well_list).value_counts()
    shared_wells = list(combined_count[list(pd.Series(combined_well_list).value_counts()==len(analytes))].index)

    # Vectorize data
    vectorized_df = pd.DataFrame(columns=analytes, index = shared_wells)

    for analyte, num in zip(analytes, range(len(analytes))):
        for well in shared_wells:
            analyte_data_full = analyte_data[num][well].fillna(analyte_data[num][well].mean())
            vectorized_df.at[well, analyte] = analyte_data_full[start_index:].values

    dates = ana_data_resample[start_index:].index
    return vectorized_df, dates


def add_dist_to_source(self, XX, source_coordinate=[436642.70,3681927.09], col_name='dist_to_source'):
    """adds column to data with the distance of a record to the source coordinate

    Args:
        XX (pd.DataFrame): data with coordinate information
        source_coordinate (list, optional): source coordinate. Defaults to [436642.70,3681927.09].
        col_name (str, optional): name to assign new column. Defaults to 'dist_to_source'.

    Returns:
        pd.DataFrame: returns original data with additional column with the distance.
    """
    x1,y1 = source_coordinate
    distances = []
    for i in range(XX.shape[0]):
        x2,y2 = XX.iloc[i][0], XX.iloc[i][1]
        distances.append(self.dist([x1,y1],[x2,y2]))
    XX[col_name] = distances
    return XX