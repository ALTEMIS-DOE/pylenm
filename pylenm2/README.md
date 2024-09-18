# Pylenm2

Pylenm2 is intended to be used as a replacement for Pylenm.
This version organizes the code in a more semantic way.


### Demo
The demo notebooks are shown in the [notebooks2](../notebooks2/) directory.

### Organization hirarchy of package is shown below:

`pylenm2`
- `data`
    1. `data_module.py`
        1. CLASS: PylenmDataModule
            1. __init__(...)
            2. __has_columns(...)
            3. __is_valid(...)
            4. is_valid_data(...)
            5. is_valid_construction_data(...)
            6. __set_units(...)
            7. get_unit(...)
            8. set_data(...)
            9. set_construction_data(...)
            10. set_jointData(...)
            11. jointData(...)
            12. is_set_jointData(...)
            13. get_data(...)
            14. get_construction_data(...)
    2. `fetchers.py`
        1. _get_individual_analyte_df(...)
        2. _getLagDate(...)
        3. getCleanData(...)
        4. getCommonDates(...)
        5. getJointData(...)
        6. get_analyte_details(...)
        7. get_data_summary(...)
        8. get_station_analytes(...)
    3. `filters.py`
        1. simplify_data(...)
        2. filter_by_column(...)
        3. filter_stations(...)
        4. query_data(...)
    4. `transformation.py`
        1. interpolate_station_data(...)
        2. interpolate_stations_by_analyte(...)
        3. _transform_time_series(...)
        4. add_dist_to_source(...)
        5. time_average_all_stations(...)

- `stats`
    1. cluster.py
        1. __cluster_data_OLD(...)
        2. cluster_data(...)
    2. gp.py
        1. get_Best_GP(...)
        2. fit_gp(...)
        3. interpolate_topo(...)
    3. metrics.py
        1. dist(...)
        2. mse(...)
    4. preprocess.py
        1. _custom_analyte_sort(...)
        2. get_MCL(...)
        3. remove_outliers(...)
        4. remove_outliers_lowess(...)
        5. _get_Best_Station(...)
        6. get_Best_Stations(...)

- `utils`
    1. constants.py
    2. custom_exceptions.py

- `visualization`
    1. `correlation.py`
        1. _format_ticks_and_labels(...)
        2. plot_corr_by_station(...)
        3. plot_all_corr_by_station(...)
        4. plot_corr_by_date_range(...)
        5. plot_corr_by_year(...)
    2. `heatmap.py`
        1. plot_correlation_heatmap(...)
        2. plot_all_correlation_heatmap(...)
    3. `map_plot.py`
        1. plot_coordinates_to_map(...)
    4. `pca.py`
        1. plot_PCA_by_date(...)
        2. plot_PCA_by_year(...)
        3. plot_PCA_by_station(...)
    5. `plots.py`
        1. _plotUpperHalf(...)
        2. plot_data(...)
        3. plot_all_data(...)
        4. plot_MCL(...)
        5. _plot_data_xOutliers_fit(...)
        6. plot_data_rollAvg(...)
        7. plot_data_lowess(...)
    6. `timeseries.py`
        1. plot_all_time_series_simple(...)
        2. plot_all_time_series(...)