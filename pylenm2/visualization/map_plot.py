from ipywidgets import HTML
from ipyleaflet import (Map, basemaps, Marker,FullScreenControl, 
                        Popup, AwesomeIcon) 


def plot_coordinates_to_map(
        # self, 
        gps_data, 
        center=[33.271459, -81.675873], 
        zoom=14,
    ) -> Map:
    """Plots the station locations on an interactive map given coordinates.

    Args:
        gps_data (pd.DataFrame): Data frame with the following column names: station_id, latitude, longitude, color. If the color column is not passed, the default color will be blue.
        center (list, optional): latitude and longitude coordinates to center the map view. Defaults to [33.271459, -81.675873].
        zoom (int, optional): value to determine the initial scale of the map. Defaults to 14.

    Returns:
        ipyleaflet.Map
    """

    # center = center
    # zoom = 14

    # Create the basemap
    m = Map(
        basemap=basemaps.Esri.WorldImagery, 
        center=center, 
        zoom=zoom,
    )

    m.add_control(FullScreenControl())
    for (index,row) in gps_data.iterrows():

        # Create icons for stations locations
        if('color' in gps_data.columns):
            icon = AwesomeIcon(
                name='tint',
                marker_color=row.loc['color'],
                icon_color='black',
                spin=False
            )
        else:
            icon = AwesomeIcon(
                name='tint',
                marker_color='blue',
                icon_color='black',
                spin=False
            )

        loc = [row.loc['LATITUDE'],row.loc['LONGITUDE']]
        station = HTML(value=row.loc['STATION_ID'])

        # Create the Popup
        popup = Popup(
            location=loc,
            child=station,
            close_button=True,
            auto_close=False,
            max_height=1,
            close_on_escape_key=False,
        )
        
        # Create the Marker and add the Popup to it
        marker = Marker(location=loc,
                        icon=icon,
                        draggable=False,
                        popup=popup,
                    )

        m.add_layer(marker)


        # marker.popup = popup

    return m

