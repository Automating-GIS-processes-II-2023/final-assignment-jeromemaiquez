"""
Functions for Urban Assessment Tool.

List of included functions:
    - filter_elem_type()
    - drop_nan_cols()
    - get_buffer_from_place()
    - parse_street_lanes()

Usage:
    ./urban_indicators_scripts.py

Author:
    Jerome Maiquez - xx.xx.2024
"""

def filter_elem_type(gdf, elem_type):
    """
    Function for filtering an OSM gdf by element type.
    Element type could be either node, way, or relation.

    Parameters
    ----------
    gdf: <geopandas.geodataframe.GeoDataFrame>
        GeoDataFrame to undergo subsetting.
    elem_type: <str> or <list> [ "node" | "way" | "relation" ]
        OSM element type(s) to retain in GeoDataFrame.

    Returns
    -------
    <geopandas.geodataframe.GeoDataFrame>
        Same geoDataFrame, whose rows are only those with matching element type. 
    """

    # import geopandas
    import geopandas as gpd

    # assertions to weed out wrong input types
    assert type(gdf) == gpd.geodataframe.GeoDataFrame, "df must be geodataframe"
    assert type(elem_type) in (str, list), "elem type must be a str or list of strs"

    # Reset index of gdf (which is multi-indexed by default)
    gdf.reset_index(inplace=True)

    # Select only the rows with the specified element type
    if type(elem_type) == str:
        assert elem_type in ("node", "way", "relation"), "elem type must be one of: node, way, or relation"
        return gdf.loc[gdf["element_type"] == elem_type].copy()

    elif type(elem_type) == list:
        assert all([i in ("node", "way", "relation") for i in elem_type]), "elem types must be: node, way, and/or relation"
        return gdf.loc[gdf["element_type"].isin(elem_type)].copy()


def drop_nan_cols(gdf, na_cutoff_percent):
    """
    Function for subsetting only columns with enough data.
    Drops columns with NaN counts above a supplied cutoff value.
    
    Parameters
    ----------
    gdf: <geopandas.geodataframe.GeoDataFrame>
        GeoDataFrame to undergo subsetting.
    na_cutoff_percent: <float> (from 0.0 to 1.0)
        Minimum percent of rows with NaN values in a column
        for it to be retained in the GeoDataFrame.

    Returns
    -------
    <geopandas.geodataframe.GeoDataFrame>
        Same geoDataFrame, but with only relevant columns retained.
    """

    # import geopandas
    import geopandas as gpd
    
    # assertions to weed out wrong input types
    assert type(gdf) == gpd.geodataframe.GeoDataFrame, "df must be geodataframe"
    assert 0.0 <= na_cutoff_percent <= 1.0, "NaN cutoff percent must be between 0 and 1"
    # assert target_geom_type in ["Point", "LineString", "Polygon"], "target geom type must be valid type of geometry"

    # Define NaN count cutoff value
    na_cutoff = len(gdf) * na_cutoff_percent

    # Create empty list for columns to be dropped
    drop_cols = []

    # For each column in gdf...
    for col in gdf:
        # If this column has more NaN rows than the cutoff value...
        if gdf[col].isna().sum() >= na_cutoff:
            # Add column to list of columns to be dropped
            drop_cols.append(col)
    
    # Drop all columns in drop_cols list
    return gdf.drop(columns=drop_cols)


def get_buffer_from_place(query, dist):
    """
    Function for generating a circle of radius `dist` around
    a point corresponding to a geocoded place name or address.

    Usage
    -----
    `dist` must be provided in meters.

    Parameters
    ----------
    query: <str>
        Place name or address to be geocoded.
    dist: <int> or <float>
        Distance (in meters) of buffer around point of interest.

    Returns
    -------
    <geopandas.geodataframe.GeoDataFrame>
        Single-row geodataframe containing place name and polygon geometry.
        Returned GeoDataFrame CRS is EPSG: 4326.
    """

    # import packages
    import geopandas as gpd
    import osmnx as ox
    import shapely as shp

    # assertions to weed out wrong input data types
    assert type(query) == str, "Place name/address must be a str"
    assert type(dist) in (int, float), "Distance must be numeric"

    # Geocode place name as point
    place = ox.geocode(query)

    # Parse place name out of query
    place_name = query.split(", ")[0] if ", " in query else query

    # Create geodataframe with two columns:
    buffer = gpd.GeoDataFrame(
        {
            # 1) name of place
            "name": [place_name],
            # 2) point geometry for geocoded lat-lon tuple
            # tuple must be reversed so that x comes first
            "geometry": [shp.Point(reversed(place))]
        },
        # Set CRS as WGS 84 first, to be changed later
        crs="EPSG:4326"
    ).to_crs("EPSG:3857") # <-- change CRS here

    # Set geometry to buffer of radius `dist` around point
    buffer["geometry"] = buffer["geometry"].buffer(dist)

    # Reproject back to WGS84 (for OSMnx), then return buffered GDF
    return buffer.to_crs("EPSG:4326")

def parse_street_lanes(streets_gdf):
    """
    Function to clean up and parse the `lanes` column
    of a street GeoDataFrame (edges of an OSM street graph).

    Parameters
    ----------
    streets_gdf: <geopandas.geodataframe.GeoDataFrame>
        Street network GDF whose `lanes` column is to be cleaned.
    
    Returns
    -------
    <geopandas.geodataframe.GeoDataFrame>
        Same GeoDataFrame, but with a parsed `lanes` column.
    """

    # import packages
    import geopandas as gpd

    # assertions to weed out wrong input types
    assert type(streets_gdf) == gpd.geodataframe.GeoDataFrame, "gdf must be geodataframe"
    assert all(x in streets_gdf.columns for x in ["osmid", "lanes", "geometry"]), "streets gdf must be from OSM"

    # Set all NaN values in `lanes` column to 1
    streets_gdf.loc[streets_gdf["lanes"].isna(), "lanes"] = 1

    # If value in `lanes` column is a list, set to its max value
    streets_gdf["lanes"] = streets_gdf["lanes"].apply(
        lambda x: max(x) if isinstance(x, list) == True else x
    )

    # Convert lanes column to type int
    streets_gdf["lanes"] = streets_gdf["lanes"].astype(int)

    return streets_gdf