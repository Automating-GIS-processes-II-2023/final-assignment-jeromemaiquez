"""
Functions for Urban Assessment Tool.

List of included functions:
    - cleanup_gdf()
    - ...

Usage:
    ./urban_indicators_scripts.py

Author:
    Jerome Maiquez - xx.xx.2024
"""

def cleanup_gdf(gdf, na_cutoff_percent, target_geom_type):
    """
    Function for cleaning up geodataframe, performing the following:
        - drops columns with NaN counts above cutoff
        - only retains records with target geometry type
    
    Parameters
    ----------
    gdf: <geopandas.geodataframe.GeoDataFrame>
        GeoDataFrame to undergo cleaning.
    na_cutoff_percent: <float> (from 0.0 to 1.0)
        Minimum percent of rows with NaN values in a column
        for it to be retained in the GeoDataFrame.
    target_geom_type: <shapely.Geometry>
        Geometry type to be retained.

    Returns
    -------
    <geopandas.geodataframe.GeoDataFrame>
        GeoDataFrame with only relevant columns and correct geometry type.
    """

    # import geopandas
    import geopandas as gpd
    
    # assertions to weed out wrong input types
    assert type(gdf) == gpd.geodataframe.geodataFrame, "df must be geodataframe"
    assert 0.0 <= na_cutoff_percent <= 1.0, "NaN cutoff percent must be between 0 and 1"
    assert target_geom_type in ["Point", "LineString", "Polygon"], "target geom type must be valid type of geometry"