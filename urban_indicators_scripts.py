"""
Functions for Urban Assessment Tool.

List of included functions:
    - filter_elem_type()
    - drop_nan_cols()

Usage:
    ./urban_indicators_scripts.py

Author:
    Jerome Maiquez - xx.xx.2024
"""

def filter_elem_type(gdf, elem_type):
    """
    Function for filtering an OSM gdf by element type.
    Element type could be either node, way, or relation.
    """

    # import geopandas
    import geopandas as gpd

    # assertions to weed out wrong input types
    assert type(gdf) == gpd.geodataframe.GeoDataFrame, "df must be geodataframe"
    assert elem_type in ("node", "way", "relation"), "elem type must be one of: node, way, or relation"

    # Reset index of gdf (which is multi-indexed by default)
    gdf.reset_index(inplace=True)

    # Select only the rows with the specified element type
    return gdf.loc[gdf["element_type"] == elem_type].copy()


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
        GeoDataFrame with only relevant columns retained.
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