"""
Functions for Urban Assessment Tool.

List of included functions:
    - filter_elem_type()
    - drop_nan_cols()
    - get_buffer_from_place()
    - parse_street_lanes()
    - plot_streets()
    - plot_heatmap()

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


def plot_streets(streets_gdf, ring_gdf, title=None, subtitle=None, width_factor=0.25, color_scheme="Greys_r"):
    """
    Function to plot street network GDF in a minimalist style.

    Parameters
    ----------
    streets_gdf: <geopandas.geodataframe.GeoDataFrame>
        Street network GeoDataFrame to be plotted.
    ring_gdf: <geopandas.geodataframe.GeoDataFrame>
        Buffer ring GeoDataFrame, whose bounds will be the basis for the axes.
    width_factor: <float>
        Scaling factor between number of lanes and line width in plot.
    color_scheme: <str> or <tuple>
        Matplotlib color map to use for edge and face colors.
    title: <str>
        Main title of the plot. If `None`, gets the query name.
    subtitle: <str>
        Subtitle of the plot. If `None`, displays "Road Network..."
    
    Returns
    -------
    <matplotlib.pyplot.figure>
        Figure containing the street network plot.
    <matplotlib.axes.Axes>
        Axes containing the street network plot and text.
    <tuple>
        Four-element tuple containing bounds of buffer ring.
    """
    
    # Import packages
    import geopandas as gpd
    import matplotlib.pyplot as plt
    from matplotlib import colormaps as cmaps

    # Assertions to weed out incorrect data types
    assert type(streets_gdf) == gpd.geodataframe.GeoDataFrame, "streets gdf must be geodataframe"
    assert all(x in streets_gdf.columns for x in ["osmid", "lanes", "geometry"]), "streets gdf must be from OSM"
    assert type(ring_gdf) == gpd.geodataframe.GeoDataFrame, "buffer gdf must be geodataframe"
    assert all(x in ring_gdf.columns for x in ["name", "geometry"]), "buffer gdf must have name & geometry"
    
    assert type(width_factor) == float, "width factor must be a float"
    assert 0.0 <= width_factor <= 1.0, "width factor must be between 0 and 1"
    assert type(color_scheme) in (str, tuple), "color_scheme must follow matplotlib color formats"

    if title != None:
        assert type(title) == str, "if supplied, title must be a str"
    if subtitle != None:
        assert type(subtitle) == str, "if supplied, subtitle must be a str"
    
    # Get min and max number of lanes
    min_lanes = streets_gdf["lanes"].min()
    max_lanes = streets_gdf["lanes"].max()

    # Assign line width per number of lanes
    linewidths = {lane: lane * width_factor for lane in range(min_lanes, max_lanes + 1)}

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 10))

    # Assign edge & face colors from specified color scheme
    cmap = cmaps[color_scheme].resampled(5)
    edgecolor = cmap(1.0)
    facecolor = cmap(0.0)
    ringcolor = cmap(0.8)

    # For each lane number in streets gdf...
    for lane, data in streets_gdf.groupby("lanes"):
        # Plot all streets with x no. of lanes...
        data.plot(
            color=edgecolor,
            ax=ax,
            label=lane,
            linewidth=linewidths[lane] # ...with the correct line width
        )
    
    # Get bounding box & radius of buffer ring...
    proj_ring_gdf = ring_gdf.to_crs(streets_gdf.crs)
    minx, miny, maxx, maxy = proj_ring_gdf.total_bounds
    ring_radius = round(proj_ring_gdf.minimum_bounding_radius()[0], 0)
    ax_pad = ring_radius * 0.10

    # ...to set x- and y-limits of axis (to ensure it is square)
    ax.set_xlim(left=minx - ax_pad, right=maxx + ax_pad)
    ax.set_ylim(bottom=miny - ax_pad, top=maxy + ax_pad)

    # Configure plot title and subtitle
    plot_title = ring_gdf.loc[0, "name"] if title == None else title
    plot_subtitle = f"Road Network (within a {ring_radius/1_000:.1f} km buffer)" if subtitle == None else subtitle

    # Set other axis properties
    ax.set_facecolor(facecolor)
    ax.set_title(
        plot_title,
        loc="left",
        color=edgecolor,

        fontsize=45.0,
        fontfamily="Garamond",
        fontstretch="expanded",

        y=0.09, # 8% of figure height (from bottom)
        x=0.03, # 3% of figure width (from left)
        pad=3.0
    )
    ax.text(
        x=minx-(ax_pad*0.25),
        y=miny-(ax_pad*(-0.25)),
        s=plot_subtitle,
        color=edgecolor,
        fontsize=20.0,
        fontfamily="Garamond",
        fontstretch="expanded"
    )
    ax.text(
        x=minx-(ax_pad*0.25),
        y=miny-(ax_pad*0.35),
        s="(c) OpenStreetMap Contributors",
        color=edgecolor,
        fontsize=10.0,
        fontfamily="Garamond",
        fontstretch="expanded"
    )

    # Add buffer ring to plot
    proj_ring_gdf.plot(
        edgecolor=ringcolor,
        facecolor="none",
        linewidth=0.75,
        linestyle="dotted",
        ax=ax
    )

    # Set colors of figure
    fig.set_facecolor(facecolor)
    fig.set_edgecolor(facecolor)

    # Remove axis labels for photo finish
    ax.set_axis_off()

    # Store buffer bounds in 1 variable
    ring_bounds = (minx, maxx, miny, maxy)

    # Show plot
    return fig, ax, ring_bounds


def plot_heatmap(gdf, weights, ax, bounds, bins=100, smoothing=1.0, alpha=0.4, color_scheme="Greys_r"):
    """
    Function that plots a heatmap of a given point GeoDataFrame.
    Must be run in the same cell as plot_streets(), for it assumes a given ax & bounds.

    Parameters
    ----------
    gdf: <geopandas.geodataframe.GeoDataFrame>
        Point GeoDataFrame to be plotted into a heatmap.
    weights: <str>
        Label of column in `gdf` whose values are to be used as heatmap weights.
    ax: <matplotlib.axes.Axes>
        Axes on which the heatmap will be plotted.
    bounds: <tuple>
        Four-element tuple containing min and max x- and y-coordinates for the heatmap.
    bins: <int>
        Number of pixels for the heatmap image on both dimensions.
    smoothing: <float>
        Smoothing factor for Gaussian filter on heatmap image.
    alpha: <float>
        Opacity of heatmap (0.0 = transparent, 1.0 = opaque).
    color_scheme: <str> or <tuple>
        Matplotlib color map to use for heatmap visualization.
    
    Returns
    -------
    <matplotlib.image>
        2D image of the generated heatmap.
    """

    # Import packages
    import numpy as np
    from scipy import ndimage
    import matplotlib.pyplot as plt
    import geopandas as gpd

    # Assertions to weed out incorrect data types
    assert type(gdf) == gpd.geodataframe.GeoDataFrame, "input gdf must be geodataframe"
    assert all([x in ["Point", "Polygon"] for x in gdf.geom_type]), "geometries in geodataframe must be point"
    assert weights in gdf.columns.values, "weights must be a column in gdf"

    assert type(smoothing) == float, "smoothing factor must be a float"
    assert type(bins) == int, "number of bins must be an int"
    assert type(color_scheme) in (str, tuple), "color_scheme must follow matplotlib color formats"

    # If geometries are polygons, convert to centroids
    if all([x == "Polygon" for x in gdf.geom_type]):
        gdf = gdf.set_geometry(gdf.centroid)

    # Get x- and y-coordinates of points
    coords = gdf.geometry.get_coordinates()

    # Create 2D weighted histogram
    heatmap, _, _ = np.histogram2d(
        x = coords.y,
        y = coords.x,
        bins =  bins,
        range = [list(bounds[2:]), list(bounds[0:2])],
        density = False,
        weights = gdf[weights]
    )

    # Convert 2D histogram into natural logarithm...
    logheatmap = np.log2(heatmap)
    logheatmap[np.isneginf(logheatmap)] = 0
    
    # ... and apply Gaussian filter
    logheatmap = ndimage.filters.gaussian_filter(
        logheatmap,
        smoothing,
        mode = "nearest",
        radius = 10
    )

    # Show logarithmic heatmap
    im = ax.imshow(
        # heatmap,
        # np.ma.masked_where(heatmap < 0, heatmap),
        np.ma.masked_where(logheatmap < 0.01, logheatmap),
        cmap = color_scheme,
        alpha = alpha,
        origin = "lower",
        extent = bounds
    )

    # Show plot...
    return im