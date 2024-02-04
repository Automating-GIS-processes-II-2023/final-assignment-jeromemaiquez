[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/6KLiFij3)
# Final Assignment

### Status

Once you are finished with the final assignment, edit this readme and add "x" to the correct box:

* [x] Submitted

* [ ] I'm still working on my final assignment. 


*If you have done this assignement in pair, write your groupmate's name here:* ----

## Topic: Urban Analysis Indicators

The aim of this project was to develop an analytics tool focusing
on urban form. The tool currently only supports two indicators: POI density (currently only skyscrapers) and betweenness centrality of street intersections. The tool is tested for the area surrounding the three biggest central business districts (CBDs) in Metro Manila: Makati, Bonifacio Global City (BGC), and Ortigas Center.

Much of the functionality in this tool is defined in functions in a separate script file. Future improvements (which I can't promise will come) include the ability to select from a pre-defined, bigger set of POIs to map and analyze, cleaned-up code (by wrapping them in script functions), and better documentation overall.

This project was developed as the final assignment for the [AutoGIS](https://autogis-site.readthedocs.io/en/latest/final-assignment/final-assignment.html#urban-indicators) course provided for free by the University of Helsinki. Much thanks to them for emphasizing the importance of open science.

Major functionalities:
- fetch & create maps of street networks and POIs within a buffer of specified radius around a geocoded place name (`get_buffer_from_place()`)
- create heatmaps of relevant point geometries (e.g. skyscrapre POIs, street intersections) within this same buffer (`plot_heatmap()`)
- create plots with a very specific style (circular map, text in the bottom left, line width proportional to lane number, etc.) (`parse_street_lanes()` and `plot_streets()`)

Supported urban indicators:
- street network (geometry)
- street intersection (heatmap)
- skyscraper POI density (heatmap)
- intersection betweenness centrality (heatmap)

Potential urban indicators to add:
- [education facility](https://data.humdata.org/dataset/hotosm_phl_education_facilities) POI density (heatmap)
- [health facility](https://data.humdata.org/dataset/hotosm_phl_health_facilities) POI density (heatmap)
- [financial services](https://data.humdata.org/dataset/hotosm_phl_financial_services) POI density (heatmap)
- commercial amenities POI density (heatmap)
- public institutions POI density (heatmap)
- isochrone from central point in map (geometry?)

### Structure of this repository:
- `osm_street_circle.ipynb`: the final Jupyter notebook
- `final_exercise.ipynb`: a previous attempt... (feel free to ignore)
- `urban_indicators_scripts.py`: script file of all custom functions

### Input data:
- OpenStreetMap street network (via OSMnx)
- OpenStreetMap POIs (via OSMnx)

### Analysis steps:
Part A: Street Network
- create buffer of specified radius from geocoded place name
- extract street network within this buffer

Part B: POI Heatmap
- extract all POIs of specified tags and value within this buffer
- create a heatmap (log 2D histogram) of POIs

Part C: Betweenness Centrality
- calculate betweenness centrality for all street intersections
- assign BC to all nodes and create heatmap

Part D: Plotting it All
- plot street network and heatmap with a specific style
- set title and subtitle of plot

### Results:
- shows spatial distribution and density of skyscraper POIs in Metro Manila, organically highlighting its 3 biggest CBDs (Makati, Bonifacio Global City, Ortigas)
- shows most important intersections in street network in the same area via calculating betweenness centrality for a sample of all nodes in the network
- creates plot of street network with specific style (variable line widths, circular mask, text in bottom-left, etc.)

### References:
- [Top Central Business Districts in Metro Manila](https://santosknightfrank.com/blogs/top-central-business-districts-in-metro-manila/) -- Quick overview of Metro Manila's many CBDs, basic facts and figures, and their histories
- [Notes on graph theory --- Centrality Measures](https://towardsdatascience.com/notes-on-graph-theory-centrality-measurements-e37d2e49550a) -- Overview of common centrality measures for networks/graphs, including betweennes centrality

### Feedback
- `Matplotlib` get very difficult to work with, very quickly -- customizing text and axes positions is very finnicky, and image coordinate systems are confusing
- Why can't `pyplot`'s `hist2d` function take in an `ax` parameter like most other functionalities in the package? Feels like an oversight...
- The difficulty in working with matplotlib emphasizes for me the continuing importance of open-source desktop GIS such as QGIS for geovisualization
- Calculating betweenness centrality via `networkx` gets very slow very soon, so it is very important to set a low `k` parameter to ensure respectable runtimes