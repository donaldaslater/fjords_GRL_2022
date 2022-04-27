# fjords_GRL_2022

This directory contains code and data required to reproduce the results in the paper

Characteristic depths, fluxes and timescales for Greenland's tidewater glacier fjords from subglacial discharge-driven upwelling during summer
D. A. Slater, D. Carroll, H. Oliver, M. J. Hopwood, F. Straneo, M. Wood, J. K. Willis, M. Morlighem
Geophysical Research Letters, 2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
It contains the following code files (all MATLAB):

process_calvingfronts.m - defines calving fronts
process_fjords.m - defines fjords
process_runoff.m - gets subglacial runoff from Mankoff 2020
process_glaciers.m - defines glaciers
link_data.m - brings datasets together (calving fronts, fjords, glaciers, CTDs)
process_plume.m - runs plume model for all glaciers
process_seasonality.m - investigates plume seasonality

rho.m - plume model equation of state
run_plume.m - plume model top level code

makeplots.m - makes plots used in the manuscript
maketable.m - makes the csv data table in the SI
fjord_plots.m - makes fjord-by-fjord map plots in the SI

polarstereo_inv.m - converts polar stereographic to lat/lon
latlon2utm.m - converts lat/lon to polar stereographic
custom_legend.m - simple code for custom plot legend
errplot.m - a simple custom plotting function called in makeplots.m
plotbox.m - another simple custom plotting function called in makeplots.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
It contains the following data files (all .mat):

1. calvingfronts.mat - these are the calving fronts created by process_calvingfronts.m
2. fjords.mat - this is the fjord mask created by process_fjords.m
3. runoff.mat - subglacial runoff from Mankoff 2020
4. glaciers.mat - this is the glaciers structure created by process_glaciers.m

5. twglaciers.mat - this is the main output, which is a matlab structure containing the 136 tidewater glaciers and their associated plume properties etc.
6. ice_ocean_sectors.mat - sector polygons derived in Slater 2019 (doi:10.5194/tc-13-2489-2019)

7. all_OMG_CTD_from_csv_ver_1.3.mat - OMG CTD data from https://omg.jpl.nasa.gov/portal/browse/OMGEV-AXCTD/
8. carroll_2016_CTDs.mat - CTD data from Carroll 2016 (doi:10.1002/2016GL070170)
9. mortensen_meire.mat - CTD data from Mortensen & Meire 2020 (doi:10.1594/PANGAEA.921991)
10. sarqardleq.mat - CTD data from Mankoff 2016


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The following files are needed but are freely available and therefore not included in this directory:

1. BedMachine Greenland v4
2. Plotting colormaps from https://www.mathworks.com/matlabcentral/fileexchange/28943-color-palette-tables-cpt-for-matlab (just colormaps, not 100% necessary)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Running order:

To reproduce the results, run the scripts in the order
process_calvingfronts.m
process_fjords.m
process_runoff.m
process_glaciers.m
link_data.m
process_plume.m
process_seasonality.m

To make the plots in the paper and SI, run makeplots.m and/or fjord_plots.m. To make the data table run maketable.m
