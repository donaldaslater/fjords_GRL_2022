% script to create a list of glaciers for Slater 2022 GRL
clear; close all;

% the list of glaciers is available as Table S1 of Morlighem 2017
% https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2017GL074954
ts1 = '~/OneDrive - University of Edinburgh/data/master2/Morlighem2017_majorglaciers.xlsx';

% read table
[a,b] = xlsread(ts1);

for ii=1:max(a(:,1)),
    glaciers(ii).lat = a(ii,3);
    glaciers(ii).lon = a(ii,4);
    glaciers(ii).gldepth = a(ii,12);
    glaciers(ii).name = b{ii+2,1};
    glaciers(ii).morlighem_id = ii;
    [glaciers(ii).x,glaciers(ii).y] = latlon2utm(glaciers(ii).lat,glaciers(ii).lon);
end

save glaciers.mat glaciers

