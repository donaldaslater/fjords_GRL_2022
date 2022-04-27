% script to write out Slater 2022 fjord stats results to csv
clear; close all;

% load dataset
load twglaciers.mat

% put into useful format
for i=1:length(twglaciers),
    names{i} = twglaciers(i).name;
    x(i) = twglaciers(i).x;
    y(i) = twglaciers(i).y;
    lon(i) = twglaciers(i).lon;
    lat(i) = twglaciers(i).lat;
    fjordarea(i) = round(twglaciers(i).fjord.area);
    fjordvol(i) = round(twglaciers(i).fjord.vol);
    fjordmeandepth(i) = twglaciers(i).fjord.meandepth;
    fjordlength(i) = round(twglaciers(i).fjord.length);
    fjordwidth(i) = twglaciers(i).fjord.width;
    gldepth(i) = round(twglaciers(i).gldepth);
    effdepth(i) = round(twglaciers(i).fjord.effdepth);
    nbdepth(i) = round(twglaciers(i).plume(1).znb);
    mhdepth(i) = round(twglaciers(i).plume(1).sol.z(end));
    Q0(i) = round(twglaciers(i).Qsg0);
    Qnb(i) = round(twglaciers(i).plume(1).Qnb);
end

% matlab does not like duplicate glacier names
% so put an underscore in some
inds = find(strcmp(names,'Sermeq Avannarleq'));
names{inds(1)} = 'Sermeq_Avannarleq';
inds = find(strcmp(names,'Nordenskiöld Gletscher'));
names{inds(1)} = 'Nordenskiöld_Gletscher';

T = table(x',y',lon',lat',fjordarea',fjordvol',fjordlength',gldepth',...
    effdepth',nbdepth',mhdepth',Q0',Qnb','RowNames',names);
T.Properties.DimensionNames{1} = 'Glacier Name';
T.Properties.VariableNames{1} = 'UTM x (m)';
T.Properties.VariableNames{2} = 'UTM y (m)';
T.Properties.VariableNames{3} = 'Longitude';
T.Properties.VariableNames{4} = 'Latitude';
T.Properties.VariableNames{5} = 'Fjord area (km^2)';
T.Properties.VariableNames{6} = 'Fjord volume (km^3)';
T.Properties.VariableNames{7} = 'Fjord length (km)';
T.Properties.VariableNames{8} = 'Grounding line depth (m)';
T.Properties.VariableNames{9} = 'Fjord effective depth (m)';
T.Properties.VariableNames{10} = 'Plume neutral buoyancy depth (m)';
T.Properties.VariableNames{11} = 'Plume maximum height (m)';
T.Properties.VariableNames{12} = 'Mean 2010-2020 JJA subglacial discharge (m^3/s)';
T.Properties.VariableNames{13} = 'Plume upwelling flux (m^3/s)';
writetable(T,'all_fjords.csv','writerownames',true);