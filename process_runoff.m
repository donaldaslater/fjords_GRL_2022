% script to process Mankoff 2020 runoff for Slater 2022 GRL
clear; close all;

% Mankoff 2020 runoff is the dataset described here
% https://essd.copernicus.org/articles/12/2811/2020/
% data files are in the format RACMO_xxxx.nc or MAR_xxxx.nc
% where xxxx is the year

% local location of Mankoff 2020 runoff dataset files
p2r = '~/OneDrive - University of Edinburgh/data/mankoff_runoff/';

% read 2010-2019 runoff
yrs = [2010:2019];
t0_RACMO = []; t0_MAR = []; q0_RACMO = []; q0_MAR = [];
% static data - location and elevation where runoff enters the ocean
% note this may differ from where it leaves the ice sheet for
% land-terminating outlets
lat = ncread([p2r,'RACMO_',num2str(yrs(1)),'.nc'],'coast_lat');
lon = ncread([p2r,'RACMO_',num2str(yrs(1)),'.nc'],'coast_lon');
z = ncread([p2r,'RACMO_',num2str(yrs(1)),'.nc'],'coast_alt');
[x,y] = latlon2utm(lat,lon);
% read data
for ii=1:length(yrs),
    % time in matlab serial date number
    t0_RACMO = [t0_RACMO,double(datenum(yrs(ii),1,1) + ...
        ncread([p2r,'RACMO_',num2str(yrs(ii)),'.nc'],'time')')];
    t0_MAR = [t0_MAR,double(datenum(yrs(ii),1,1) + ...
        ncread([p2r,'MAR_',num2str(yrs(ii)),'.nc'],'time')')];
    % runoff in m3/s
    q0_RACMO = [q0_RACMO,ncread([p2r,'RACMO_',num2str(yrs(ii)),'.nc'],'runoff')'];
    q0_MAR = [q0_MAR,ncread([p2r,'MAR_',num2str(yrs(ii)),'.nc'],'runoff')'];
end

% remove those with NaNs
inds = find(isnan(q0_MAR(:,1)));
x(inds) = [];
y(inds) = [];
z(inds) = [];
lat(inds) = [];
lon(inds) = [];
q0_RACMO(inds,:) = [];
q0_MAR(inds,:) = [];

% convert runoff to monthly and time to decimal year
% note MAR has all of 2019
t_MAR = []; q_MAR = [];
dv_MAR = datevec(t0_MAR);
yrs = unique(dv_MAR(:,1));
for ii=1:length(yrs),
    mnths = unique(dv_MAR(find(dv_MAR(:,1)==yrs(ii)),2))';
    t_MAR = [t_MAR,yrs(ii)+(mnths-0.5)/12];
    for jj=1:length(mnths),
        q_MAR(:,end+1) = mean(q0_MAR(:,find(dv_MAR(:,1)==yrs(ii) & dv_MAR(:,2)==mnths(jj))),2);
    end
end
%  whereas RACMO only goes up to August 2019
t_RACMO = []; q_RACMO = [];
dv_RACMO = datevec(t0_RACMO);
yrs = unique(dv_RACMO(:,1));
for ii=1:length(yrs),
    mnths = unique(dv_RACMO(find(dv_RACMO(:,1)==yrs(ii)),2))';
    t_RACMO = [t_RACMO,yrs(ii)+(mnths-0.5)/12];
    for jj=1:length(mnths),
        q_RACMO(:,end+1) = mean(q0_RACMO(:,find(dv_RACMO(:,1)==yrs(ii) & dv_RACMO(:,2)==mnths(jj))),2);
    end
end

% there appear to be some outlets with exactly the same position
% where this occurs combine them into one
i = 1;
while i<=length(x),
    d = (x-x(i)).^2+(y-y(i)).^2;
    if sum(d==0)>1,
        inds = find(d==0);
        q_MAR(inds(1),:) = sum(q_MAR(inds,:));
        q_RACMO(inds(1),:) = sum(q_RACMO(inds,:));
        x(inds(2:end)) = [];
        y(inds(2:end)) = [];
        z(inds(2:end)) = [];
        lat(inds(2:end)) = [];
        lon(inds(2:end)) = [];
        q_MAR(inds(2:end),:) = [];
        q_RACMO(inds(2:end),:) = [];
    end
    i = i+1;
end

% also save mean of RACMO and MAR during overlapping time period
t = t_RACMO;
q = 0.5*(q_RACMO+q_MAR(:,1:length(t)));

% place into structure and save as final runoff structure
runoff.x = x;
runoff.y = y;
runoff.z = z;
runoff.lat = lat;
runoff.lon = lon;
runoff.t = t;
runoff.q = q;
runoff.t_MAR = t_MAR;
runoff.q_MAR = q_MAR;
runoff.t_RACMO = t_RACMO;
runoff.q_RACMO = q_RACMO;
save runoff.mat runoff

% % optional figure to see where outlets are
% plotgreenland;
% hold on;
% plot(runoff.x,runoff.y,'ko');