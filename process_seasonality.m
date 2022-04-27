% script to investigate seasonality for Slater 2022 GRL
clear; close all;

% load master dataset
load twglaciers.mat

% load CTD profiles from Carroll 2016
% https://onlinelibrary.wiley.com/doi/abs/10.1002/2016GL070170
load carroll_2016_CTDs.mat

% smoothing window for profiles
sm = 20;

% default plume width and entrainment values
w0 = twglaciers(1).plume(1).w;
E0 = twglaciers(1).plume(1).E;

% rink
gid = 5;

% summer
ff = 'rink_sep_2013';
% profile to use (chosen as close to ~4 km from front to match OMG)
ind = 1;
% cycle over profiles and run plume model  
zi = [twglaciers(gid).gldepth:1:0];
Ti = interp1(-ctd.(ff).depth,ctd.(ff).temp(:,ind),zi,'nearest','extrap');
Si = interp1(-ctd.(ff).depth,ctd.(ff).salt(:,ind),zi,'nearest','extrap'); 
Ti = movmean(Ti,sm); Si = movmean(Si,sm);
sol = run_plume(double(zi),double(0*zi),double(Ti),double(Si),...
    double(0*zi),double(twglaciers(gid).Qsg0/w0),E0);
twglaciers(gid).seasonality.summer.nb = sol.zNB;
twglaciers(gid).seasonality.summer.Ta = sol.Ta;
twglaciers(gid).seasonality.summer.Sa = sol.Sa;
twglaciers(gid).seasonality.summer.z = sol.z;
twglaciers(gid).seasonality.summer.sta = rho(sol.Ta,sol.Sa,0)-1000;
twglaciers(gid).seasonality.summer.st = rho(sol.T,sol.S,0)-1000;
twglaciers(gid).seasonality.summer.name = ff;

% winter
ff = 'rink_july_2013';
% profile to use (chosen as close to ~4 km from front to match OMG)
ind = 1;
% cycle over profiles and run plume model
zi = [twglaciers(gid).gldepth:1:0];
Ti = interp1(-ctd.(ff).depth,ctd.(ff).temp(:,ind),zi,'nearest','extrap');
Si = interp1(-ctd.(ff).depth,ctd.(ff).salt(:,ind),zi,'nearest','extrap');
Ti = movmean(Ti,sm); Si = movmean(Si,sm);
sol = run_plume(double(zi),double(0*zi),double(Ti),double(Si),...
    double(0*zi),double(twglaciers(gid).Qsg0/w0),E0);
twglaciers(gid).seasonality.winter.nb = sol.zNB;
twglaciers(gid).seasonality.winter.Ta = sol.Ta;
twglaciers(gid).seasonality.winter.Sa = sol.Sa;
twglaciers(gid).seasonality.winter.z = sol.z;
twglaciers(gid).seasonality.winter.sta = rho(sol.Ta,sol.Sa,0)-1000;
twglaciers(gid).seasonality.winter.st = rho(sol.T,sol.S,0)-1000;
twglaciers(gid).seasonality.winter.name = ff;

% interannual
ffs = {'rink_july_2013','rink_july_2014','rink_july_2015'};
% profiles to use (chosen as close to ~4 km from front to match OMG)
inds = [1,1,1];
for j=1:length(ffs),
    ff = ffs{j};
    ind = inds(j);
    % cycle over profiles and run plume model  
    zi = [twglaciers(gid).gldepth:1:0];
    Ti = interp1(-ctd.(ff).depth,ctd.(ff).temp(:,ind),zi,'nearest','extrap');
    Si = interp1(-ctd.(ff).depth,ctd.(ff).salt(:,ind),zi,'nearest','extrap');
    Ti = movmean(Ti,sm); Si = movmean(Si,sm);
    sol = run_plume(double(zi),double(0*zi),double(Ti),double(Si),...
        double(0*zi),double(twglaciers(gid).Qsg0/w0),E0);
    twglaciers(gid).interannual(j).nb = sol.zNB;
    twglaciers(gid).interannual(j).Ta = sol.Ta;
    twglaciers(gid).interannual(j).Sa = sol.Sa;
    twglaciers(gid).interannual(j).z = sol.z;
    twglaciers(gid).interannual(j).sta = rho(sol.Ta,sol.Sa,0)-1000;
    twglaciers(gid).interannual(j).st = rho(sol.T,sol.S,0)-1000;
    twglaciers(gid).interannual(j).name = ff;
end

% helheim
gid = 3;

% summer
ff = 'sermilik_aug_2010';
% profile to use (chosen as close to ~25 km from front to match OMG)
ind = 1;
% cycle over profiles and run plume model  
zi = [twglaciers(gid).gldepth:1:0];
Ti = interp1(-ctd.(ff).depth,ctd.(ff).temp(:,ind),zi,'nearest','extrap');
Si = interp1(-ctd.(ff).depth,ctd.(ff).salt(:,ind),zi,'nearest','extrap');
Ti = movmean(Ti,sm); Si = movmean(Si,sm);
sol = run_plume(double(zi),double(0*zi),double(Ti),double(Si),...
    double(0*zi),double(twglaciers(gid).Qsg0/w0),E0);
twglaciers(gid).seasonality.summer.nb = sol.zNB;
twglaciers(gid).seasonality.summer.Ta = sol.Ta;
twglaciers(gid).seasonality.summer.Sa = sol.Sa;
twglaciers(gid).seasonality.summer.z = sol.z;
twglaciers(gid).seasonality.summer.sta = rho(sol.Ta,sol.Sa,0)-1000;
twglaciers(gid).seasonality.summer.st = rho(sol.T,sol.S,0)-1000;
twglaciers(gid).seasonality.summer.name = ff;

% winter
ff = 'sermilik_mar_2010';
% profile to use (chosen as close to ~25 km from front to match OMG)
ind = 2;
% cycle over profiles and run plume model
zi = [twglaciers(gid).gldepth:1:0];
Ti = interp1(-ctd.(ff).depth,ctd.(ff).temp(:,ind),zi,'nearest','extrap');
Si = interp1(-ctd.(ff).depth,ctd.(ff).salt(:,ind),zi,'nearest','extrap');
Ti = movmean(Ti,sm); Si = movmean(Si,sm);
sol = run_plume(double(zi),double(0*zi),double(Ti),double(Si),...
    double(0*zi),double(twglaciers(gid).Qsg0/w0),E0);
twglaciers(gid).seasonality.winter.nb = sol.zNB;
twglaciers(gid).seasonality.winter.Ta = sol.Ta;
twglaciers(gid).seasonality.winter.Sa = sol.Sa;
twglaciers(gid).seasonality.winter.z = sol.z;
twglaciers(gid).seasonality.winter.sta = rho(sol.Ta,sol.Sa,0)-1000;
twglaciers(gid).seasonality.winter.st = rho(sol.T,sol.S,0)-1000;
twglaciers(gid).seasonality.winter.name = ff;

% interannual
ffs = {'sermilik_sep_2008','sermilik_aug_2009','sermilik_aug_2010',...
       'sermilik_aug_2011','sermilik_sep_2012'};
% profiles to use (chosen as close to ~25 km from front to match OMG)
inds = [1,4,1,2,1];
for j=1:length(ffs),
    ff = ffs{j};
    ind = inds(j);
    % cycle over profiles and run plume model  
    zi = [twglaciers(gid).gldepth:1:0];
    Ti = interp1(-ctd.(ff).depth,ctd.(ff).temp(:,ind),zi,'nearest','extrap');
    Si = interp1(-ctd.(ff).depth,ctd.(ff).salt(:,ind),zi,'nearest','extrap');
    Ti = movmean(Ti,sm); Si = movmean(Si,sm);
    sol = run_plume(double(zi),double(0*zi),double(Ti),double(Si),...
        double(0*zi),double(twglaciers(gid).Qsg0/w0),E0);
    twglaciers(gid).interannual(j).nb = sol.zNB;
    twglaciers(gid).interannual(j).Ta = sol.Ta;
    twglaciers(gid).interannual(j).Sa = sol.Sa;
    twglaciers(gid).interannual(j).z = sol.z;
    twglaciers(gid).interannual(j).sta = rho(sol.Ta,sol.Sa,0)-1000;
    twglaciers(gid).interannual(j).st = rho(sol.T,sol.S,0)-1000;
    twglaciers(gid).interannual(j).name = ff;
end

%% save data

save twglaciers.mat twglaciers