% script to run plume model for Slater 2022 GRL
clear; close all;

% load master dataset
load twglaciers.mat

% parameter sweep for plume model
w0 = [100,250,500]; % plume widths
E0 = [0.05,0.1,0.15]; % entrainment coefficients

%% loop over glaciers and run plume model

for i=1:length(twglaciers),

    i

    % parameters result in 7 sets of simulations
    % 1 - central width, mean runoff, central entrainment
    % 2&3 - low and high width
    % 4&5 - low and high runoff
    % 6&7 - low and high entrainment
    w = [w0(2),w0(1),w0(3),w0(2),w0(2),w0(2),w0(2)];
    E = [E0(2),E0(2),E0(2),E0(2),E0(2),E0(1),E0(3)];
    Q = [twglaciers(i).Qsg0,twglaciers(i).Qsg0,...
         twglaciers(i).Qsg0,twglaciers(i).Qsgmin,...
         twglaciers(i).Qsgmax,twglaciers(i).Qsg0,...
         twglaciers(i).Qsg0];

    % run plume model in 7 cases
    for j=1:length(w),
        
        sol = run_plume(double(twglaciers(i).profile.z),0*double(twglaciers(i).profile.z),...
            double(twglaciers(i).profile.T),double(twglaciers(i).profile.S),...
            double(0*twglaciers(i).profile.z),double(Q(j)/w(j)),E(j));

        % grounding line water properties
        [~,gl_ind] = min(twglaciers(i).profile.z);           
        twglaciers(i).plume(j).Tgl = twglaciers(i).profile.T(gl_ind);
        twglaciers(i).plume(j).Sgl = twglaciers(i).profile.S(gl_ind);
        twglaciers(i).plume(j).Qgl = Q(j);

        % neutral buoyancy water properties
        twglaciers(i).plume(j).Tnb = sol.TNB;
        twglaciers(i).plume(j).Snb = sol.SNB;
        twglaciers(i).plume(j).znb = sol.zNB;
        twglaciers(i).plume(j).Qnb = sol.QNB*w(j);

        % save parameters
        twglaciers(i).plume(j).E = E(j);
        twglaciers(i).plume(j).w = w(j);

        % save full plume solution
        twglaciers(i).plume(j).sol = sol;
        
    end

end

%% get volume below neutral buoyancy and renewal time per plume

load fjords.mat
dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
da = dx*dy;

for i=1:length(twglaciers),
    fid(i) = twglaciers(i).fjord.id;
end

% cycle over fjords
for i=1:length(fid),
    % get glaciers discharging into this fjord
    inds = find(fid==fid(i));
    % for these glaciers, get the shallower neutral buoyancy and sum flux
    nb = [];
    flux = 0;
    maxflux = 0;
    minflux = 0;
    for j=1:length(inds),
        nb = [nb,twglaciers(inds(j)).plume(1).znb];
        flux = flux+86400*twglaciers(inds(j)).plume(1).Qnb;
        maxflux = maxflux+86400*max([twglaciers(inds(j)).plume.Qnb]);
        minflux = minflux+86400*min([twglaciers(inds(j)).plume.Qnb]);
    end
    nb = max(nb);
    % get volume below plume NB
    bed = b(twglaciers(inds(1)).fjord.inds);
    bed0 = bed; bed0(bed<nb) = NaN;
    nbvol = (nansum(-bed)-nansum(-bed0))*da/1e9; % cubic km
    for j=1:length(inds),
        twglaciers(inds(j)).fjord.vol_below_NB = nbvol;
        twglaciers(inds(j)).fjord.renewaltime = nbvol*1e9/flux;
        twglaciers(inds(j)).fjord.renewaltimemin = nbvol*1e9/maxflux;
        twglaciers(inds(j)).fjord.renewaltimemax = nbvol*1e9/minflux;
    end
end


%% save data

save twglaciers.mat twglaciers