% script to link fjords, runoff and glacier data for Slater 2022 GRL
clear; close all;

% load datasets prepared so far
load fjords.mat
load runoff.mat
load glaciers.mat
load calvingfronts.mat

% load bedmachine (for purpose of plotting)
bmfile = '~/Documents/BedMachineGreenland-2021-04-20.nc';
bm.x = single(ncread(bmfile,'x'));
bm.y = single(ncread(bmfile,'y'));
bm.t = single(ncread(bmfile,'mask'))';

% load OMG CTD casts
load all_OMG_CTD_from_csv_ver_1.3.mat

%% link glaciers to calving fronts

% when linking glaciers to calving front by closest distance
% it is necessary to move some glaciers to get the right match
% this is done based on the plot below
glaciers(16).x = -389075; glaciers(16).y = -1473720;
glaciers(144).x = -82175; glaciers(144).y = -910175;
glaciers(216).x = 548125; glaciers(216).y = -1329420;
glaciers(161).x = 659575; glaciers(161).y = -2271580;
glaciers(102).x = 664225; glaciers(102).y = -2266480;
glaciers(121).x = 370375; glaciers(121).y = -2562880;
glaciers(122).x = 115675; glaciers(122).y = -3120280;
glaciers(128).x = 129775; glaciers(128).y = -3152820;
glaciers(108).x = 85525; glaciers(108).y = -3214920;

% link glaciers to calving fronts
for i=1:length(glaciers),
    d = sqrt(([cfs.x]-glaciers(i).x).^2+([cfs.y]-glaciers(i).y).^2);
    [~,ind] = min(d);
    twglaciers(i).x = cfs(ind).x;
    twglaciers(i).y = cfs(ind).y;
    twglaciers(i).gldepth = cfs(ind).gldepth;
end

% add metadata
for i=1:length(glaciers),
    twglaciers(i).morlighem_id = glaciers(i).morlighem_id;
    twglaciers(i).name = glaciers(i).name;
    [twglaciers(i).lat,twglaciers(i).lon] = polarstereo_inv(twglaciers(i).x,twglaciers(i).y);
end

% % optional plot to check the linking
% figure(); hold on;
% imagesc(bm.x,bm.y,bm.t); axis xy equal;
% for i=1:length(twglaciers),
%     plot([twglaciers(i).x,glaciers(i).x],[twglaciers(i).y,glaciers(i).y],'ko-');
%     text(glaciers(i).x,glaciers(i).y,num2str(i),'fontsize',16);
% end

% manually remove those that are land-terminating using the plot
% 81 is a special case - the calving front is the same as 33
inds = [164,205,172,178,197,173,189,182,207,176,165,148,213,206,196,...
        222,231,203,200,187,81];
disp(['Removing ',num2str(length(inds)),' land-terminating glaciers']);
twglaciers(inds) = [];

% remove glaciers with grounding lines shallower than 50 m
inds = find([twglaciers.gldepth]>-50);
disp(['Removing ',num2str(length(inds)),' glaciers with grounding line shallower than 50 m']);
twglaciers(inds) = [];

%% link glaciers to runoff

% we will not consider glaciers with mean annual runoff < 2.5 m3/s
inds = find(mean(runoff.q_MAR,2)<2.5);
runoff.x(inds) = [];
runoff.y(inds) = [];
runoff.z(inds) = [];
runoff.lat(inds) = [];
runoff.lon(inds) = [];
runoff.q(inds,:) = [];
runoff.q_MAR(inds,:) = [];
runoff.q_RACMO(inds,:) = [];

% for each glacier, find the closest large runoff outlet
for i=1:length(twglaciers),
    d = sqrt(([runoff.x]-twglaciers(i).x).^2+([runoff.y]-twglaciers(i).y).^2);
    [~,id] = min(d);
    twglaciers(i).runoff.t = runoff.t;
    twglaciers(i).runoff.q = runoff.q(id,:);
    twglaciers(i).runoff.t_MAR = runoff.t_MAR;
    twglaciers(i).runoff.q_MAR = runoff.q_MAR(id,:);
    twglaciers(i).runoff.t_RACMO = runoff.t_RACMO;
    twglaciers(i).runoff.q_RACMO = runoff.q_RACMO(id,:);
    twglaciers(i).runoff.x = runoff.x(id);
    twglaciers(i).runoff.y = runoff.y(id);
end

% as before, unfortunately some matches require manual intervention
% perhaps due to differences in mask between glaciers definition and
% runoff routing
% alison glacier
gi(1) = find([twglaciers.morlighem_id]==17);
[~,id(1)] = min((runoff.x-(-332650)).^2+(runoff.y-(-1642550)).^2);
% nordenskiold glacier
gi(2) = find([twglaciers.morlighem_id]==126);
[~,id(2)] = min((runoff.x-(-378150)).^2+(runoff.y-(-1497750)).^2);
% heilprin glacier
gi(3) = find([twglaciers.morlighem_id]==69);
[~,id(3)] = min((runoff.x-(-488550)).^2+(runoff.y-(-1271850)).^2);
% ryder glacier
gi(4) = find([twglaciers.morlighem_id]==144);
[~,id(4)] = min((runoff.x-(-85050)).^2+(runoff.y-(-879250)).^2);
% kangerlussuaq glacier
gi(5) = find([twglaciers.morlighem_id]==2);
[~,id(5)] = min((runoff.x-(497750)).^2+(runoff.y-(-2293350)).^2);
% polaric glacier
gi(6) = find([twglaciers.morlighem_id]==61);
[~,id(6)] = min((runoff.x-(525450)).^2+(runoff.y-(-2368950)).^2);
% unnamed danell fjord
gi(7) = find([twglaciers.morlighem_id]==108);
[~,id(7)] = min((runoff.x-(85950)).^2+(runoff.y-(-3214550)).^2);
% frederiksborg glacier
gi(8) = find([twglaciers.morlighem_id]==38);
[~,id(8)] = min((runoff.x-(555150)).^2+(runoff.y-(-2316550)).^2);
% 79N
gi(9) = find([twglaciers.morlighem_id]==74);
[~,id(9)] = min((runoff.x-(464450)).^2+(runoff.y-(-1000150)).^2);
% Peterman
gi(10) = find([twglaciers.morlighem_id]==78);
[~,id(10)] = min((runoff.x-(-279750)).^2+(runoff.y-(-924650)).^2);
% hagen brae
gi(11) = find([twglaciers.morlighem_id]==111);
[~,id(11)] = min((runoff.x-(286850)).^2+(runoff.y-(-874950)).^2);
% waltershausen
gi(12) = find([twglaciers.morlighem_id]==215);
[~,id(12)] = min((runoff.x-(629450)).^2+(runoff.y-(-1646050)).^2);
% daugard-jensen
gi(13) = find([twglaciers.morlighem_id]==8);
[~,id(13)] = min((runoff.x-(559550)).^2+(runoff.y-(-1894850)).^2);
% jaette
gi(14) = find([twglaciers.morlighem_id]==185);
[~,id(14)] = min((runoff.x-(545550)).^2+(runoff.y-(-1720950)).^2);
% update runoff outlet for these systems
for jj=1:length(gi),
    twglaciers(gi(jj)).runoff.t = runoff.t;
    twglaciers(gi(jj)).runoff.q = runoff.q(id(jj),:);
    twglaciers(gi(jj)).runoff.t_MAR = runoff.t_MAR;
    twglaciers(gi(jj)).runoff.q_MAR = runoff.q_MAR(id(jj),:);
    twglaciers(gi(jj)).runoff.t_RACMO = runoff.t_RACMO;
    twglaciers(gi(jj)).runoff.q_RACMO = runoff.q_RACMO(id(jj),:);
    twglaciers(gi(jj)).runoff.x = runoff.x(id(jj));
    twglaciers(gi(jj)).runoff.y = runoff.y(id(jj));
end
% glaciers with no significant runoff outlet according to routing
inds_to_remove = [144,148,129,130,140,114,151];
disp(['Removing ',num2str(length(inds_to_remove)),' glaciers without a significant runoff outlet']);
twglaciers(inds_to_remove) = [];

% optional plot to check the linking
% figure(); hold on;
% imagesc(bm.x,bm.y,bm.t); axis xy equal;
% scatter([runoff.x],[runoff.y],mean(runoff.q_MAR,2),'r','filled');
% for i=1:length(twglaciers),
%     plot([twglaciers(i).x,twglaciers(i).runoff.x],[twglaciers(i).y,twglaciers(i).runoff.y],'ko-');
%     text(twglaciers(i).x,twglaciers(i).y,num2str(twglaciers(i).morlighem_id),'fontsize',16);
% end

% calculate headline runoff stats
tlims = [2010,2020];
for i=1:length(twglaciers),
    tq = twglaciers(i).runoff.t;
    tq_inds = find(tq>tlims(1) & tq<tlims(2) & tq-floor(tq)>0.4 & tq-floor(tq)<0.7);
    % mean JJA runoff during tlims
    twglaciers(i).Qsg0 = mean(twglaciers(i).runoff.q(tq_inds));
    % min summer monthly runoff during tlims
    twglaciers(i).Qsgmin = min(twglaciers(i).runoff.q(tq_inds));
    % max summer monthly runoff during tlims
    twglaciers(i).Qsgmax = max(twglaciers(i).runoff.q(tq_inds));
end

%% link glaciers to fjords

% get individual fjords
fi = bwconncomp(fb,4);
% number of unique fjords
disp([num2str(fi.NumObjects),' unique fjords']);

% for each glacier, find the closest fjord
[X,Y] = meshgrid(x,y);
for i=1:length(twglaciers),
    for j=1:length(fi.PixelIdxList),
        xfjord = X(fi.PixelIdxList{j});
        yfjord = Y(fi.PixelIdxList{j});
        dfjord = sqrt((xfjord-twglaciers(i).x).^2+(yfjord-twglaciers(i).y).^2);
        d(j) = min(dfjord);
    end
    [~,id] = min(d);
    twglaciers(i).fjord.id = id;
    twglaciers(i).fjord.inds = fi.PixelIdxList{id};
    twglaciers(i).fjord.x = mean(X(twglaciers(i).fjord.inds));
    twglaciers(i).fjord.y = mean(Y(twglaciers(i).fjord.inds));
end

% % optional plot for manual check
% figure(); hold on;
% imagesc(x,y,f); axis xy equal;
% for i=1:length(twglaciers),
%     plot([twglaciers(i).x,twglaciers(i).fjord.x],...
%          [twglaciers(i).y,twglaciers(i).fjord.y],'ro-');
%     plot(twglaciers(i).x,twglaciers(i).y,'ko');
%     text(twglaciers(i).x,twglaciers(i).y,num2str(twglaciers(i).morlighem_id),'fontsize',16);
% end

% calculate basic fjord properties

% resolution
dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
da = dx*dy;

% area, volume, mean and max depth
for i=1:length(twglaciers),
    twglaciers(i).fjord.area = length(twglaciers(i).fjord.inds)*da/1e6; % sq km
    twglaciers(i).fjord.vol = -sum(b(twglaciers(i).fjord.inds))*da/1e9; % cubic km
    twglaciers(i).fjord.meandepth = 1000*twglaciers(i).fjord.vol/twglaciers(i).fjord.area; % m
    twglaciers(i).fjord.maxdepth = -min(b(twglaciers(i).fjord.inds));
end

% fjord length
for i=1:length(twglaciers);
    % find pixels that are adjacent to the shelf
    twglaciers(i).fjord.b_inds = [];
    inds = twglaciers(i).fjord.inds;
    [inds_i,inds_j]=ind2sub(size(b),inds);
    for j=1:length(inds),
        kern = f(inds_i(j)-1:inds_i(j)+1,inds_j(j)-1:inds_j(j)+1);
        kern(1,1)=NaN; kern(3,3)=NaN; kern(1,3)=NaN; kern(3,1)=NaN;
        if any(kern(:)==0),
            twglaciers(i).fjord.b_inds = [twglaciers(i).fjord.b_inds,inds(j)];
        end
    end
    % find glacier pixel
    dfjord = sqrt((twglaciers(i).x-X(inds)).^2+(twglaciers(i).y-Y(inds)).^2);
    [~,id] = min(dfjord);
    twglaciers(i).fjord.g_ind = inds(id(1));
    % make mask that is only 1 for this fjord
    fjordmask = 0*f;
    fjordmask(inds) = 1;
    % get distances through fjord from each point to glacier
    D = (dx/1000)*bwdistgeodesic(logical(fjordmask),twglaciers(i).fjord.g_ind,'quasi-euclidean');
    % length is the minimum of D on the points adjacent to the shelf
    twglaciers(i).fjord.length = min(D(twglaciers(i).fjord.b_inds));
end

% fjord effective depth
for i=1:length(twglaciers),
    twglaciers(i).fjord.effdepth = ed(twglaciers(i).fjord.g_ind);
end

% remove glaciers with effective depth shallower than 25 m but before
% doing that, intervene for Sarqardleq - bathymetry data from
% Stevens et al. (2016) shows an effective depth of 70 m but this data is
% not in BedMachine
i_sar = find([twglaciers.morlighem_id]==168);
twglaciers(i_sar).fjord.effdepth = -70;
% now remove shallow effective depths
inds_effdepth = [];
for i=1:length(twglaciers),
    if twglaciers(i).fjord.effdepth>-25,
        inds_effdepth = [inds_effdepth,i];
    end
end
twglaciers(inds_effdepth) = [];
disp(['Removing ',num2str(length(inds_effdepth)),' glaciers with effective depth shallower than 25 m']);

% if fjord assigned to glacier has no shelf-adjacent points
% with ss=5 in define_fjords.m script this shouldn't result in anything
for i=1:length(twglaciers),
    if isempty(twglaciers(i).fjord.b_inds),
        twglaciers(i).fjord.id = NaN;
        twglaciers(i).fjord.inds = NaN;
        twglaciers(i).fjord.area = NaN;
        twglaciers(i).fjord.vol = NaN;
        twglaciers(i).fjord.meandepth = NaN;
        twglaciers(i).fjord.maxdepth = NaN;
        twglaciers(i).fjord.b_inds = NaN;
        twglaciers(i).fjord.g_ind = NaN;
        twglaciers(i).fjord.length = NaN;
    end
end
        
% fjord width
for i=1:length(twglaciers),
    twglaciers(i).fjord.width = twglaciers(i).fjord.area/twglaciers(i).fjord.length;
end

% how many fjords for all the glaciers?
for i=1:length(twglaciers),
    fjords(i) = twglaciers(i).fjord.id;
end
disp([num2str(length(twglaciers)),' twglaciers drain into ',...
    num2str(length(unique(fjords))),' fjords']);

%% link glaciers to CTDs

% process cast locations
for i=1:length(cast),
    omglat(i) = cast{i}.lat(1);
    omglon(i) = cast{i}.lon(1);
end
[omgx,omgy] = latlon2utm(omglat,omglon);

% a few casts need bad data removed manually
[~,id]=min((omgx-(-323552)).^2+(omgy-(-1801987)).^2);
inds=[23:28]; cast{id}.depth(inds)=[]; cast{id}.temp(inds)=[]; cast{id}.sal(inds)=[];
[~,id]=min((omgx-(-486218)).^2+(omgy-(-1250280)).^2);
inds=[13]; cast{id}.depth(inds)=[]; cast{id}.temp(inds)=[]; cast{id}.sal(inds)=[];
[~,id]=min((omgx-(-374524)).^2+(omgy-(-1509113)).^2);
inds=[6:10,18]; cast{id}.depth(inds)=[]; cast{id}.temp(inds)=[]; cast{id}.sal(inds)=[];
[~,id]=min((omgx-(-209721)).^2+(omgy-(-2097227)).^2);
inds=[51:53]; cast{id}.depth(inds)=[]; cast{id}.temp(inds)=[]; cast{id}.sal(inds)=[];
[~,id]=min((omgx-(-457226)).^2+(omgy-(-1414614)).^2);
inds=[15,19:22]; cast{id}.depth(inds)=[]; cast{id}.temp(inds)=[]; cast{id}.sal(inds)=[];
[~,id]=min((omgx-(-306533)).^2+(omgy-(-1827521)).^2);
inds=[47:50]; cast{id}.depth(inds)=[]; cast{id}.temp(inds)=[]; cast{id}.sal(inds)=[];
[~,id]=min((omgx-(-433061)).^2+(omgy-(-1432208)).^2);
inds=[34:42]; cast{id}.depth(inds)=[]; cast{id}.temp(inds)=[]; cast{id}.sal(inds)=[];
[~,id]=min((omgx-(-523376)).^2+(omgy-(-1393094)).^2);
inds=[25:28]; cast{id}.depth(inds)=[]; cast{id}.temp(inds)=[]; cast{id}.sal(inds)=[];

% associate casts to glaciers based on proximity
for i=1:length(twglaciers),
    d = (omgx-twglaciers(i).x).^2+(omgy-twglaciers(i).y).^2;
    [~,id] = min(d);
    twglaciers(i).profile.x = omgx(id);
    twglaciers(i).profile.y = omgy(id);
    twglaciers(i).profile.t = cast{id}.time;
    twglaciers(i).profile.z = [twglaciers(i).gldepth:1:0];
    % interpolate profiles onto z
    twglaciers(i).profile.T = interp1(-cast{id}.depth,cast{id}.temp,twglaciers(i).profile.z,'linear',NaN);
    twglaciers(i).profile.S = interp1(-cast{id}.depth,cast{id}.sal,twglaciers(i).profile.z,'linear',NaN);
    % replace possible NaNs at top or bottom with nearest value
    if isnan(twglaciers(i).profile.T(1)),
        j=1;
        while isnan(twglaciers(i).profile.T(j+1)), j=j+1; end
        twglaciers(i).profile.T(1:j) = twglaciers(i).profile.T(j+1);
        twglaciers(i).profile.S(1:j) = twglaciers(i).profile.S(j+1);
    end
    if isnan(twglaciers(i).profile.T(end)),
        j=length(twglaciers(i).profile.T);
        while isnan(twglaciers(i).profile.T(j-1)), j=j-1; end
        twglaciers(i).profile.T(j:end) = twglaciers(i).profile.T(j-1);
        twglaciers(i).profile.S(j:end) = twglaciers(i).profile.S(j-1);
    end
end

% as usual, a number need manual intervention to link the best cast
clearvars gi id
gi(1) = find([twglaciers.morlighem_id]==66);
[~,id(1)] = min((omgx-(-64364)).^2+(omgy-(-3210540)).^2);
gi(2) = find([twglaciers.morlighem_id]==134);
[~,id(2)] = min((omgx-(-221837)).^2+(omgy-(-2131880)).^2);
gi(3) = find([twglaciers.morlighem_id]==74);
[~,id(3)] = min((omgx-(487259)).^2+(omgy-(-1024170)).^2);
gi(4) = find([twglaciers.morlighem_id]==119);
[~,id(4)] = min((omgx-(-256318)).^2+(omgy-(-2102420)).^2);
% gi(5) = find([twglaciers.morlighem_id]==191);
% [~,id(5)] = min((omgx-(-277878)).^2+(omgy-(-2052460)).^2);
gi(5) = find([twglaciers.morlighem_id]==16);
[~,id(5)] = min((omgx-(-392847)).^2+(omgy-(-1478040)).^2);
gi(6) = find([twglaciers.morlighem_id]==114);
[~,id(6)] = min((omgx-(-491461)).^2+(omgy-(-1409740)).^2);
gi(7) = find([twglaciers.morlighem_id]==146);
[~,id(7)] = min((omgx-(-544276)).^2+(omgy-(-1409590)).^2);
gi(8) = find([twglaciers.morlighem_id]==219);
[~,id(8)] = min((omgx-(-584212)).^2+(omgy-(-1333190)).^2);
% gi(10) = find([twglaciers.morlighem_id]==240);
% [~,id(10)] = min((omgx-(612002)).^2+(omgy-(-1146980)).^2);
% gi(11) = find([twglaciers.morlighem_id]==154);
% [~,id(11)] = min((omgx-(644144)).^2+(omgy-(-2316710)).^2);
gi(9) = find([twglaciers.morlighem_id]==14);
[~,id(9)] = min((omgx-(194847)).^2+(omgy-(-2747640)).^2);
% gi(13) = find([twglaciers.morlighem_id]==20);
% [~,id(13)] = min((omgx-(-64365)).^2+(omgy-(-3210540)).^2);
gi(10) = find([twglaciers.morlighem_id]==3);
[~,id(10)] = min((omgx-(320705)).^2+(omgy-(-2579780)).^2);

for jj=1:length(gi),
    twglaciers(gi(jj)).profile.x = omgx(id(jj));
    twglaciers(gi(jj)).profile.y = omgy(id(jj));
    twglaciers(gi(jj)).profile.t = cast{id(jj)}.time;
    twglaciers(gi(jj)).profile.T = interp1(-cast{id(jj)}.depth,cast{id(jj)}.temp,twglaciers(gi(jj)).profile.z,'nearest','extrap');
    twglaciers(gi(jj)).profile.S = interp1(-cast{id(jj)}.depth,cast{id(jj)}.sal,twglaciers(gi(jj)).profile.z,'nearest','extrap');
end

% a number require datasets other than OMG

% kns, narssap sermia and qooqqup sermia
% data from mortensen (2020) available at https://doi.pangaea.de/10.1594/PANGAEA.921991
load ~/'OneDrive - University of Edinburgh'/biogeochem/Mortensen-Meire_2020/datasets/mortensen_meire.mat
i_kns = find([twglaciers.morlighem_id]==36);
[twglaciers(i_kns).profile.x,twglaciers(i_kns).profile.y] = latlon2utm(kns.lat,kns.lon);
twglaciers(i_kns).profile.t = kns.t;
twglaciers(i_kns).profile.T = interp1(kns.z,kns.T,twglaciers(i_kns).profile.z,'nearest','extrap');
twglaciers(i_kns).profile.S = interp1(kns.z,kns.S,twglaciers(i_kns).profile.z,'nearest','extrap');
i_nars = find([twglaciers.morlighem_id]==63);
[twglaciers(i_nars).profile.x,twglaciers(i_nars).profile.y] = latlon2utm(nars.lat,nars.lon);
twglaciers(i_nars).profile.t = nars.t;
twglaciers(i_nars).profile.T = interp1(nars.z,nars.T,twglaciers(i_nars).profile.z,'nearest','extrap');
twglaciers(i_nars).profile.S = interp1(nars.z,nars.S,twglaciers(i_nars).profile.z,'nearest','extrap');
i_qooq = find([twglaciers.morlighem_id]==66);
[twglaciers(i_qooq).profile.x,twglaciers(i_qooq).profile.y] = latlon2utm(qooq.lat,qooq.lon);
twglaciers(i_qooq).profile.t = qooq.t;
twglaciers(i_qooq).profile.T = interp1(qooq.z,qooq.T,twglaciers(i_qooq).profile.z,'nearest','extrap');
twglaciers(i_qooq).profile.S = interp1(qooq.z,qooq.S,twglaciers(i_qooq).profile.z,'nearest','extrap');

% sarqardleq
% data from Fiamma Straneo. 2019. Temperature and salinity profiles...
% adjacent to a tidewater glacier in Sarqardleq Fjord, West Greenland, ...
% collected during July 2013. Arctic Data Center. doi:10.18739/A2B853H78
load ~/'OneDrive - University of Edinburgh'/biogeochem/sarqardleq/sarqardleq.mat
i_sar = find([twglaciers.morlighem_id]==168);
twglaciers(i_sar).profile.x = sar.x;
twglaciers(i_sar).profile.y = sar.y;
twglaciers(i_sar).profile.t = sar.t;
twglaciers(i_sar).profile.T = interp1(sar.z,sar.T,twglaciers(i_sar).profile.z,'nearest','extrap');
twglaciers(i_sar).profile.S = interp1(sar.z,sar.S,twglaciers(i_sar).profile.z,'nearest','extrap');

% optional plot to check the linking
% figure(); hold on;
% imagesc(bm.x,bm.y,bm.t); axis xy equal;
% for i=1:length(twglaciers),
%     plot([twglaciers(i).x,twglaciers(i).profile.x],[twglaciers(i).y,twglaciers(i).profile.y],'ko-');
%     text(twglaciers(i).x,twglaciers(i).y,num2str(twglaciers(i).morlighem_id),'fontsize',16);
% end
% plot(omgx,omgy,'ro');

% optional plots to check casts are sensible
% figure();
% for i=1:length(twglaciers),
%     subplot(1,2,1); cla; hold on;
%     plot(twglaciers(i).profile.S,twglaciers(i).profile.z);
%     plot(smooth(twglaciers(i).profile.S,10),twglaciers(i).profile.z);
%     title(num2str(i));
%     subplot(1,2,2); cla; hold on;
%     plot(twglaciers(i).profile.T,twglaciers(i).profile.z);
%     plot([0,1],twglaciers(i).gldepth*[1,1],'k--');
%     title(num2str(i));
%     waitforbuttonpress;
% end

% smooth profiles with 10 m moving average
wd = 10;
for i=1:length(twglaciers),
    twglaciers(i).profile.T = smooth(twglaciers(i).profile.T,wd);
    twglaciers(i).profile.S = smooth(twglaciers(i).profile.S,wd);
end

%% save structure

save twglaciers.mat twglaciers