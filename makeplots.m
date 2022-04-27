% script to make manuscript plots
function [] = makeplots(plotnum)

% plotnum = 2 for fig. 2
% plotnum = 2.1 for fig. 2 ordered clockwise (fig. S4)
% plotnum = 3 for fig. 3 (old version)
% plotnum = 3.1 for new version of fig. 3 with upwelling sensitivity panel
% plotnum = 4 for fig. 4
% plotnum = 5 for fig. S1
% plotnum = 6 for fig. S2
% plotnum = 7 for fig. S3
% plotnum = 8 for fig. S5
% plotnum = 9 for fig. S6

close all; clearvars -except plotnum
% plotting params
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fs = 8;
bw = 50;
cols = [0.000,0.447,0.741;
        0.850,0.325,0.098;
        0.929,0.694,0.125;     
        0.494,0.184,0.556;
        0.466,0.674,0.188;
        0.301,0.745,0.933;
        0.635,0.078,0.184];
load twglaciers.mat; % load key data file

%% characteristic depths - fig. 2
if plotnum == 2,

lspace = 0.06;
rspace = 0.01;
bspace = 0.12;
tspace = 0.03;
hspace1 = 0.01;
hspace2 = 0.09;
pw1 = 0.52;
pw2 = 0.13;
pw3 = 1-lspace-rspace-hspace1-hspace2-pw1-pw2;
ph = 1-bspace-tspace;

% get data
for i=1:length(twglaciers),
    fjordvol(i) = twglaciers(i).fjord.vol;
    fjordmeandepth(i) = twglaciers(i).fjord.meandepth;
    fjordlength(i) = twglaciers(i).fjord.length;
    fjordwidth(i) = twglaciers(i).fjord.width;
    gldepth(i) = twglaciers(i).gldepth;
    effdepth(i) = twglaciers(i).fjord.effdepth;
    for j=1:7,
        nbdepth(i,j) = twglaciers(i).plume(j).znb;
        mhdepth(i,j) = twglaciers(i).plume(j).sol.z(end);
    end
end

% sorting by gldepth
[~,inds] = sort(gldepth);
gldepth = gldepth(inds);
effdepth = effdepth(inds);
for j=1:7,
    nbdepth(:,j) = nbdepth(inds,j);
    mhdepth(:,j) = mhdepth(inds,j);
end

% density profiles
z = [-899:0];
zextrap = -400;
sigma = [];
sigmaplume = [];
for i=1:length(twglaciers),
    Z = twglaciers(i).profile.z;
    T = twglaciers(i).profile.T;
    S = twglaciers(i).profile.S;
    % ambient
    dens1 = rho(T,S,0)-1000;
    if Z(1)<zextrap,
        dens2 = interp1(Z,dens1,z,'linear','extrap');
    else
        dens2 = interp1(Z,dens1,z,'linear',NaN);
    end
    sigma = [sigma;dens2];
    % plume
    dens1 = rho(twglaciers(i).plume(1).sol.T,twglaciers(i).plume(1).sol.S,0);
    dens2 = interp1(twglaciers(i).plume(1).sol.z,dens1,z,'linear',NaN);
    sigmaplume = [sigmaplume;dens2];
end

% density by sector
load ice_ocean_sectors.mat; % sector boundaries
for i=1:length(regions),
    regions(i).inds = find(inpolygon([twglaciers.x],[twglaciers.y],regions(i).ocean.x,regions(i).ocean.y));
    regions(i).sigmamean = nanmean(sigma(regions(i).inds,:));
    regions(i).sigmastd = nanstd(sigma(regions(i).inds,:));
end

% characteristic depths
ms = 5;
lw0 = 1;
a1 = axes('position',[lspace,bspace,pw1,ph]); hold on;
p(1)=plot(gldepth,'.-','linewidth',lw0,'markersize',ms,'markerfacecolor',cols(1,:));
p(2)=plot(effdepth,'.-','linewidth',lw0,'markersize',ms,'markerfacecolor',cols(2,:));
p(4)=plot(mhdepth(:,1),'.-','linewidth',lw0,'markersize',ms,'markerfacecolor',cols(3,:));
patch([1:length(twglaciers),fliplr([1:length(twglaciers)])],[min(mhdepth'),fliplr(max(mhdepth'))],cols(3,:),'facealpha',0.5,'edgecolor','none');
p(3)=plot(nbdepth(:,1),'.-','linewidth',lw0,'markersize',ms,'markerfacecolor',cols(4,:));
patch([1:length(twglaciers),fliplr([1:length(twglaciers)])],[min(nbdepth'),fliplr(max(nbdepth'))],cols(4,:),'facealpha',0.5,'edgecolor','none');
xlabel('tidewater glaciers, ordered by grounding line depth','fontsize',fs); ylabel('depth (m)','fontsize',fs);
set(gca,'box','on','fontsize',fs,'xtick',[0:25:150],'xgrid','off','ygrid','on');
legend(p,{'grounding line','fjord effective depth','plume neutral buoyancy','plume maximum height'},'location','southeast','fontsize',fs);
xlim([0 length(twglaciers)+1]); ylim([-1020 0]);
set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'});

% density profiles
a2 = axes('position',[lspace+pw1+hspace1,bspace,pw2,ph]); hold on;
sigmamean = nanmean(sigma);
sigmastd = nanstd(sigma);
sigmastd(1:find(z==zextrap)) = sigmastd(find(z==zextrap));
patch([sigmamean+sigmastd,fliplr(sigmamean-sigmastd)],...
    [z,fliplr(z)],'k','facealpha',0.25,'edgecolor','none');
plot(sigmamean,z,'k','linewidth',0.5);
plot(regions(1).sigmamean(find(z>-500)),z(find(z>-500)),'linewidth',1,'color',cols(5,:));
plot(regions(4).sigmamean(find(z>-500)),z(find(z>-500)),'linewidth',1,'color',cols(6,:));
plot(regions(6).sigmamean(find(z>-500)),z(find(z>-500)),'linewidth',1,'color',cols(7,:));
set(gca,'box','on','fontsize',fs,'ytick',[-900:100:0],'yticklabel',[]);
xlim([25 28]); ylim([-900,0]); grid on;
xlabel('$\sigma_{\theta}$ (kg$\,$m$^{-3}$)','fontsize',fs,'interpreter','latex');
% manual legend
tx = 25.2;
ty = -700;
dy = 50;
text(tx,ty-0*dy,'All','fontsize',fs,'color','k');
text(tx,ty-1*dy,'SE','fontsize',fs,'color',cols(5,:));
text(tx,ty-2*dy,'CW','fontsize',fs,'color',cols(6,:));
text(tx,ty-3*dy,'NW','fontsize',fs,'color',cols(7,:));

% uncertainty on neutral buoyancy for two glaciers
id = [3,5];
sep = 0.15;
a2 = axes('position',[lspace+pw1+hspace1+pw2+hspace2,bspace,pw3,ph]); hold on;
for i=1:length(id),
    % discharge
    yvals = [twglaciers(id(i)).plume(1).znb,twglaciers(id(i)).plume(4).znb,twglaciers(id(i)).plume(5).znb];
    errplot(a2,i-2*sep,yvals,cols(4,:),'p',0.5);
    plot(i-2*sep,twglaciers(id(i)).plume(1).znb,'o','color',cols(4,:),'markersize',4);
    % seasonal stratification
    yvals = [twglaciers(id(i)).seasonality.summer.nb,twglaciers(id(i)).seasonality.winter.nb];
    errplot(a2,i+sep,yvals,cols(4,:),'d',0.5);
    % interannual stratification
    yvals = [twglaciers(id(i)).plume(1).znb];
    for j=1:length(twglaciers(id(i)).interannual), yvals = [yvals,twglaciers(id(i)).interannual(j).nb]; end
    errplot(a2,i+2*sep,yvals,cols(4,:),'^',0.5);
    % hydrology
    yvals = [twglaciers(id(i)).plume(1).znb,twglaciers(id(i)).plume(2).znb,twglaciers(id(i)).plume(3).znb];
    errplot(a2,i-sep,yvals,cols(4,:),'s',0.5);
    plot(i-sep,twglaciers(id(i)).plume(1).znb,'o','color',cols(4,:),'markersize',4);
    % entrainment
    yvals = [twglaciers(id(i)).plume(1).znb,twglaciers(id(i)).plume(6).znb,twglaciers(id(i)).plume(7).znb];
    errplot(a2,i,yvals,cols(4,:),'x',0.5);
    plot(i,twglaciers(id(i)).plume(1).znb,'o','color',cols(4,:),'markersize',4);
end
set(gca,'box','on','fontsize',fs,'xtick',[1,2],'xticklabel',{'Helheim','Kangilliup'});
xtickangle(0);
set(gca,'ytick',[-400:100:-100],'yticklabel',{'400','300','200','100'});
% legend(q,{'discharge','hydrology'},'location','southeast','fontsize',fs);
xlim([0.5,2.5]); ylim([-350,-60]);
ylabel('plume neutral buoyancy depth (m)','fontsize',fs);
text(1-2*sep,-190,'subglacial discharge','fontsize',fs,'rotation',90,'verticalalignment','middle','horizontalalignment','right');
text(1,-190,'entrainment coefficient','fontsize',fs,'rotation',90,'verticalalignment','middle','horizontalalignment','right');
text(1+sep,-190,'seasonal stratification','fontsize',fs,'rotation',90,'verticalalignment','middle','horizontalalignment','right');
text(1+2*sep,-190,'interannual stratification','fontsize',fs,'rotation',90,'verticalalignment','middle','horizontalalignment','right');
text(1-sep,-190,'plume width','fontsize',fs,'rotation',90,'verticalalignment','middle','horizontalalignment','right');

% plot labels
annotation('textbox','position',[lspace,bspace+ph+0.01,0,0],'string','\textbf{a}','edgecolor','none','interpreter','latex','fontsize',fs+2);
annotation('textbox','position',[lspace+hspace1+pw1+pw2-0.03,bspace+ph+0.01,0,0],'string','\textbf{b}','edgecolor','none','interpreter','latex','fontsize',fs+2);
annotation('textbox','position',[0.96,bspace+ph+0.01,0,0],'string','\textbf{c}','edgecolor','none','interpreter','latex','fontsize',fs+2);

% save
saveplot(17,7,300,'fig2.png');
saveplot_pdf(17,7,300,'fig2.pdf');
close all;
    
end

%% characteristic depths - fig. 2 ordered clockwise
if plotnum == 2.1,

lspace = 0.06;
rspace = 0.02;
bspace = 0.12;
tspace = 0.03;
pw1 = 1-lspace-rspace;
ph = 1-bspace-tspace;

% get data
for i=1:length(twglaciers),
    x(i) = twglaciers(i).x;
    y(i) = twglaciers(i).y;
    fjordvol(i) = twglaciers(i).fjord.vol;
    fjordmeandepth(i) = twglaciers(i).fjord.meandepth;
    fjordlength(i) = twglaciers(i).fjord.length;
    fjordwidth(i) = twglaciers(i).fjord.width;
    gldepth(i) = twglaciers(i).gldepth;
    effdepth(i) = twglaciers(i).fjord.effdepth;
    for j=1:7,
        nbdepth(i,j) = twglaciers(i).plume(j).znb;
        mhdepth(i,j) = twglaciers(i).plume(j).sol.z(end);
    end
end

% sort clockwise around ice sheet from north
% separate into east and west
x1 = x(find(x>=0));
y1 = y(find(x>=0));
x2 = x(find(x<0));
y2 = y(find(x<0));
% sort east side from north to south
[~,inds] = sort(y1);
inds = fliplr(inds);
x1 = x1(inds); y1 = y1(inds);
% sort west side from south to north
[~,inds] = sort(y2);
x2 = x2(inds); y2 = y2(inds);
% put back together
x0 = [x1,x2];
y0 = [y1,y2];
% get final indices
clearvars inds;
for i=1:length(x0),
    inds(i) = find(x==x0(i) & y==y0(i));
end
% sort data
x = x(inds);
y = y(inds);
gldepth = gldepth(inds);
effdepth = effdepth(inds);
for j=1:7,
    nbdepth(:,j) = nbdepth(inds,j);
    mhdepth(:,j) = mhdepth(inds,j);
end

% on this basis, the indices correspond to sectors as follows
% NO: 132 -> 3
% NE: 4 -> 14
% CE: 15 -> 25
% SE: 26 -> 59
% SW: 60 -> 62
% CW: 63 -> 77
% NW: 78 -> 131
xt = [3.5,14.5,25.5,59.5,62.5,77.5,131.5];

% characteristic depths
ms = 5;
lw0 = 1;
a1 = axes('position',[lspace,bspace,pw1,ph]); hold on;
p(1)=plot(gldepth,'.-','linewidth',lw0,'markersize',ms,'markerfacecolor',cols(1,:));
p(2)=plot(effdepth,'.-','linewidth',lw0,'markersize',ms,'markerfacecolor',cols(2,:));
p(4)=plot(mhdepth(:,1),'.-','linewidth',lw0,'markersize',ms,'markerfacecolor',cols(3,:));
patch([1:length(twglaciers),fliplr([1:length(twglaciers)])],[min(mhdepth'),fliplr(max(mhdepth'))],cols(3,:),'facealpha',0.5,'edgecolor','none');
p(3)=plot(nbdepth(:,1),'.-','linewidth',lw0,'markersize',ms,'markerfacecolor',cols(4,:));
patch([1:length(twglaciers),fliplr([1:length(twglaciers)])],[min(nbdepth'),fliplr(max(nbdepth'))],cols(4,:),'facealpha',0.5,'edgecolor','none');
xlabel('tidewater glaciers, ordered clockwise from north','fontsize',fs); ylabel('depth (m)','fontsize',fs);
set(gca,'box','on','fontsize',fs,'xtick',[0:25:150],'xgrid','off','ygrid','on');
legend(p,{'grounding line','fjord effective depth','plume neutral buoyancy','plume maximum height'},'location','southeast','fontsize',fs);
xlim([0 length(twglaciers)+1]); ylim([-1020 0]);
set(gca,'ytick',[-1000:200:0],'yticklabel',{'1000','800','600','400','200','0'});
set(gca,'xtick',xt,'xticklabel',{' '});
text(1.5,-1060,'NO','fontsize',fs,'HorizontalAlignment','center');
text(0.5*(xt(1)+xt(2)),-1060,'NE','fontsize',fs,'HorizontalAlignment','center');
text(0.5*(xt(2)+xt(3)),-1060,'CE','fontsize',fs,'HorizontalAlignment','center');
text(0.5*(xt(3)+xt(4)),-1060,'SE','fontsize',fs,'HorizontalAlignment','center');
text(0.5*(xt(4)+xt(5)),-1060,'SW','fontsize',fs,'HorizontalAlignment','center');
text(0.5*(xt(5)+xt(6)),-1060,'CW','fontsize',fs,'HorizontalAlignment','center');
text(0.5*(xt(6)+xt(7)),-1060,'NW','fontsize',fs,'HorizontalAlignment','center');
text(134,-1060,'NO','fontsize',fs,'HorizontalAlignment','center');

% save
saveplot(17,7,300,'figS4.png');
close all;
    
end

%% fluxes (including dilution factors) and renewal times - fig. 3
if plotnum == 3,

lspace = 0.07;
rspace = 0.01;
hspace1 = 0.09;
hspace2 = 0.04;
bspace1 = 0.1;
bspace2 = 0.15;
vspace = 0.1;
tspace = 0.025;
pw1 = 0.53;
pw2 = (1-lspace-rspace-hspace1-hspace2-pw1)/2;
ph1 = (1-tspace-bspace1-vspace)/2;
ph2 = 1-bspace2-tspace;
    
% get data
for i=1:length(twglaciers),
    Q0(i) = twglaciers(i).Qsg0;
    Qmax(i) = twglaciers(i).Qsgmin;
    Qmin(i) = twglaciers(i).Qsgmax;
    gldepth(i) = twglaciers(i).gldepth;
    Qnb(i) = twglaciers(i).plume(1).Qnb;
    znb(i) = twglaciers(i).plume(1).znb;
    V(i) = twglaciers(i).fjord.vol_below_NB;
    T(i) = twglaciers(i).fjord.renewaltime;
    Tmin(i) = twglaciers(i).fjord.renewaltimemin;
    Tmax(i) = twglaciers(i).fjord.renewaltimemax;
    for j=1:5
        Qnball(i,j) = twglaciers(i).plume(j).Qnb;
    end
end

% put min 1 day renewal timescale
T(T<1) = 1;
Tmin(Tmin<1) = 1;
Tmax(Tmax<1) = 1;

% sorting by gldepth
[~,inds] = sort(gldepth);
gldepth = gldepth(inds);
Q0 = Q0(inds);
Qmin = Qmin(inds);
Qmax = Qmax(inds);
Qnb = Qnb(inds);
V = V(inds);
T = T(inds);
Tmin = Tmin(inds);
Tmax = Tmax(inds);
for j=1:5
    Qnball(:,j) = Qnball(inds,j);
end

% binning by depths
dbins = [-1020:10:0];
dbinshalf = dbins+0.5*(dbins(2)-dbins(1));
Q0bins = 0*dbins;
Qnbbins = 0*dbins;
for i=1:length(dbins)-1,
    inds = find(gldepth>=dbins(i) & gldepth<dbins(i+1));
    Q0bins(i) = sum(Q0(inds))/(dbins(2)-dbins(1));
    inds = find(znb>=dbins(i) & znb<dbins(i+1));
    Qnbbins(i) = sum(Qnb(inds))/(dbins(2)-dbins(1));
end
% make vector with correct abundance so that can fit kernel
glvec = [];
nbvec = [];
for i=1:length(dbins),
    glvec = [glvec,dbinshalf(i)*ones(1,round(Q0bins(i)))];
    nbvec = [nbvec,dbinshalf(i)*ones(1,round(Qnbbins(i)))];
end
    
% find those with renewal in less than 1 season
inds = find(T<=100);

figure();
ms = 5;
lw0 = 1;

% fluxes
a1 = axes('position',[lspace,bspace1+1*(ph1+vspace),pw1,ph1]); hold on;
% sgd
patch([1:length(twglaciers),fliplr([1:length(twglaciers)])],[Qmin,fliplr(Qmax)]/10^3,cols(1,:),'facealpha',0.5,'edgecolor','none');
plot(Q0/10^3,'.-','color',cols(1,:),'linewidth',lw0,'markersize',ms,'markerfacecolor',cols(1,:));
patch([1:length(twglaciers),fliplr([1:length(twglaciers)])],[min(Qnball'),fliplr(max(Qnball'))]/10^3,cols(4,:),'facealpha',0.5,'edgecolor','none');
plot(Qnb/10^3,'.-','color',cols(4,:),'linewidth',lw0,'markersize',ms,'markerfacecolor',cols(4,:));
ylabel('flux (mSv = 10$^3$ m$^3$/s)','fontsize',fs);
set(gca,'box','on','fontsize',fs,'xtick',[0:25:150],'yscale','log',...
    'ytick',logspace(-2,2,5),'yticklabel',{'0.01','0.1','1','10','100'});
xlim([0 length(twglaciers)+1]); ylim([0.01 100]);
text(30,1,'\textbf{subglacial discharge}','fontsize',fs,'color',cols(1,:));
text(70,25,'\textbf{upwelling flux}','fontsize',fs,'color',cols(4,:));

% discharge versus grounding line depth
a2 = axes('position',[lspace+pw1+hspace1,bspace2,pw2,ph2]); hold on;
fa = 0.3;
pd1 = fitdist(glvec','kernel','bandwidth',bw);
xvals = [-1020:0];
yvals1 = pdf(pd1,xvals);
% renormalize to get correct weighting (factors of 10^3 result in mSv)
yvals1 = yvals1*sum(Q0)/trapz(xvals,yvals1);
patch([0,yvals1,0]/10^3,[-1020,xvals,0],cols(1,:),'edgecolor','none','facealpha',fa);
plot(yvals1/10^3,xvals,'color',cols(1,:),'linewidth',1);
ylim([-1020,0]); xlim([0 0.04]); grid on;
set(gca,'box','on','fontsize',fs,'ytick',[-1000:200:0],...
    'yticklabel',{'1000','800','600','400','200','0'},...
    'xtick',[0:0.02:0.06]);
xlabel({'subglacial';'discharge (mSv/m)'},'fontsize',fs);
ylabel('grounding line depth (m)','fontsize',fs);
text(0.005,-350,['\textbf{',num2str(0.01*round(100*sum(Q0(:))/10^6)),' Sv}'],...
    'fontsize',fs,'color',cols(1,:));

% Qnb versus neutral buoyancy depth
a3 = axes('position',[lspace+pw1+hspace1+pw2+hspace2,bspace2,pw2,ph2]); hold on;
fa = 0.3;
pd2 = fitdist(nbvec','kernel','bandwidth',bw);
xvals = [-1020:0];
yvals2 = pdf(pd2,xvals);
% renormalize to get correct weighting (factors of 10^3 result in mSv)
yvals2 = yvals2*sum(Qnb)/trapz(xvals,yvals2);
patch([0,yvals2,0]/10^3,[-1020,xvals,0],cols(4,:),'edgecolor','none','facealpha',fa);
plot(yvals2/10^3,xvals,'color',cols(4,:),'linewidth',1);
ylim([-1020,0]); xlim([0 6.5]); grid on;
set(gca,'box','on','fontsize',fs,'yticklabel',[],'xtick',[0:2:6]);
xlabel({'flux after';'upwelling (mSv/m)'},'fontsize',fs);
ylabel('plume neutral buoyancy depth (m)','fontsize',fs);
text(1.3,-100,['\textbf{',num2str(0.01*round(100*sum(Qnb(:))/10^6)),' Sv}'],...
    'fontsize',fs,'color',cols(4,:));

% renewal timescale
a4 = axes('position',[lspace,bspace1+0*(ph1+vspace),pw1,ph1]); hold on;
plot([0,length(twglaciers)+1],100*[1,1],'--','color',0.5*[1,1,1],'linewidth',1);
patch([1:length(twglaciers),fliplr([1:length(twglaciers)])],[Tmin,fliplr(Tmax)],'k','facealpha',0.5,'edgecolor','none');
plot(T,'k.-','linewidth',lw0,'markersize',ms,'markerfacecolor','k');
ylabel('plume renewal time (days)','fontsize',fs);
set(gca,'box','on','fontsize',fs,'xtick',[0:25:150],'yscale','log',...
    'ytick',[10^0,10^1,10^2,10^3],'yticklabel',{'1','10','100','1000'});
xlim([0 length(twglaciers)+1]); ylim([0.8,3*10^3]);
xlabel('tidewater glaciers, ordered by grounding line depth','fontsize',fs);

% plot labels
annotation('textbox','position',[lspace+pw1-0.03,bspace1+1*(ph1+vspace)+ph1,0,0],'string','\textbf{a}','edgecolor','none','interpreter','latex','fontsize',fs+2);
annotation('textbox','position',[lspace+pw1+hspace1+0.1,bspace2+ph2,0,0],'string','\textbf{b}','edgecolor','none','interpreter','latex','fontsize',fs+2);
annotation('textbox','position',[lspace+pw1+hspace1+pw2+hspace2+0.1,bspace2+ph2,0,0],'string','\textbf{c}','edgecolor','none','interpreter','latex','fontsize',fs+2);
annotation('textbox','position',[lspace,bspace1+0*(ph1+vspace)+ph1,0,0],'string','\textbf{d}','edgecolor','none','interpreter','latex','fontsize',fs+2);

saveplot(17,8.5,300,'fig3_old.png');
close all;
    
    
end

%% fig. 3 with sensitivity panel
if plotnum == 3.1,

lspace = 0.07;
rspace = 0.01;
hspace1 = 0.08;
hspace2 = 0.04;
bspace1 = 0.1;
bspace2 = 0.5;
vspace = 0.1;
tspace = 0.025;
pw1 = 0.54;
pw2 = (1-lspace-rspace-hspace1-hspace2-pw1)/2;
ph1 = (1-tspace-bspace1-vspace)/2;
ph2 = 1-bspace2-tspace;
    
% get data
for i=1:length(twglaciers),
    Q0(i) = twglaciers(i).Qsg0;
    Qmax(i) = twglaciers(i).Qsgmin;
    Qmin(i) = twglaciers(i).Qsgmax;
    gldepth(i) = twglaciers(i).gldepth;
    Qnb(i) = twglaciers(i).plume(1).Qnb;
    znb(i) = twglaciers(i).plume(1).znb;
    V(i) = twglaciers(i).fjord.vol_below_NB;
    T(i) = twglaciers(i).fjord.renewaltime;
    Tmin(i) = twglaciers(i).fjord.renewaltimemin;
    Tmax(i) = twglaciers(i).fjord.renewaltimemax;
    for j=1:7,
        Qnball(i,j) = twglaciers(i).plume(j).Qnb;
    end
end

% put min 1 day renewal timescale
T(T<1) = 1;
Tmin(Tmin<1) = 1;
Tmax(Tmax<1) = 1;

% sorting by gldepth
[~,inds] = sort(gldepth);
gldepth = gldepth(inds);
Q0 = Q0(inds);
Qmin = Qmin(inds);
Qmax = Qmax(inds);
Qnb = Qnb(inds);
V = V(inds);
T = T(inds);
Tmin = Tmin(inds);
Tmax = Tmax(inds);
for j=1:7,
    Qnball(:,j) = Qnball(inds,j);
end

% binning by depths
dbins = [-1020:10:0];
dbinshalf = dbins+0.5*(dbins(2)-dbins(1));
Q0bins = 0*dbins;
Qnbbins = 0*dbins;
for i=1:length(dbins)-1,
    inds = find(gldepth>=dbins(i) & gldepth<dbins(i+1));
    Q0bins(i) = sum(Q0(inds))/(dbins(2)-dbins(1));
    inds = find(znb>=dbins(i) & znb<dbins(i+1));
    Qnbbins(i) = sum(Qnb(inds))/(dbins(2)-dbins(1));
end
% make vector with correct abundance so that can fit kernel
glvec = [];
nbvec = [];
for i=1:length(dbins),
    glvec = [glvec,dbinshalf(i)*ones(1,round(Q0bins(i)))];
    nbvec = [nbvec,dbinshalf(i)*ones(1,round(Qnbbins(i)))];
end
    
% find those with renewal in less than 1 season
inds = find(T<=100);

figure();
ms = 5;
lw0 = 1;

% fluxes
a1 = axes('position',[lspace,bspace1+1*(ph1+vspace),pw1,ph1]); hold on;
% sgd
patch([1:length(twglaciers),fliplr([1:length(twglaciers)])],[Qmin,fliplr(Qmax)]/10^3,cols(1,:),'facealpha',0.5,'edgecolor','none');
plot(Q0/10^3,'.-','color',cols(1,:),'linewidth',lw0,'markersize',ms,'markerfacecolor',cols(1,:));
patch([1:length(twglaciers),fliplr([1:length(twglaciers)])],[min(Qnball'),fliplr(max(Qnball'))]/10^3,cols(4,:),'facealpha',0.5,'edgecolor','none');
plot(Qnb/10^3,'.-','color',cols(4,:),'linewidth',lw0,'markersize',ms,'markerfacecolor',cols(4,:));
ylabel('flux (mSv = 10$^3$ m$^3$/s)','fontsize',fs);
set(gca,'box','on','fontsize',fs,'xtick',[0:25:150],'yscale','log',...
    'ytick',logspace(-2,2,5),'yticklabel',{'0.01','0.1','1','10','100'});
xlim([0 length(twglaciers)+1]); ylim([0.01 100]);
text(30,1,'\textbf{subglacial discharge}','fontsize',fs,'color',cols(1,:));
text(70,25,'\textbf{upwelling flux}','fontsize',fs,'color',cols(4,:));

% discharge versus grounding line depth
a2 = axes('position',[lspace+pw1+hspace1,bspace2,pw2,ph2]); hold on;
fa = 0.3;
pd1 = fitdist(glvec','kernel','bandwidth',bw);
xvals = [-1020:0];
yvals1 = pdf(pd1,xvals);
% renormalize to get correct weighting (factors of 10^3 result in mSv)
yvals1 = yvals1*sum(Q0)/trapz(xvals,yvals1);
patch([0,yvals1,0]/10^3,[-1020,xvals,0],cols(1,:),'edgecolor','none','facealpha',fa);
plot(yvals1/10^3,xvals,'color',cols(1,:),'linewidth',1);
ylim([-1020,0]); xlim([0 0.04]); grid on;
set(gca,'box','on','fontsize',fs,'ytick',[-1000:200:0],...
    'yticklabel',{'1000','800','600','400','200','0'},...
    'xtick',[0:0.02:0.06]); xtickangle(0);
xlabel({'subglacial';'discharge (mSv/m)'},'fontsize',fs);
ylabel('grounding line depth (m)','fontsize',fs);
text(0.005,-350,['\textbf{',num2str(0.01*round(100*sum(Q0(:))/10^6)),' Sv}'],...
    'fontsize',fs,'color',cols(1,:));

% Qnb versus neutral buoyancy depth
a3 = axes('position',[lspace+pw1+hspace1+pw2+hspace2,bspace2,pw2,ph2]); hold on;
fa = 0.3;
pd2 = fitdist(nbvec','kernel','bandwidth',bw);
xvals = [-1020:0];
yvals2 = pdf(pd2,xvals);
% renormalize to get correct weighting (factors of 10^3 result in mSv)
yvals2 = yvals2*sum(Qnb)/trapz(xvals,yvals2);
patch([0,yvals2,0]/10^3,[-1020,xvals,0],cols(4,:),'edgecolor','none','facealpha',fa);
plot(yvals2/10^3,xvals,'color',cols(4,:),'linewidth',1);
ylim([-1020,0]); xlim([0 6.5]); grid on;
set(gca,'box','on','fontsize',fs,'yticklabel',[],'xtick',[0:2:6]);
xlabel({'flux after';'upwelling (mSv/m)'},'fontsize',fs);
ylabel('plume neutral buoyancy depth (m)','fontsize',fs);
text(1.1,-110,['\textbf{',num2str(0.01*round(100*sum(Qnb(:))/10^6)),' Sv}'],...
    'fontsize',fs,'color',cols(4,:));

% renewal timescale
a4 = axes('position',[lspace,bspace1+0*(ph1+vspace),pw1,ph1]); hold on;
plot([0,length(twglaciers)+1],100*[1,1],'--','color',0.5*[1,1,1],'linewidth',1);
patch([1:length(twglaciers),fliplr([1:length(twglaciers)])],[Tmin,fliplr(Tmax)],'k','facealpha',0.5,'edgecolor','none');
plot(T,'k.-','linewidth',lw0,'markersize',ms,'markerfacecolor','k');
ylabel('plume-driven renewal time (days)','fontsize',fs);
set(gca,'box','on','fontsize',fs,'xtick',[0:25:150],'yscale','log',...
    'ytick',[10^0,10^1,10^2,10^3],'yticklabel',{'1','10','100','1000'});
xlim([0 length(twglaciers)+1]); ylim([0.8,3*10^3]);
xlabel('tidewater glaciers, ordered by grounding line depth','fontsize',fs);

% total upwelling sensitivity panel
ms = 6;
a5 = axes('position',[lspace+pw1+hspace1,bspace1,2*pw2+hspace2,0.23]); hold on;
xscaling = linspace(0,3,100);
yscaling = sum(Qnball(:,1))*xscaling.^(2/3)/1e6;
qq(3) = plot(xscaling,yscaling,'--','color',0.5*[1,1,1],'linewidth',0.75);
plot(1,sum(Qnball(:,1))/1e6,'ks','markersize',ms);
plot(1,sum(Qnball(:,1))/1e6,'kx','markersize',ms);
qq(1) = plot(100/250,sum(Qnball(:,2))/1e6,'ks','markersize',ms);
plot(500/250,sum(Qnball(:,3))/1e6,'ks','markersize',ms);
qq(2) = plot(0.05/0.1,sum(Qnball(:,6))/1e6,'kx','markersize',ms);
plot(0.15/0.1,sum(Qnball(:,7))/1e6,'kx','markersize',ms);
set(gca,'box','on','fontsize',fs);
xlabel('ratio of parameter to default','fontsize',fs);
ylabel('total upwelling (Sv)','fontsize',fs);
xlim([0.3,2.1]); ylim([0.5,1.8]);
% custom legend
customlegend(a5,qq,1.42,0.65,0.1,0.2,0.15,...
    {'plume width','entrainment','scaling'},fs);

% plot labels
annotation('textbox','position',[lspace+pw1-0.03,bspace1+1*(ph1+vspace)+ph1,0,0],'string','\textbf{a}','edgecolor','none','interpreter','latex','fontsize',fs+2);
annotation('textbox','position',[lspace+pw1+hspace1+0.1,bspace2+ph2,0,0],'string','\textbf{b}','edgecolor','none','interpreter','latex','fontsize',fs+2);
annotation('textbox','position',[lspace+pw1+hspace1+pw2+hspace2+0.102,bspace2+ph2+0.01,0,0],'string','\textbf{c}','edgecolor','none','interpreter','latex','fontsize',fs+2);
annotation('textbox','position',[lspace,bspace1+0*(ph1+vspace)+ph1,0,0],'string','\textbf{e}','edgecolor','none','interpreter','latex','fontsize',fs+2);
annotation('textbox','position',[lspace+pw1+hspace1,bspace1+0.23,0,0],'string','\textbf{d}','edgecolor','none','interpreter','latex','fontsize',fs+2);

saveplot(17,8.5,300,'fig3.png');
saveplot_pdf(17,8.5,300,'fig3.pdf');
close all;
    
    
end



%% regional differences plot - fig. 4
if plotnum == 4,

% ice-ocean sectors    
load ice_ocean_sectors.mat
names = {'SE','SW','CE','CW','NE','NW','NO'};
for i=1:length(regions),
    regions(i).inds = find(inpolygon([twglaciers.x],[twglaciers.y],regions(i).ocean.x,regions(i).ocean.y));
end

for i=1:length(twglaciers),
    gldepth(i) = twglaciers(i).gldepth; % GL depth
    fedepth(i) = twglaciers(i).fjord.effdepth; % fe depth
    nbdepth(i) = twglaciers(i).plume(1).znb; % NB depth
    sgd(i) = twglaciers(i).Qsg0; % SGD
    upwellingflux(i) = twglaciers(i).plume(1).Qnb; % upwelled flux
    fjordlength(i) = twglaciers(i).fjord.length; % length
    renewaltime(i) = twglaciers(i).fjord.renewaltime; % renewal time
end

for i=1:length(regions),
    regions(i).gldepth = gldepth(regions(i).inds);
    regions(i).fedepth = fedepth(regions(i).inds);
    regions(i).nbdepth = nbdepth(regions(i).inds);
    regions(i).q = sgd(regions(i).inds);
    regions(i).qnb = upwellingflux(regions(i).inds);
end
    

N = [];
G = [];
F = [];
Z = [];
L = [];
T = [];
sectororder = [5,3,1,2,4,6,7];
for i=1:length(sectororder),
    j = sectororder(i);
    N = [N;repmat(names(j),length(gldepth(regions(j).inds)),1)];
    G = [G;gldepth(regions(j).inds)'];
    F = [F;fedepth(regions(j).inds)'];
    Z = [Z;nbdepth(regions(j).inds)'];
    L = [L;fjordlength(regions(j).inds)'];
    T = [T;renewaltime(regions(j).inds)'];
    Q0(j) = sum(sgd(regions(j).inds));
    Qnb(j) = sum(upwellingflux(regions(j).inds));
    Gm(j) = median(gldepth(regions(j).inds));
    Zm(j) = median(nbdepth(regions(j).inds));
end
T(T<1) = 1;
    
lspace1 = 0.065;
lspace2 = 0.08;
rspace1 = 0.01;
rspace2 = 0.09;
hspace = 0.14;
bspace = 0.05;
tspace = 0.04;
vspace = 0.1;
pw1 = 1-lspace1-rspace1;
pw2 = (1-lspace2-rspace2-hspace)/2;
ph = (1-tspace-bspace-vspace)/2;
figure();

% characteristic depths
a1 = axes('position',[lspace1,bspace+1*(ph+vspace),pw1,ph]); hold on;
w = 0.2;
lw = 1;
dw = 0.25;
for i=1:length(sectororder),
    j = sectororder(i);
    plotbox(a1,i-dw,regions(j).gldepth,w,cols(1,:),lw);
    plotbox(a1,i,regions(j).fedepth,w,cols(2,:),lw);
    plotbox(a1,i+dw,regions(j).nbdepth,w,cols(4,:),lw);
end

set(gca,'box','on','fontsize',fs,'ticklabelinterpreter','latex',...
    'ytick',[-1000:200:0],'yticklabel',...
    {'1000','800','600','400','200','0'});
xstring = names(sectororder);
for i=1:length(xstring),
    j = sectororder(i);
    xstring{i} = [xstring{i},' (',num2str(length(regions(j).inds)),')'];
end
set(gca,'xtick',[1:7],'xticklabel',xstring);
ylabel('depth (m)','fontsize',fs);
xlim([1-0.5,7+0.5]); ylim([-1020,0]);
tx = 1;
ty = 50;
dx = 1;
text(3.1,ty,'grounding line depth','fontsize',fs,'interpreter','latex','color',cols(1,:),...
    'horizontalalignment','left');
text(4.5,ty,'fjord effective depth','fontsize',fs,'interpreter','latex','color',cols(2,:),...
    'horizontalalignment','left');
text(5.8,ty,'plume neutral buoyancy depth','fontsize',fs,'interpreter','latex','color',cols(4,:),...
    'horizontalalignment','left');

% subglacial discharge
a2 = axes('position',[lspace2,bspace,pw2,ph]); hold on;
% yyaxis left
for i=1:length(sectororder),
    j = sectororder(i);
    plotbox(a2,i,regions(j).q/1000,0.5,cols(1,:),lw);
end
ylabel({'subglacial discharge','by glacier (mSv)'},'fontsize',fs);
set(gca,'box','on','fontsize',fs,'ticklabelinterpreter','latex');
set(gca,'xtick',[1:7],'xticklabel',names(sectororder));
xlim([0.5,7.5]); ylim([0,1.3]);

% total fluxes
a3 = axes('position',[lspace2+1*(pw2+hspace),bspace,pw2,ph]); hold on;
yyaxis left
bar([1:7]-1/6,Q0(sectororder)/1000,1/3,'facecolor',cols(1,:));
ylabel({'subglacial discharge','by region (mSv)'},'fontsize',fs);
yyaxis right
bar([1:7]+1/6,Qnb(sectororder)/1000,1/3,'facecolor',cols(4,:));
set(gca,'Ycolor',cols(4,:));
ylabel({'upwelling flux','by region (mSv)'},'fontsize',fs);
set(gca,'box','on','fontsize',fs,'ticklabelinterpreter','latex');
xlim([0.5,7.5]);
set(gca,'xtick',[1:7],'xticklabel',names(sectororder));

% plot labels
dx = -0.05;
dy = 0.04;
annotation('textbox','position',[lspace1+dx,bspace+2*ph+vspace+dy,0,0],'string','\textbf{a}','edgecolor','none','interpreter','latex','fontsize',fs+2);
annotation('textbox','position',[lspace2+dx-0.01,bspace+ph+dy,0,0],'string','\textbf{b}','edgecolor','none','interpreter','latex','fontsize',fs+2);
annotation('textbox','position',[lspace2+1*(pw2+hspace)+dx,bspace+ph+dy,0,0],'string','\textbf{c}','edgecolor','none','interpreter','latex','fontsize',fs+2);

saveplot(17,10,300,'fig4.png');
saveplot_pdf(17,10,300,'fig4.pdf');
close all;
    
    
end

%% distance of casts from fronts and timings - fig. S1
if plotnum == 5,
    
    lspace = 0.07;
    hspace = 0.11;
    rspace = 0.02;
    pw1 = 0.32;
    pw2 = 1-lspace-rspace-hspace-pw1;
    bspace = 0.17;
    tspace = 0.03;
    ph = 1-bspace-tspace;
    
    
    % distances
    x = [twglaciers.x]; y = [twglaciers.y];
    for i=1:length(twglaciers),
        xp(i) = twglaciers(i).profile.x;
        yp(i) = twglaciers(i).profile.y;
    end
    d = sqrt((x-xp).^2+(y-yp).^2)/10^3;
    figure();
    axes('position',[lspace,bspace,pw1,ph]);
    histogram(d,[0:5:150],'facecolor',cols(1,:));
    set(gca,'box','on','fontsize',fs);
    xlabel('distance from CTD cast to glacier (km)','fontsize',fs);
    ylabel('number of glaciers/casts','fontsize',fs);

    % timings
    for i=1:length(twglaciers),
        tstr = datevec(twglaciers(i).profile.t);
        yrnum(i) = tstr(1);
        mnthnum(i) = tstr(2);
    end
    
    % years span 2013-2020, months span July-October
    yrunique = unique(yrnum);
    mnthunique = unique(mnthnum);
    % add a dummy month to make nice gap in bar graph
    mnthunique = [mnthunique,20];
    % make useful vector for bar graph
    tvec = zeros(1,length(yrunique)*length(mnthunique));
    for i=1:length(yrunique),
        for j=1:length(mnthunique),
            tvec((i-1)*length(mnthunique)+j) = sum(yrnum==yrunique(i) & mnthnum==mnthunique(j));
            if mnthunique(j)<10,
                lab{(i-1)*length(mnthunique)+j} = [num2str(yrunique(i)),'-0',num2str(mnthunique(j))];
            elseif mnthunique(j)<13,
                lab{(i-1)*length(mnthunique)+j} = [num2str(yrunique(i)),'-',num2str(mnthunique(j))];
            else
                lab{(i-1)*length(mnthunique)+j} = '';
            end
        end
    end
    
    axes('position',[lspace+pw1+hspace,bspace,pw2,ph]);
    bar(tvec,'facecolor',cols(1,:));
    set(gca,'xtick',[1:length(tvec)],'xticklabel',lab);
    xtickangle(gca,90);
    set(gca,'fontsize',fs);
    ylabel('number of glaciers/casts','fontsize',fs);
    xlim([0,length(tvec)-1]);
    
    % plot labels
    annotation('textbox','position',[lspace+0*(pw1+hspace)-0.07,bspace+ph+0.02,0,0],'string','\textbf{a}','edgecolor','none','interpreter','latex','fontsize',fs+4);
    annotation('textbox','position',[lspace+1*(pw1+hspace)-0.07,bspace+ph+0.02,0,0],'string','\textbf{b}','edgecolor','none','interpreter','latex','fontsize',fs+4);
    
    saveplot(15,7,300,'figS1.png');
    close all;  
    
end

%% kolmogorov-smirnov tests on regions vs all - fig. S6
if plotnum == 9,
    
load ice_ocean_sectors.mat
names = {'SE','SW','CE','CW','NE','NW','NO'};

% reorder into required order
sectororder = [5,3,1,2,4,6,7];
names = names(sectororder);
regions2 = regions;
for i=1:length(sectororder),
    regions2(i) = regions(sectororder(i));
end
regions = regions2;

for i=1:length(regions),
    regions(i).inds = find(inpolygon([twglaciers.x],[twglaciers.y],regions(i).ocean.x,regions(i).ocean.y));
end

for i=1:length(twglaciers),
    gldepth(i) = twglaciers(i).gldepth; % GL depth
    fedepth(i) = twglaciers(i).fjord.effdepth; % fe depth
    nbdepth(i) = twglaciers(i).plume(1).znb; % NB depth
    sgd(i) = twglaciers(i).Qsg0; % SGD
    upwellingflux(i) = twglaciers(i).plume(1).Qnb; % upwelled flux
    fjordlength(i) = twglaciers(i).fjord.length; % length
    renewaltime(i) = twglaciers(i).fjord.renewaltime; % renewal time
end

for i=1:length(regions),
    regions(i).gldepth = gldepth(regions(i).inds);
    regions(i).fedepth = fedepth(regions(i).inds);
    regions(i).nbdepth = nbdepth(regions(i).inds);
    regions(i).renewaltime = renewaltime(regions(i).inds);
    regions(i).sgd = sgd(regions(i).inds);
    regions(i).upw = upwellingflux(regions(i).inds);
end

for i=1:length(regions),
    glallother = [];
    feallother = [];
    nballother = [];
    renewalallother = [];
    sgdallother = [];
    upwallother = [];
    for j=1:length(regions),
        if i~=j, glallother = [glallother,regions(j).gldepth]; end
        if i~=j, feallother = [feallother,regions(j).fedepth]; end
        if i~=j, nballother = [nballother,regions(j).nbdepth]; end
        if i~=j, renewalallother = [renewalallother,regions(j).renewaltime]; end
        if i~=j, sgdallother = [sgdallother,regions(j).sgd]; end
        if i~=j, upwallother = [upwallother,regions(j).upw]; end
    end
    [~,p(1,i)] = kstest2(regions(i).gldepth,glallother);
    [~,p(2,i)] = kstest2(regions(i).fedepth,feallother);
    [~,p(3,i)] = kstest2(regions(i).nbdepth,nballother);
    [~,p(4,i)] = kstest2(regions(i).renewaltime,renewalallother);
    [~,p(5,i)] = kstest2(regions(i).sgd,sgdallother);
    [~,p(6,i)] = kstest2(regions(i).upw,upwallother);
end

reds = cbrewer('seq','Reds',9);
reds = [1,1,1;reds];
reds = flipud(reds);
lspace = 0.24;
hspace = 0.07;
rspace = 0.1;
pw = 1-lspace-rspace;
bspace = 0.08;
tspace = 0.04;
ph = 1-bspace-tspace;
cbspace = 0.02;
cbw = 0.02;
cbh = 0.6;
fs2 = 6;

figure();
a1 = axes('position',[lspace,bspace,pw,ph]);
imagesc(p); colormap(reds); caxis([0 1]);
axis xy;
set(gca,'xtick',[1:7],'xticklabel',names,'fontsize',fs);
set(gca,'ytick',[1:6],'yticklabel',...
    {'grounding line depth','fjord effective depth','plume neutral buoyancy depth',...
    'plume renewal time','subglacial discharge','upwelling flux'});
for i=1:length(regions),
    for j=1:6,
        if p(j,i)>=0.05,
            text(i,j,num2str(0.01*round(100*p(j,i))),'fontsize',fs2,...
                'horizontalalignment','center','verticalalignment','middle');
        else
            text(i,j,num2str(0.01*round(100*p(j,i))),'fontsize',fs2,'color','y',...
                'horizontalalignment','center','verticalalignment','middle');
        end
    end
end

h = colorbar('position',[lspace+pw+cbspace,bspace+(ph-cbh)/2,cbw,cbh]);
set(h,'fontsize',fs,'ticklabelinterpreter','latex');
ylabel(h,'$p$-value','fontsize',fs,'interpreter','latex');

saveplot(17,10,300,'figS6.png');
close all;

end

%% neutral buoyancy sensitivity plots - fig. S5
if plotnum == 8,
    
    lspace = 0.065;
    rspace = 0.015;
    hspace = 0.04;
    tspace = 0.05;
    bspace = 0.11;
    pw = (1-lspace-rspace-3*hspace)/4;
    ph = 1-bspace-tspace;
    lw1 = 1;
    lw2 = 0.5;
    
    figure();
    
    gid = 3;
    a1 = axes('position',[lspace+0*(pw+hspace),bspace,pw,ph]); hold on;
    plot([25,28],twglaciers(gid).gldepth*[1,1],'linewidth',lw1,'color',0.5*[1,1,1]);
    for i=1:length(twglaciers(gid).interannual),
        q(i) = plot(twglaciers(gid).interannual(i).sta,twglaciers(gid).interannual(i).z,'linewidth',lw1,'color',cols(i,:));
        plot(twglaciers(gid).interannual(i).st,twglaciers(gid).interannual(i).z,'linewidth',lw2,'color',cols(i,:));
        plot([25,28],twglaciers(gid).interannual(i).nb*[1,1],'--','linewidth',lw2,'color',cols(i,:));
    end
%     z = twglaciers(gid).plume(1).sol.z;
%     sta = rho(twglaciers(gid).plume(1).sol.Ta,twglaciers(gid).plume(1).sol.Sa,0)-1000;
%     plot(sta,z,'k','linewidth',lw1);
%     st = rho(twglaciers(gid).plume(1).sol.T,twglaciers(gid).plume(1).sol.S,0)-1000;
%     plot(st,z,'linewidth',lw2,'color','k');
%     plot([25,28],twglaciers(gid).plume(1).znb*[1,1],'--','linewidth',lw2,'color','k');
    q(length(twglaciers(gid).interannual)+1) = plot(twglaciers(gid).seasonality.winter.sta,twglaciers(gid).seasonality.winter.z,...
        'linewidth',lw1,'color',cols(length(twglaciers(gid).interannual)+1,:));
    plot(twglaciers(gid).seasonality.winter.st,twglaciers(gid).seasonality.winter.z,...
        'linewidth',lw2,'color',cols(length(twglaciers(gid).interannual)+1,:));
    plot([25,28],twglaciers(gid).seasonality.winter.nb*[1,1],'--',...
        'linewidth',lw2,'color',cols(length(twglaciers(gid).interannual)+1,:));

    xlabel('$\sigma_{\theta}$ (kg$\,$m$^{-3}$)','fontsize',fs,'interpreter','latex');
    ylabel('depth (m)','fontsize',fs,'interpreter','latex');
%     legend(q,'a','b','c','d','e','f','location','southwest','fontsize',fs,'interpreter','latex');
    xlim([25,28]); ylim([-900,0]);
    set(gca,'box','on','fontsize',fs,'ytick',[-800:200:0],...
        'yticklabel',{'800','600','400','200','0'});
%     title('Helheim: stratification','fontsize',fs-1,'interpreter','latex');
    
    a2 = axes('position',[lspace+1*(pw+hspace),bspace,pw,ph]); hold on;
    plot([25,28],twglaciers(gid).gldepth*[1,1],'linewidth',lw1,'color',0.5*[1,1,1]);
    z = twglaciers(gid).plume(1).sol.z;
    sta = rho(twglaciers(gid).plume(1).sol.Ta,twglaciers(gid).plume(1).sol.Sa,0)-1000;
    plot(sta,z,'k','linewidth',lw1);
    for i=1:7,
        z = twglaciers(gid).plume(i).sol.z;
        st = rho(twglaciers(gid).plume(i).sol.T,twglaciers(gid).plume(i).sol.S,0)-1000;
        plot(st,z,'linewidth',lw2,'color',cols(i,:));
        plot([25,28],twglaciers(gid).plume(i).znb*[1,1],'--','linewidth',lw2,'color',cols(i,:));
    end
    xlabel('$\sigma_{\theta}$ (kg$\,$m$^{-3}$)','fontsize',fs,'interpreter','latex');
    xlim([25,28]); ylim([-900,0]);
    set(gca,'box','on','fontsize',fs,'ytick',[-800:200:0],'yticklabel',{});
%     title('Helheim: discharge/hydrology','fontsize',fs-1,'interpreter','latex');
    
    gid = 5;
    a3 = axes('position',[lspace+2*(pw+hspace),bspace,pw,ph]); hold on;
    plot([25,28],twglaciers(gid).gldepth*[1,1],'linewidth',lw1,'color',0.5*[1,1,1]);
    for i=1:length(twglaciers(gid).interannual),
        plot(twglaciers(gid).interannual(i).sta,twglaciers(gid).interannual(i).z,'linewidth',lw1,'color',cols(i,:));
        plot(twglaciers(gid).interannual(i).st,twglaciers(gid).interannual(i).z,'linewidth',lw2,'color',cols(i,:));
        plot([25,28],twglaciers(gid).interannual(i).nb*[1,1],'--','linewidth',lw2,'color',cols(i,:));
    end
    plot(twglaciers(gid).seasonality.summer.sta,twglaciers(gid).seasonality.summer.z,...
        'linewidth',lw1,'color',cols(length(twglaciers(gid).interannual)+1,:));
    plot(twglaciers(gid).seasonality.summer.st,twglaciers(gid).seasonality.summer.z,...
        'linewidth',lw2,'color',cols(length(twglaciers(gid).interannual)+1,:));
    plot([25,28],twglaciers(gid).seasonality.summer.nb*[1,1],'--',...
        'linewidth',lw2,'color',cols(length(twglaciers(gid).interannual)+1,:));
    xlabel('$\sigma_{\theta}$ (kg$\,$m$^{-3}$)','fontsize',fs,'interpreter','latex');
    xlim([25,28]); ylim([-900,0]);
    set(gca,'box','on','fontsize',fs,'ytick',[-800:200:0],'yticklabel',{});
%     title('Kangilliup: stratification','fontsize',fs-1,'interpreter','latex');
        
    a4 = axes('position',[lspace+3*(pw+hspace),bspace,pw,ph]); hold on;
    plot([25,28],twglaciers(gid).gldepth*[1,1],'linewidth',lw1,'color',0.5*[1,1,1]);
    z = twglaciers(gid).plume(1).sol.z;
    sta = rho(twglaciers(gid).plume(1).sol.Ta,twglaciers(gid).plume(1).sol.Sa,0)-1000;
    plot(sta,z,'k','linewidth',lw1);
    for i=1:7,
        z = twglaciers(gid).plume(i).sol.z;
        st = rho(twglaciers(gid).plume(i).sol.T,twglaciers(gid).plume(i).sol.S,0)-1000;
        plot(st,z,'linewidth',lw2,'color',cols(i,:));
        plot([25,28],twglaciers(gid).plume(i).znb*[1,1],'--','linewidth',lw2,'color',cols(i,:));
    end
    xlabel('$\sigma_{\theta}$ (kg$\,$m$^{-3}$)','fontsize',fs,'interpreter','latex');
    xlim([25,28]); ylim([-900,0]);
    set(gca,'box','on','fontsize',fs,'ytick',[-800:200:0],'yticklabel',{});
%     title('Kangilliup: discharge/hydrology','fontsize',fs-1,'interpreter','latex');
    
    saveplot(17,8,300,'figS5.png');
    close all;
    
end

%% histograms of characteristic depths - fig. S2
if plotnum == 6,
    
% get data
for i=1:length(twglaciers),
    fjordvol(i) = twglaciers(i).fjord.vol;
    fjordmeandepth(i) = twglaciers(i).fjord.meandepth;
    fjordlength(i) = twglaciers(i).fjord.length;
    fjordwidth(i) = twglaciers(i).fjord.width;
    gldepth(i) = twglaciers(i).gldepth;
    effdepth(i) = twglaciers(i).fjord.effdepth;
    for j=1:5
        nbdepth(i,j) = twglaciers(i).plume(j).znb;
        mhdepth(i,j) = twglaciers(i).plume(j).sol.z(end);
    end
end

lspace = 0.06;
hspace = 0.05;
rspace = 0.01;
pw = (1-lspace-rspace-3*hspace)/4;
bspace = 0.14;
tspace = 0.05;
ph = 1-bspace-tspace;
xlims = [0,1050];
ylims = [0,60];
lw = 1.5;
cw = 0.5*[1,1,1];
tx = 950;
ty = 56;

figure();
a1 = axes('position',[lspace+0*(pw+hspace),bspace,pw,ph]); hold on;
histogram(abs(gldepth),[0:50:1050],'facecolor',cols(1,:));
plot(prctile(abs(gldepth),25)*[1,1],ylims,':','color',cw,'linewidth',1);
plot(prctile(abs(gldepth),50)*[1,1],ylims,'color',cw,'linewidth',1);
plot(prctile(abs(gldepth),75)*[1,1],ylims,':','color',cw,'linewidth',1);
set(gca,'fontsize',fs,'box','on');
xlabel('grounding line depth (m)','fontsize',fs);
ylabel('number of glaciers','fontsize',fs);
xlim(xlims); ylim(ylims);
text(tx,ty,'\textbf{a}','fontsize',fs+2);
set(gca,'xtick',[0:300:1000]);

a2 = axes('position',[lspace+1*(pw+hspace),bspace,pw,ph]); hold on;
histogram(abs(effdepth),[0:50:1050],'facecolor',cols(2,:));
plot(prctile(abs(effdepth),25)*[1,1],ylims,':','color',cw,'linewidth',1);
plot(prctile(abs(effdepth),50)*[1,1],ylims,'color',cw,'linewidth',1);
plot(prctile(abs(effdepth),75)*[1,1],ylims,':','color',cw,'linewidth',1);
set(gca,'fontsize',fs,'box','on','yticklabel',[]);
xlabel('fjord effective depth (m)','fontsize',fs);
xlim(xlims); ylim(ylims);
text(tx,ty,'\textbf{b}','fontsize',fs+2);
set(gca,'xtick',[0:300:1000]);

a3 = axes('position',[lspace+2*(pw+hspace),bspace,pw,ph]); hold on;
histogram(abs(squeeze(nbdepth(:,1))),[0:50:1050],'facecolor',cols(4,:));
plot(prctile(abs(squeeze(nbdepth(:,1))),25)*[1,1],ylims,':','color',cw,'linewidth',1);
plot(prctile(abs(squeeze(nbdepth(:,1))),50)*[1,1],ylims,'color',cw,'linewidth',1);
plot(prctile(abs(squeeze(nbdepth(:,1))),75)*[1,1],ylims,':','color',cw,'linewidth',1);
set(gca,'fontsize',fs,'box','on','yticklabel',[]);
xlabel('plume neutral buoyancy depth (m)','fontsize',fs);
xlim(xlims); ylim(ylims);
text(tx,ty,'\textbf{c}','fontsize',fs+2);
set(gca,'xtick',[0:300:1000]);

a4 = axes('position',[lspace+3*(pw+hspace),bspace,pw,ph]); hold on;
histogram(abs(squeeze(mhdepth(:,1))),[0:50:1050],'facecolor',cols(3,:));
plot(prctile(abs(squeeze(mhdepth(:,1))),25)*[1,1],ylims,':','color',cw,'linewidth',1);
plot(prctile(abs(squeeze(mhdepth(:,1))),50)*[1,1],ylims,'color',cw,'linewidth',1);
plot(prctile(abs(squeeze(mhdepth(:,1))),75)*[1,1],ylims,':','color',cw,'linewidth',1);
set(gca,'fontsize',fs,'box','on','yticklabel',[]);
xlabel('plume maximum height (m)','fontsize',fs);
xlim(xlims); ylim(ylims);
text(tx,ty,'\textbf{d}','fontsize',fs+2);
set(gca,'xtick',[0:300:1000]);

saveplot(17,6,300,'figS2.png');
close all;
        
end

%% relationships between neutral buoyancy and grounding line depth - fig. S3
if plotnum == 7,
    
% get data
for i=1:length(twglaciers),
    gldepth(i) = abs(twglaciers(i).gldepth);
    nbdepth(i) = abs(twglaciers(i).plume(1).znb);
end

max(gldepth)

lspace = 0.13;
rspace = 0.01;
pw = 1-lspace-rspace;
bspace = 0.14;
tspace = 0.04;
ph = 1-bspace-tspace;
xlims = [0,1050];
ylims = [0,400];
tx = 30;
ty = 370;
dy = 25;

p = polyfit(gldepth,nbdepth,1);
r = corrcoef(gldepth,nbdepth);
r2 = r(1,2)^2;

figure();
a1 = axes('position',[lspace,bspace,pw,ph]); hold on;
plot(gldepth,nbdepth,'k.','markersize',8);
plot([min(gldepth),max(gldepth)],p(1)*[min(gldepth),max(gldepth)]+p(2),'--','color',cols(7,:),'linewidth',1);
text(tx,ty,['r$^2$ = ',num2str(0.01*round(100*r2))],'fontsize',fs,'color',cols(7,:));
text(tx,ty-dy,['slope of fit = ',num2str(0.01*round(100*p(1)))],'fontsize',fs,'color',cols(7,:));
set(gca,'fontsize',fs,'box','on');
xlabel('grounding line depth (m)','fontsize',fs);
ylabel('plume neutral buoyancy depth (m)','fontsize',fs);
xlim(xlims); ylim(ylims);

saveplot(8,6,300,'figS3.png');
close all;
    
    
end
        
end