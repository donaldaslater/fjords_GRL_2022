%% script to make plot for each fjord
% provided as part of Slater 2022 GRL
clear; close all;

load fjords.mat
load twglaciers.mat
[X,Y] = meshgrid(x,y);

%% plots of every system

cmap_bed = cptcmap('GMT_relief');
pad = 20e3;
lspace = 0.01;
rspace = 0.01;
tspace = 0.01;
bspace = 0.3;
pw = 1-lspace-rspace;
ph = 1-bspace-tspace;

for i=1:length(twglaciers),
% for i=1:1,
    
    % plot dimensions
    minx = min(X(twglaciers(i).fjord.inds));
    maxx = max(X(twglaciers(i).fjord.inds));
    miny = min(Y(twglaciers(i).fjord.inds));
    maxy = max(Y(twglaciers(i).fjord.inds));
    % check if cast is inside these bounds
    xp = twglaciers(i).profile.x;
    yp = twglaciers(i).profile.y;
    if xp>minx & xp<maxx & yp>miny & yp<maxy,
        disp('Cast is inside plot');
    else
        disp('Cast is outside plot: fixing...');
        if xp<minx,
            minx = xp;
            if yp<miny, miny = yp;
            elseif yp>maxy, maxy = yp;
            end
        elseif xp>maxx,
            maxx = xp;
            if yp<miny, miny = yp;
            elseif yp>maxy, maxy = yp;
            end
        elseif yp<miny,
            miny = yp;
        elseif yp>maxy,
            maxy = yp;
        end
    end
    dx = maxx-minx;
    dy = maxy-miny;
    [dmax,coord] = max([dx,dy]);
    if coord==1,
        xlims=[minx-pad,maxx+pad];
        ylims=0.5*(miny+maxy)+(dmax/2+pad)*[-1,1];
    else
        xlims=0.5*(minx+maxx)+(dmax/2+pad)*[-1,1];
        ylims=[miny-pad,maxy+pad];
    end
    fb2 = fb;
    fb2(twglaciers(i).fjord.inds)=2;
    
    % plot
    figure('Visible','off');
    a1 = axes('position',[lspace,bspace,pw,ph]); hold on;
    pcolor2(x,y,b); axis xy;
    colormap(a1,cmap_bed); caxis([-1000,1000]);
    xlim(xlims); ylim(ylims);
    set(a1,'box','on','xtick',[],'ytick',[]);
    a2 = axes('position',get(a1,'position')); hold on;
    pcolor2(x,y,m); axis xy;
    colormap(a2,0.75*[1,1,1]);
    q(1)=plot(xb,yb,'k--','linewidth',1);
    contour(x,y,fb2,[2,2],'r','linewidth',1);
    q(3)=plot(NaN,NaN,'r','linewidth',1);
    q(5)=plot([twglaciers.x],[twglaciers.y],'kp','markerfacecolor','y','markersize',8);
    q(4)=plot(twglaciers(i).x,twglaciers(i).y,'kp','markerfacecolor','r','markersize',13);
    q(2)=plot(X(twglaciers(i).fjord.b_inds),Y(twglaciers(i).fjord.b_inds),'k^',...
        'markerfacecolor','r','markersize',5);
    q(6)=plot(twglaciers(i).profile.x,twglaciers(i).profile.y,'ko',...
        'markerfacecolor','r','markersize',7);
    set(a2,'box','off','xtick',[],'ytick',[],'visible','off','color','none');
    xlim(xlims); ylim(ylims);
    l=legend(q,{'fjords-shelf boundary','fjord mouth','fjord extent','glacier','other glaciers','CTD cast'});
    lp = l.Position;
    l.Position = [0.75,bspace-0.13,0,0];
    % stats
    dx = 0.02;
    dy = 0.02;
    d1 = 0.01;
    
    annotation('textbox',[lspace+dx,bspace-0*dy-1*d1,1,0],'string',['Glacier: ',twglaciers(i).name,' (',num2str(0.01*round(100*twglaciers(i).lat)),char(176),'N, ',num2str(0.01*round(100*abs(twglaciers(i).lon))),char(176),'W)'],'fontsize',12,'interpreter','tex','edgecolor','none');
    
    annotation('textbox',[lspace+dx,bspace-1*dy-3*d1,1,0],'string',['Area = ',num2str(round(twglaciers(i).fjord.area)),' km$^2$'],'fontsize',12,'interpreter','latex','edgecolor','none');
    annotation('textbox',[lspace+dx,bspace-2*dy-3*d1,1,0],'string',['Volume = ',num2str(round(twglaciers(i).fjord.vol)),' km$^3$'],'fontsize',12,'interpreter','latex','edgecolor','none');
    annotation('textbox',[lspace+dx,bspace-3*dy-3*d1,1,0],'string',['Length = ',num2str(round(twglaciers(i).fjord.length)),' km'],'fontsize',12,'interpreter','latex','edgecolor','none');
    
    annotation('textbox',[lspace+dx,bspace-4*dy-4*d1,1,0],'string',['Mean depth = ',num2str(round(twglaciers(i).fjord.meandepth)),' m'],'fontsize',12,'interpreter','latex','edgecolor','none');
    annotation('textbox',[lspace+dx,bspace-5*dy-4*d1,1,0],'string',['Max depth = ',num2str(round(twglaciers(i).fjord.maxdepth)),' m'],'fontsize',12,'interpreter','latex','edgecolor','none');
    annotation('textbox',[lspace+dx,bspace-6*dy-4*d1,1,0],'string',['GL depth = ',num2str(round(abs(twglaciers(i).gldepth))),' m'],'fontsize',12,'interpreter','latex','edgecolor','none');
    annotation('textbox',[lspace+dx,bspace-7*dy-4*d1,1,0],'string',['Effective depth = ',num2str(round(abs(twglaciers(i).fjord.effdepth))),' m'],'fontsize',12,'interpreter','latex','edgecolor','none');
    annotation('textbox',[lspace+dx,bspace-8*dy-4*d1,1,0],'string',['Plume NB depth = ',num2str(round(abs(twglaciers(i).plume(1).znb))),' m'],'fontsize',12,'interpreter','latex','edgecolor','none');
    
    annotation('textbox',[lspace+dx,bspace-9*dy-5*d1,1,0],'string',['Mean SGD flux = ',num2str(round(twglaciers(i).plume(1).Qgl)),' m$^3$/s = ',num2str(0.01*round(100*twglaciers(i).plume(1).Qgl/1000)),' mSv'],'fontsize',12,'interpreter','latex','edgecolor','none');
    annotation('textbox',[lspace+dx,bspace-10*dy-5*d1,1,0],'string',['Upwelled flux = ',num2str(round(twglaciers(i).plume(1).Qnb)),' m$^3$/s = ',num2str(0.01*round(100*twglaciers(i).plume(1).Qnb/1000)),' mSv'],'fontsize',12,'interpreter','latex','edgecolor','none');
   
    if i<10, savenum=['00',num2str(i)];
    elseif i<100, savenum=['0',num2str(i)];
    else savenum=num2str(i);
    end
    fw = 14;
    fh = fw*pw/ph;
    saveplot(fw,fh,300,['glacier_fjord_',savenum,'.png']);
    close all;
    
end