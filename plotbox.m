function [] = plotbox(ax,x,y,w,c,lw);

% inputs
% ax - axes
% x - xvalue to plot
% y - vector of values
% w - width of box
% c - color
% lw - linewidth

axes(ax); hold on;
plot(x+0.5*[-w,w],prctile(y,50)*[1,1],'color',c,'linewidth',lw);
plot(x+0.5*[-w,w],prctile(y,25)*[1,1],'color',c,'linewidth',lw);
plot(x+0.5*[-w,w],prctile(y,75)*[1,1],'color',c,'linewidth',lw);
plot((x-0.5*w)*[1,1],[prctile(y,25),prctile(y,75)],'color',c,'linewidth',lw);
plot((x+0.5*w)*[1,1],[prctile(y,25),prctile(y,75)],'color',c,'linewidth',lw);
plot(x+0.25*[-w,w],min(y)*[1,1],'color',c,'linewidth',lw);
plot(x+0.25*[-w,w],max(y)*[1,1],'color',c,'linewidth',lw);
plot(x*[1,1],[min(y),prctile(y,25)],'--','color',c,'linewidth',lw);
plot(x*[1,1],[max(y),prctile(y,75)],'--','color',c,'linewidth',lw);


