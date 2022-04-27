% function to plot error ranges
function [] = errplot(gca,x,yvals,col,s,lw);

miny = min(yvals);
maxy = max(yvals);

plot(x*[1,1],[miny,maxy],'color',col,'linewidth',lw);
plot(x,miny,s,'markerfacecolor',col,'markeredgecolor',col);
plot(x,maxy,s,'markerfacecolor',col,'markeredgecolor',col);