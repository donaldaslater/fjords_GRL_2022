% function to make custom, usually tighter legend
function [] = customlegend(ax,p,x0,y0,dx,dy,w,l,fs);

axes(ax); hold on;

for i=1:length(p),
    if strcmp(p(i).LineStyle,'none'),
        plot(ax,x0,y0+(i-1)*dy,'color',p(i).Color,'LineStyle',p(i).LineStyle,...
                'LineWidth',p(i).LineWidth,'Marker',p(i).Marker,...
                'MarkerSize',p(i).MarkerSize,'MarkerFaceColor',p(i).MarkerFaceColor);
    else
        plot(ax,x0+(w/2)*[-1,1],y0+(i-1)*dy*[1,1],'color',p(i).Color,'LineStyle',p(i).LineStyle,...
                'LineWidth',p(i).LineWidth,'Marker',p(i).Marker,...
                'MarkerSize',p(i).MarkerSize,'MarkerFaceColor',p(i).MarkerFaceColor);
    end
    text(x0+dx,y0+(i-1)*dy,l{i},'fontsize',fs,'verticalalignment','middle',...
        'horizontalalignment','left');
end