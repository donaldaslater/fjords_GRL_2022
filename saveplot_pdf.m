% function to save plot and avoid usual block of code

function [] = saveplot_pdf(fw,fh,res,filename);

    set(gcf,'PaperUnits','centimeters');
    set(gcf,'inverthardcopy','off');
    set(gcf,'PaperPosition',[0 0 fw fh]);
    set(gcf,'PaperSize',[fw fh]);
    set(gcf,'units','centimeters','position',[0 0 fw fh]);
    set(gcf,'color','w');
    res_string = ['-r',num2str(res)];
    print('-dpdf',res_string,filename);

end