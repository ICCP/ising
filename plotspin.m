function plotspin(i,j,k)
    if k < 0
        xval = double(i) - 0.5;
        yval = double(j) - 0.5;
        text(xval,yval,'\downarrow','fontsize',15,'VerticalAlignment', ...
        'middle','HorizontalAlignment','center','Color','m')
    else
        xval = double(i) - 0.5;
        yval = double(j) - 0.5;
        text(xval,yval,'\uparrow','fontsize',15,'VerticalAlignment', ...
        'middle','HorizontalAlignment','center','Color','b')
    end
end
