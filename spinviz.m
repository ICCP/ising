function tridip = spinviz
% This program visualizes the spin for Ising model for
% 2D lattice system at the end of the iteration
% BEFORE execution check ising.f for parameter setup in WRITE(31,*) file
clear all;
clc;
fileID = fopen('Spin_states.out');
%
header = textscan(fileID,'%s %d %s %d',1);
NROWS = header{2}; NCOLS = header{4};
lattice_sites = NROWS * NCOLS;
datain = textscan(fileID,'%d',lattice_sites);
spins = datain{1};
fclose(fileID);
%
% Plot grid/box for spins to be displayed
hold on
plot([0 0],[0 NROWS],'k','LineWidth',2);
plot([0 NCOLS],[0 0],'k','LineWidth',2);
%
for i = 1:NROWS
    plot([i i],[0 NROWS],'k','LineWidth',2);
end
%
for j = 1:NCOLS
    plot([0 NROWS],[j j],'k','LineWidth',2);
end
% Plot spins inside box
for I = 1:NROWS
    for J = 1:NCOLS
        spin = spins((I-1)* NROWS + J);
        plotspin(I,J,spin)
    end
end
 axis equal
 axis off
 set(gcf,'Color',[1,1,1]);
hold off
%
set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'spin_visualization', 'pdf') %Save figure
end
%
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