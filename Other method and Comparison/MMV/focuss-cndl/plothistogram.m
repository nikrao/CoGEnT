function plothistogram(data, rows, cols, numbins)
% plothistogram   - Plots the histogram of each row of the data matrix into
%                   a single figure
%
% plothistogram(data, rows, cols)
%
% data      - Data matrix, the histogram of each row is plotted
% rows      - Rows to plot
% cols      - Columns to plot
% numbins   - Number of bins in each histogram
%
% JFM    3/6/2001
% Rev:   3/10/2001

index = 1;

datarows = size(data,1);
clf;

for i = 1:rows
    for j = 1:cols
        if(index <= datarows)
            subplot(rows, cols, index);
            hist(data(index,:), numbins);
            ax = axis;
            ax(1) = -10; ax(2) = 10;
            axis(ax);
            setfonts;        
            index = index + 1;
        end
    end
end

return;

        
        