function makeScatterPlot(directory)
    addpath(['.',filesep,'MATLAB_scripts']);
    if nargin < 1
        directory = uigetdir('.','Select the directory to plot from');
    end
    % Try to load an existing aggregate.
    try
        data = load([directory,filesep,'aggregate.mat']);
        data = data.aggregate;
        disp('Existing aggregate loaded.')
    catch
        data = aggregateData(directory);
    end

    % Plot aggregated data. The field names should EXACTLY match the variable names in 'processFile.m'.
    figure
    x = data.periodAutoCorr;
    y = data.periodFourier;
    scatter(x, y, 'Marker', 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'black', 'SizeData', 100)
    xlabel('Period from autocorrelation (s)')
    ylabel('Period from Fourier analysis (s)')
    title('This is the title of the scatter plot')
    axis equal
    grid on
    set(gca,'FontSize',16)
    
    % Save the axis limits for later use.
    xlims = xlim;
    ylims = ylim;

    % Plot the line of best fit.
    hold on
    coefficients = polyfit(x, y, 1);
    plot(xlim, polyval(coefficients, xlim), 'LineWidth', 1, 'Color', 'black')
    % Restore the axis limits.
    xlim(xlims)
    ylim(ylims)

    exportgraphics(gcf,'lastScatter.png')
end