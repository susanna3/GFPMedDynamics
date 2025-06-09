function results = process_medGFP_datasets(WTData, SogData, SogHetData, ZenData, binedges, Path, DataStructure)
% Process and plot MedGFP datasets for WT, sog, soghet, and zen
% Inputs:
%   WTData, SogData, SogHetData, ZenData - matrices where col 1 = position, cols 3:end = time series
%   Live or Fixed with use a different function
% Outputs:
%   results - struct with fields for each dataset: Med, Up, Low, Smooth, Rate

    % Dataset definitions
    datasetNames = {'wt','sog','soghet','zen'};
    datasetVars = {WTData, SogData, SogHetData, ZenData};
    datasetCols = {'Blues','PuRd','Oranges','Greens'};
    results = DataStructure;

    % Process each dataset
    for i = 1:numel(datasetNames)
        name = datasetNames{i};
        data = datasetVars{i};
        colorscheme = datasetCols{i};
        cmap = cbrewer2(colorscheme, numel(binedges)+10);
        stats = compute_binned_stats(data, binedges);
        [smooth, rate] = plot_smoothed_bins(stats.Med, binedges, cmap);
        stats.Smooth = max(smooth, 0);  % clamp to 0
        stats.Rate = rate;
        results.(name) = stats;
        saveas(gcf, strcat(Path,sprintf('Presentations/PlacedImages/MedPlot_%s.svg', name)));
    end
end

function stats = compute_binned_stats(data, binedges)
    posBins = discretize(data(:,1), binedges, 'IncludedEdge','left');
    allData = [posBins, data];
    nbins = numel(binedges) - 1;
    ntime = 60;  

    BinnedMed = zeros(nbins,ntime);
    BinnedUpMed = zeros(nbins,ntime);
    BinnedLowMed = zeros(nbins,ntime);

    for binno = 1:max(posBins)
        nucinbin = allData(allData(:,1)==binno,3:end);
        for t = 1:ntime
            nucatt = nucinbin(:,t);
            avgmed = mean(nucatt);
            devmed = std(nucatt) ./ sqrt(length(nucatt));
            BinnedMed(binno,t) = avgmed;
            BinnedUpMed(binno,t) = avgmed + devmed;
            BinnedLowMed(binno,t) = avgmed - devmed;
        end
    end

    stats.Med = BinnedMed;
    stats.Up = BinnedUpMed;
    stats.Low = BinnedLowMed;
end

function [BinnedSmooth, BinnedRate] = plot_smoothed_bins(BinnedMed, binedges, cmap)
    BinnedSmooth = [];
    BinnedRate = [];
    figure; hold on;

    for x = 1:(numel(binedges)-1)
        xdata = 1:size(BinnedMed,2);
        ydata = BinnedMed(x,:);
        ft = fit(xdata', ydata', 'SmoothingSpline', 'SmoothingParam', 0.01);
        fittedmed = feval(ft, xdata);
        rateofmed = differentiate(ft, xdata);
        BinnedSmooth = [BinnedSmooth, fittedmed];
        BinnedRate = [BinnedRate, rateofmed];
        c = cmap(numel(binedges)+10-x,:);
        scatter(xdata, ydata, 'MarkerEdgeColor', c, 'MarkerFaceColor', c);
        plot(xdata, fittedmed, 'Color', c);
    end
    xlim([20,60]);
    ylim([0,0.05]);
    hold off;

end
