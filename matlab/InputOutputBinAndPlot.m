%% Set Image Colors

Path = 'C:/Users/Susanna Brantley/Brantley Lab Dropbox/Susanna Brantley/Duke/DiTaliaLab/';
%Path = 'D:/Brantley Lab Dropbox/Susanna Brantley/Duke/DiTaliaLab/';
%Path = '/Users/susannabrantley/Dropbox/Duke/DiTaliaLab/';

DataPath = 'ImageAnalysis/Datasets/';

mymagenta = '#D81B60';
myblue = '#1E88E5';
mygreen = '#004D40';
myyellow = '#FFC107';

FigurePath = strcat(Path,'Presentations/PlacedImages/');

set(0,'DefaultAxesFontName','Arial');

% load heatmap settings
mymagentamap = load('mymagentamap.mat').mymagentamap;
mybluemap = load('mybluemap.mat').mybluemap;
mygreenmap = load('mygreenmap.mat').mygreenmap;
myyellowmap = load('myyellowmap.mat').myyellowmap;

s.matlab.fonts.editor.code.Name.TemporaryValue = 'Arial';

%% Load Data

genename = "Ush";

WTData = load("AllDataWTLive.mat", 'AllTrackedNuclei');
WTData = WTData.AllTrackedNuclei;

SogHetData = load("AllDataSogHetLive.mat", 'AllTrackedNucleiSog');
SogHetData = SogHetData.AllTrackedNucleiSog;

SogData = load("AllDataSogLive.mat", 'AllTrackedNucleiSog');
SogData = SogData.AllTrackedNucleiSog;

ZenData = load("AllDataZenLive.mat","ZenAllTrackedNuclei");
ZenData = ZenData.ZenAllTrackedNuclei;

WTDataFixed = load(strcat("AllDataWT.mat"), 'AllWTNucleiData');
WTTiming = load(strcat("AllTimingWT.mat"), 'AllWTNucleiTiming');
WTDataFixed =WTDataFixed.AllWTNucleiData;
WTTiming = WTTiming.AllWTNucleiTiming;

SogDataFixed = load(strcat("AllDataSog.mat"), 'AllMutantNucleiData');
SogTiming = load(strcat("AllTimingSog.mat"), 'AllMutantNucleiTiming');
SogDataFixed = SogDataFixed.AllMutantNucleiData;
SogTiming = SogTiming.AllMutantNucleiTiming;

SogHetDataFixed = load(strcat("AllDataSogHet.mat"), 'AllMutantNucleiData');
SogHetTiming = load(strcat("AllTimingSogHet.mat"), 'AllMutantNucleiTiming');
SogHetDataFixed = SogHetDataFixed.AllMutantNucleiData;
SogHetTiming = SogHetTiming.AllMutantNucleiTiming;

ZenDataFixed = load(strcat("AllDataZen.mat"), 'AllMutantNucleiData');
ZenTiming = load(strcat("AllTimingZen.mat"), 'AllMutantNucleiTiming');
ZenDataFixed = ZenDataFixed.AllMutantNucleiData;
ZenTiming = ZenTiming.AllMutantNucleiTiming;

AllWTDataUsh = load(strcat("AllDataWT",genename,".mat"), 'AllWTNucleiData');
AllWTTimingUsh = load(strcat("AllTimingWT",genename,".mat"), 'AllWTNucleiTiming');
AllWTDataUsh = AllWTDataUsh.AllWTNucleiData;
AllWTTimingUsh = AllWTTimingUsh.AllWTNucleiTiming;

AllSogDataUsh = load(strcat("AllDataSog",genename,".mat"), 'AllMutantNucleiData');
AllSogTimingUsh = load(strcat("AllTimingSog",genename,".mat"), 'AllMutantNucleiTiming');
AllSogDataUsh = AllSogDataUsh.AllMutantNucleiData;
AllSogTimingUsh = AllSogTimingUsh.AllMutantNucleiTiming;

AllSogHetDataUsh = load(strcat("AllDataSogHet",genename,".mat"), 'AllMutantNucleiData');
AllSogHetTimingUsh = load(strcat("AllTimingSogHet",genename,".mat"), 'AllMutantNucleiTiming');
AllSogHetDataUsh = AllSogHetDataUsh.AllMutantNucleiData;
AllSogHetDataUsh = AllSogHetDataUsh';
AllSogHetTimingUsh = AllSogHetTimingUsh.AllMutantNucleiTiming;

AllZenDataUsh = load(strcat("AllDataZen",genename,".mat"), 'AllMutantNucleiData');
AllZenTimingUsh = load(strcat("AllTimingZen",genename,".mat"), 'AllMutantNucleiTiming');
AllZenDataUsh = AllZenDataUsh.AllMutantNucleiData;
AllZenTimingUsh = AllZenTimingUsh.AllMutantNucleiTiming;


%% Bin Live Med, plot binned data
binedges = 0:5:75;
timebinedges = [0,25:5:60];
results = process_medGFP_datasets(WTData, SogData, SogHetData, ZenData, binedges, Path, struct());
results = process_psmad_datasets(WTDataFixed,WTTiming,SogDataFixed,SogTiming,SogHetDataFixed,SogHetTiming,ZenDataFixed,ZenTiming, ...
    binedges, timebinedges, Path, results);

%% Plot rate of Med 

genotypeNames = fieldnames(results);
%binedges = 0:5:75;  % Make sure this matches what was used in processing
datasetCols = {'Blues','PuRd','Oranges','Greens'};

for i = 1:numel(genotypeNames)
    cmap = cbrewer2(datasetCols{i}, numel(binedges) + 10);
    genotype = genotypeNames{i};
    medRate = results.(genotype).Rate;

    figure;
    hold on;
    fontsize(14, "points")  % If this doesn't work, use set(gca, 'FontSize', 14)
    
    for x = 1:(numel(binedges)-1)  % One fewer bin
        xdata = 1:size(medRate,1);  % Timepoints
        ydata = medRate(:,x);      % Rate of change for bin x
        c = cmap(numel(binedges)+10-x,:);
        plot(xdata', ydata, 'Color', c, 'LineWidth', 2, 'LineStyle', '--');
    end

    xlabel('Time (Minutes)');
    ylabel('dMed/dt');
    ylim([0, 0.003]);
    xlim([20, 60]);
    hold off;

    saveas(gcf, strcat(Path,sprintf('Presentations/PlacedImages/RatePlot_%s.svg', genotype)));
end

%% Calculate integrated Med

genotypeNames = fieldnames(results);
%binedges = 0:5:75;
datasetCols = {'Blues','PuRd','Oranges','Greens'};

inttime = 60;  % integration window size
offtime = 0;

for i = 1:numel(genotypeNames)
    colorscheme = datasetCols{i};
    cmap = cbrewer2(colorscheme, numel(binedges) + 10);
    genotype = genotypeNames{i};
    BinnedMed = results.(genotype).Smooth;  % use smoothed data
    
    MeanCumInt = [];
    figure;
    hold on;
    set(gca, 'FontSize', 14);

    for x = 1:(numel(binedges)-1)
        avgdata = BinnedMed(:,x)';  % time vector for this bin
        int = [zeros(offtime,1); cumtrapz(avgdata(offtime+1:inttime))];
        for t = (inttime+1):60
            intwindow = t-inttime;
            timeint = max(cumtrapz(avgdata(intwindow:t)));
            int = [int, timeint];
        end
        c = cmap(numel(binedges)+10-x,:);
        scatter(1:60, int, 'MarkerEdgeColor', c,'MarkerFaceColor',c);
        MeanCumInt = [MeanCumInt, int'];
    end

    xlabel('Time (Minutes)');
    ylabel('Cumulative Intensity');
    ylim([0, 1]);
    xlim([1, 60]);
    hold off;

    saveas(gcf, strcat(Path,sprintf('Presentations/PlacedImages/IntPlot_%s.svg', genotype)));

    results.(genotype).CumInt = MeanCumInt;
end

%% Bin nuclei from fixed data and get average pSMAD

genotypeNames = fieldnames(results);
RNADataNames = {'WT','Sog','SogHet','Zen'};

% Define coarse bin edges
spatial_edges = binedges;
nSpatial = numel(spatial_edges) - 1;

for g = 1:numel(genotypeNames)
    genotype = genotypeNames{g};
    dataname = RNADataNames{g};
    colorscheme = datasetCols{g};
    cmap = cbrewer2(colorscheme, numel(binedges)+10);
    fprintf('Processing RNA for %s...\n', genotype);

    % Construct data variable names
    dataVar = ['All', dataname, 'DataUsh'];
    RNAdata = evalin('base', dataVar);

    numEmbryos = size(RNAdata, 1);

    % Initialize storage across embryos
    FractionOnRaw = NaN(numEmbryos, nSpatial);
    MeanPSMADRaw = NaN(numEmbryos, nSpatial);

    % Loop through embryos
    for embryo = 1:numEmbryos
        embryodata = RNAdata{embryo};

        % Filter AP range if needed
        %embryodata = embryodata(embryodata(:,5) > 75 & embryodata(:,5) < 150, :);

        position = embryodata(:,1);       % spatial x
        psmad = embryodata(:,2);          % pSMAD level
        % figure
        % hold on
        % scatter(position,psmad)
        % ylim([0,0.5])
        % hold off
        expression = embryodata(:,3) > 0; % ON/OFF call

        % Bin in space only
        [~, ~, spaceBin] = histcounts(position, spatial_edges);

        for s = 1:nSpatial
            inBin = (spaceBin == s);
            if any(inBin)
                FractionOnRaw(embryo, s) = sum(expression(inBin)) / sum(inBin);
                MeanPSMADRaw(embryo, s) = mean(psmad(inBin));
            end
        end
    end

    % Store individual and averaged results
    results.(genotype).FractionOnRaw_all = FractionOnRaw;   % embryo × bin
    results.(genotype).MeanPSMADRaw_all = MeanPSMADRaw;
end

%% Bin gene expression but first smooth across distance for each embryo

genotypeNames = fieldnames(results);
%binedges = 0:5:75;
binCenters = (binedges(1:end-1) + binedges(2:end))/2;
datasetCols = {'Blues','PuRd','Oranges','Greens'};
RNADataNames = {'WT','Sog','SogHet','Zen'};

for g = 1:numel(genotypeNames)
    genotype = genotypeNames{g};
    dataname = RNADataNames{g};
    colorscheme = datasetCols{g};
    cmap = cbrewer2(colorscheme, numel(binedges)+10);
    fprintf('Processing RNA for %s...\n', genotype);

    % Construct data variable names
    dataVar = ['All', dataname, 'DataUsh'];
    timeVar = ['All', dataname, 'TimingUsh'];

    % Load RNA data and timings
    RNAdata = evalin('base', dataVar);
    RNAtimes = evalin('base', timeVar);

    numEmbryos = size(RNAdata, 1);
    fittedResults = cell(numEmbryos, 1);

    % First: Fit fraction ON at each bin for each embryo
    figure
    hold on
    cmap = cbrewer2('Blues',numEmbryos);
    for embryo = 1:numEmbryos
        embryodata = RNAdata{embryo};

        % Filter AP range (optional)
        %embryodata = embryodata(embryodata(:,5) > 75 & embryodata(:,5) < 150, :);
        position = embryodata(:,1);
        expression = embryodata(:,3) > 0;

        [~, ~, binIdx] = histcounts(position, binedges);
        fractionExpressing = zeros(1, numel(binedges)-1);

        for i = 1:(numel(binedges)-1)
            nucleiInBin = (binIdx == i);
            nInBin = sum(nucleiInBin);
            if nInBin > 0
                fractionExpressing(i) = sum(expression(nucleiInBin)) / nInBin;
            else
                fractionExpressing(i) = NaN;
            end
        end

        % Smooth along position
        validIdx = ~isnan(fractionExpressing);
        fit_x = binCenters(validIdx);
        fit_y = fractionExpressing(validIdx);

        if numel(fit_x) >= 2
            [fitresult, ~] = fit(fit_x.', fit_y.', 'smoothingspline', 'SmoothingParam', 0.005);
            eval_y = feval(fitresult, binCenters);
            eval_y = min(max(eval_y, 0), 1);  % Clip between [0,1]
        else
            eval_y = NaN(size(binCenters));
        end

        fittedResults{embryo} = [binCenters(:), eval_y(:)];
        scatter(fit_x,fit_y)
        plot(fit_x,feval(fitresult,fit_x),'Color',cmap(embryo,:))
    end
    hold off

    % Second: For each bin, fit logistic curve over time
    FractionOn = [];
    for pos_idx = 1:length(binCenters)
        position = binCenters(pos_idx);
        collectedTimes = RNAtimes(:,1);  % Embryo times
        collectedFracs = cellfun(@(data) data(pos_idx,2), fittedResults);

        % Optionally filter out early times
        keepIdx = collectedTimes > 25;
        collectedTimes = collectedTimes(keepIdx);
        collectedFracs = collectedFracs(keepIdx);

        if numel(collectedTimes) >= 3  % Enough points to fit
            try
                [logitCoef, ~] = glmfit(collectedTimes, collectedFracs, 'binomial', 'logit');
                logitFit = glmval(logitCoef, 1:60, 'logit');
            catch
                logitFit = nan(1,60);
            end
        else
            logitFit = nan(1,60);
        end

        FractionOn(:, pos_idx) = logitFit(:);
    end

    % Store final results in the struct
    results.(genotype).RNAfraction = FractionOn;
end

%% Plot probability v. med levels

useCumInt = true;  % Set true to use CumInt instead of SmoothMed
lagtime = 0;
%binedges = 0:5:75;
genotypeNames = fieldnames(results);
datasetCols = {'Blues','PuRd','Oranges','Greens'};
datasetCols2 = {myblue,mymagenta,myyellow,mygreen};

poscutoff = 9;

figure
hold on
legendEntries = {};

for g = 1:numel(genotypeNames)
    colorscheme = datasetCols{g};
    genotype = genotypeNames{g};
    cmap = cbrewer2(colorscheme, numel(binedges)-1);
    data = results.(genotype);
    
    % Select X-axis variable
    if useCumInt
        XAxisMat = data.CumInt;
        xlabelstr = 'Integrated Med';
    else
        XAxisMat = data.Smooth;
        xlabelstr = 'Instantaneous Med';
    end
    
    % Matrices: 60 x Nbins
    Y = data.RNAfraction;
    X = [zeros(lagtime, size(XAxisMat,2)); XAxisMat(1:(end-lagtime), :)];

    % Flatten for fit
    xdata = X(:);
    ydata = Y(:);

    % Remove NaNs or Infs
    valid = isfinite(xdata) & isfinite(ydata);
    xdata = xdata(valid);
    ydata = ydata(valid);
    
    % Fit Hill function
    startPoints = [1, 1];
    s = fitoptions('Method','NonlinearLeastSquares','Robust','LAR',...
        'Lower',[0 0],'Upper',[Inf Inf],'Startpoint',startPoints);
    HillEqn = fittype('x.^a1 ./ (x.^a1 + a2.^a1)', 'options', s);

    [fitobj, gof] = fit(xdata, ydata, HillEqn);

    % Plot scatter points
    thisColor = datasetCols2{g};
    scatter(xdata, ydata, 10, 'MarkerEdgeColor', thisColor, 'MarkerFaceColor', thisColor, 'MarkerFaceAlpha', 0.3)

    % Plot fitted curve
    xfit = linspace(min(xdata), max(xdata), 1000);
    yfit = feval(fitobj, xfit);
    plot(xfit, yfit, 'Color', thisColor, 'LineWidth', 2)

    % Add text with error in corner
    xpos = 0.6; ypos = 0.2 - 0.05*g;
    text(xpos, ypos, sprintf('%s: RMSE=%.5f', genotype, gof.rmse), ...
        'Color', thisColor, 'FontSize', 10, 'Units', 'normalized')
    %legendEntries{end+1} = genotype;
end

xlabel(xlabelstr)
ylabel('Fraction Nuclei On')
ylim([0,1])
xlim([0, max(xfit)])
fontsize(14,"points")
hold off

if useCumInt
    filelabelstr = 'Int';
else
    filelabelstr = 'Inst';
end

saveas(gcf, strcat(Path,'Presentations/PlacedImages/HillPlots',filelabelstr,genename,'.svg'));

useCumInt = true;  % Set true to use CumInt instead of SmoothMed
lagtime = 0;
%binedges = 0:5:75;
nbins = length(binedges) - 1;
genotypeNames = fieldnames(results);

nGenos = numel(genotypeNames);
nCols = ceil(sqrt(nGenos));
nRows = ceil(nGenos / nCols);

figure
for g = 1:nGenos
    colorscheme = datasetCols{g};
    cmap = cbrewer2(colorscheme, nbins + 10);  % For spatial bins
    genotype = genotypeNames{g};
    data = results.(genotype);
    
    % Select X-axis matrix
    if useCumInt
        XAxisMat = data.CumInt;
        xlabelstr = 'Integrated Med';
    else
        XAxisMat = data.Smooth;
        xlabelstr = 'Instantaneous Med';
    end
    
    % Y = Fraction ON
    Y = data.RNAfraction;  % 60 x Nbins
    X = [zeros(lagtime, nbins); XAxisMat(1:(end-lagtime), :)];  % Shift if lagtime

    % Subplot setup
    subplot(nRows, nCols, g)
    hold on

    xdata_all = [];
    ydata_all = [];

    % Plot by spatial bin
    for b = 1:nbins
        xdata = X(:, b);
        ydata = Y(:, b);

        valid = isfinite(xdata) & isfinite(ydata);
        xdata = xdata(valid);
        ydata = ydata(valid);

        scatter(xdata, ydata, 12, 'MarkerEdgeColor', cmap(end-b,:), ...
            'MarkerFaceColor', cmap(end-b,:), 'MarkerFaceAlpha', 0.5)

        xdata_all = [xdata_all; xdata];
        ydata_all = [ydata_all; ydata];
    end

    % Fit Hill function to all pooled data
    valid = isfinite(xdata_all) & isfinite(ydata_all);
    xdata_fit = xdata_all(valid);
    ydata_fit = ydata_all(valid);

    startPoints = [1 1];
    s = fitoptions('Method','NonlinearLeastSquares','Robust','LAR', ...
        'Lower',[0 0],'Upper',[Inf Inf],'Startpoint',startPoints);
    HillEqn = fittype('x.^a1 ./ (x.^a1 + a2.^a1)', 'options', s);

    [fitobj, gof] = fit(xdata_fit, ydata_fit, HillEqn);

    % Plot fit
    xfit = linspace(min(xdata_fit), max(xdata_fit), 500);
    yfit = feval(fitobj, xfit);
    plot(xfit, yfit, 'k-', 'LineWidth', 2)

    % Annotate subplot
    text(0.6, 0.2, sprintf('RMSE=%.5f', gof.rmse), ...
        'Units', 'normalized', 'FontSize', 10)
    title(genotype)
    xlabel(xlabelstr)
    ylabel('Fraction Nuclei On')
    ylim([0 1])
    xlim([0, max(xfit)])
    ax = gca;                   % Get current axes
    colormap(ax, cmap);  % Use a specific colormap per subplot
    %colorbar(ax);  
    hold off
end

if useCumInt
    filelabelstr = 'Int';
else
    filelabelstr = 'Inst';
end

saveas(gcf, strcat(Path,'Presentations/PlacedImages/HillPlots',filelabelstr,genename,'2.svg'));
%%

useCumInt = false;  % Set true to use CumInt instead of SmoothMed
lagtime = 0;
%binedges = 0:5:75;
genotypeNames = fieldnames(results);
datasetCols = {'Blues','PuRd','Oranges','Greens'};
datasetCols2 = {myblue,mymagenta,myyellow,mygreen};

figure
hold on
legendEntries = {};

for g = 1:numel(genotypeNames)
    colorscheme = datasetCols{g};
    genotype = genotypeNames{g};
    cmap = cbrewer2(colorscheme, numel(binedges)-1);
    data = results.(genotype);
    
    % Select X-axis variable
    if useCumInt
        XAxisMat = data.CumInt;
        xlabelstr = 'Integrated Med';
    else
        XAxisMat = data.Smooth;
        xlabelstr = 'Instantaneous Med';
    end
    
    % Matrices: 60 x Nbins
    Y = data.RNAfraction;
    X = [zeros(lagtime, size(XAxisMat,2)); XAxisMat(1:(end-lagtime), :)];

    % Flatten for fit
    xdata = X(:);
    ydata = Y(:);

    % Remove NaNs or Infs
    valid = isfinite(xdata) & isfinite(ydata);
    xdata = xdata(valid);
    ydata = ydata(valid);
    
    % Fit Hill function
    startPoints = [1, 1];
    s = fitoptions('Method','NonlinearLeastSquares','Robust','LAR',...
        'Lower',[0 0],'Upper',[Inf Inf],'Startpoint',startPoints);
    HillEqn = fittype('x.^a1 ./ (x.^a1 + a2.^a1)', 'options', s);

    [fitobj, gof] = fit(xdata, ydata, HillEqn);

    % Plot scatter points
    thisColor = datasetCols2{g};
    scatter(xdata, ydata, 10, 'MarkerEdgeColor', thisColor, 'MarkerFaceColor', thisColor, 'MarkerFaceAlpha', 0.3)

    % Plot fitted curve
    xfit = linspace(min(xdata), max(xdata), 1000);
    yfit = feval(fitobj, xfit);
    plot(xfit, yfit, 'Color', thisColor, 'LineWidth', 2)

    % Add text with error in corner
    xpos = 0.6; ypos = 0.2 - 0.05*g;
    text(xpos, ypos, sprintf('%s: RMSE=%.5f', genotype, gof.rmse), ...
        'Color', thisColor, 'FontSize', 10, 'Units', 'normalized')
    legend off
    %legendEntries{end+1} = genotype;
end

xlabel(xlabelstr)
ylabel('Fraction Nuclei On')
ylim([0,1])
xlim([0, max(xfit)])
fontsize(14,"points")
hold off

if useCumInt
    filelabelstr = 'Int';
else
    filelabelstr = 'Inst';
end

saveas(gcf, strcat(Path,'Presentations/PlacedImages/HillPlots',filelabelstr,genename,'.svg'));

useCumInt = false;  % Set true to use CumInt instead of SmoothMed
lagtime = 0;
%binedges = 0:5:75;
nbins = length(binedges) - 1;
genotypeNames = fieldnames(results);

nGenos = numel(genotypeNames);
nCols = ceil(sqrt(nGenos));
nRows = ceil(nGenos / nCols);

figure
for g = 1:nGenos
    colorscheme = datasetCols{g};
    cmap = cbrewer2(colorscheme, nbins + 10);  % For spatial bins
    genotype = genotypeNames{g};
    data = results.(genotype);
    
    % Select X-axis matrix
    if useCumInt
        XAxisMat = data.CumInt;
        xlabelstr = 'Integrated Med';
        filelabelstr = 'Int';
    else
        XAxisMat = data.Smooth;
        xlabelstr = 'Instantaneous Med';
        filelabelstr = 'Inst';
    end
    
    % Y = Fraction ON
    Y = data.RNAfraction;  % 60 x Nbins
    X = [zeros(lagtime, nbins); XAxisMat(1:(end-lagtime), :)];  % Shift if lagtime

    % Subplot setup
    subplot(nRows, nCols, g)
    hold on

    xdata_all = [];
    ydata_all = [];

    % Plot by spatial bin
    for b = 1:nbins
        xdata = X(:, b);
        ydata = Y(:, b);

        valid = isfinite(xdata) & isfinite(ydata);
        xdata = xdata(valid);
        ydata = ydata(valid);

        scatter(xdata, ydata, 12, 'MarkerEdgeColor', cmap(end-b,:), ...
            'MarkerFaceColor', cmap(end-b,:), 'MarkerFaceAlpha', 0.5)

        xdata_all = [xdata_all; xdata];
        ydata_all = [ydata_all; ydata];
    end

    % Fit Hill function to all pooled data
    valid = isfinite(xdata_all) & isfinite(ydata_all);
    xdata_fit = xdata_all(valid);
    ydata_fit = ydata_all(valid);

    startPoints = [1 1];
    s = fitoptions('Method','NonlinearLeastSquares','Robust','LAR', ...
        'Lower',[0 0],'Upper',[Inf Inf],'Startpoint',startPoints);
    HillEqn = fittype('x.^a1 ./ (x.^a1 + a2.^a1)', 'options', s);

    [fitobj, gof] = fit(xdata_fit, ydata_fit, HillEqn);

    % Plot fit
    xfit = linspace(min(xdata_fit), max(xdata_fit), 500);
    yfit = feval(fitobj, xfit);
    plot(xfit, yfit, 'k-', 'LineWidth', 2)

    % Annotate subplot
    text(0.6, 0.2, sprintf('RMSE=%.5f', gof.rmse), ...
        'Units', 'normalized', 'FontSize', 10)
    title(genotype)
    xlabel(xlabelstr)
    ylabel('Fraction Nuclei On')
    ylim([0 1])
    xlim([0, max(xfit)])
    ax = gca;                   % Get current axes
    colormap(ax, cmap);  % Use a specific colormap per subplot
    %colorbar(ax);
    hold off
end

if useCumInt
    filelabelstr = 'Int';
else
    filelabelstr = 'Inst';
end

saveas(gcf, strcat(Path,'Presentations/PlacedImages/HillPlots',filelabelstr,genename,'2.svg'));

%% Plot probability for each position bin and save the Ka value

genotypeNames = fieldnames(results);
colors = {myblue,mymagenta,myyellow,mygreen};  % Distinct color per genotype

HillEqn = fittype('x.^a1./(x.^a1 + a2.^a1)', ...
    'independent', 'x', 'coefficients', {'a1', 'a2'});
fitopts = fitoptions('Method','NonlinearLeastSquares', ...
    'Startpoint', [1, 1], 'Lower', [0, 0], 'Upper', [Inf, Inf], ...
    'Robust', 'LAR');

figure
hold on
xlabel('Position (Microns)')
ylabel('Ka Value')
xlim([0,60])
ylim([0,0.5])
fontsize(14,"points")
for g = 1:3%numel(genotypeNames)

    KaValues = zeros(numel(binedges)-1,1);

    c = colors{g};
    genotype = genotypeNames{g};
    data = results.(genotype);

    for x = 1:numel(binedges) - 4
        bmp = data.CumInt;
        rna = data.RNAfraction;
        xdata = bmp(:,x);
        ydata = rna(:,x);
        [fitresult, gof] = fit(xdata, ydata, HillEqn, fitopts);

        ka = fitresult.a2;
        KaValues(x,1) = ka;
    end
    scatter(binCenters,KaValues,35,'MarkerFaceColor',c,'MarkerEdgeColor','none')
end

saveas(gcf, strcat(Path,'Presentations/PlacedImages/KaOverSpaceNoZen',genename,'.svg'));

%% Plot probability v. pSMAD levels with linear Hill fit

HillEqn = fittype('x.^a1./(x.^a1 + a2.^a1)', ...
    'independent', 'x', 'coefficients', {'a1', 'a2'});
fitopts = fitoptions('Method','NonlinearLeastSquares', ...
    'Startpoint', [1, 1], 'Lower', [0, 0], 'Upper', [Inf, Inf], ...
    'Robust', 'LAR');

genotypeNames = fieldnames(results);
colors = {myblue,mymagenta,myyellow,mygreen};  % Distinct color per genotype
figure; hold on;
fontsize(14,'points')
for g = 1:numel(genotypeNames)
    genotype = genotypeNames{g};
    fracMat = results.(genotype).FractionOnRaw_all;
    psmadMat = results.(genotype).MeanPSMADRaw_all;

    % Convert NaNs to 0 if not already done
    fracMat(isnan(fracMat)) = 0;
    psmadMat(isnan(psmadMat)) = 0;

    % Flatten matrices
    xdata = psmadMat(:);
    ydata = fracMat(:);

    % Remove any zero psmad to avoid log(0) issues in the fit (optional)
    valid = xdata > 0;
    xdata = xdata(valid);
    ydata = ydata(valid);

    % Scatter plot
    scatter(xdata, ydata, 25, 'MarkerFaceColor', colors{g}, ...
        'MarkerEdgeColor', 'none','MarkerFaceAlpha',0.5);

    % Fit Hill function
    try
        [fitresult, gof] = fit(xdata, ydata, HillEqn, fitopts);
        ci = confint(fitresult, 0.95);  % 95% confidence intervals

        % Extract confidence intervals for a1 and a2
        a1_CI = ci(:,1);
        a2_CI = ci(:,2);

        fprintf('%s:\n', genotype);
        fprintf('  Hill coefficient (a1): %.2f [%.2f, %.2f]\n', fitresult.a1, a1_CI(1), a1_CI(2));
        fprintf('  Half-max (a2):         %.4f [%.4f, %.4f]\n', fitresult.a2, a2_CI(1), a2_CI(2));
        fprintf('  R² = %.3f\n', gof.rsquare);

        results.(genotype).HillFit = fitresult;   % The fitted model
        results.(genotype).HillCI = ci;    % Struct with .a1 and .a2 confidence intervals

        % Plot fit
        xfit = linspace(min(xdata), max(xdata), 200);
        yfit = feval(fitresult, xfit);
        plot(xfit, yfit, 'Color', colors{g}, 'LineWidth', 2);

        % Annotate with R^2 in corner
        % text(0.05, 0.95 - 0.1*g, ...
        %     sprintf('%s: R^2=%.2f', genotype, gof.rsquare), ...
        %     'Units', 'normalized', 'Color', colors{g}, 'FontSize', 14);
    catch ME
        warning('Fit failed for %s: %s', genotype, ME.message);
    end
end

xlabel('Mean pSMAD level');
ylabel('Fraction On');
xlim([0, 1]);
ylim([0, 1]);

saveas(gcf, strcat(Path,'Presentations/PlacedImages/TimeIndHill',genename,'.svg'));


%%

genotypeNames = fieldnames(results);
heatmapCols = {mybluemap,mymagentamap,myyellowmap,mygreenmap};

% Define labels
xlabels = string(binedges(1,2:end));   % bin edges (length 16 for 15 bins, use centers if needed)
ylabels = string(1:60);     % time points

for i = 1:numel(genotypeNames)
    heatmapcolor = heatmapCols{i};  % should be a valid colormap (Nx3 RGB)
    genotype = genotypeNames{i};

    % Data to plot
    fixedmedRate = results.(genotype).Smooth;  % should be 60 x Nbins matrix
    fixedmedRateBinary = double(fixedmedRate > 0.01);
    figure;
    h = heatmap(xlabels, ylabels, fixedmedRateBinary);
    % Customize heatmap appearance
    h.Colormap = heatmapcolor;  
    h.XLabel = 'Microns';
    h.YLabel = 'Minutes';
    % Optional: Make axes in correct order
    h.YDisplayLabels = ylabels;
    h.XDisplayLabels = xlabels;
    h.GridVisible = false;
    h.ColorLimits = [0 1];

    saveas(gcf, strcat(Path,sprintf('Presentations/PlacedImages/FixedMadHeatmapBinary_%s.svg', genotype)));

    % Data to plot
    livemedRate = results.(genotype).Smooth;  % should be 60 x Nbins matrix
    figure;
    h = heatmap(xlabels, ylabels, livemedRate);
    % Customize heatmap appearance
    h.Colormap = heatmapcolor;  
    h.XLabel = 'Position';
    h.YLabel = 'Minutes';
    % Optional: Make axes in correct order
    h.YDisplayLabels = ylabels;
    h.XDisplayLabels = xlabels;
    h.GridVisible = false;
    h.ColorLimits = [0 0.05];
    saveas(gcf, strcat(Path,sprintf('Presentations/PlacedImages/MedHeatmap_%s.svg', genotype)));

end

%%

genotypeNames = fieldnames(results);
heatmapCols = {mybluemap,mymagentamap,myyellowmap,mygreenmap};

% Define labels
xlabels = string(binedges(1,2:end));   % bin edges (length 16 for 15 bins, use centers if needed)
ylabels = string(1:60);     % time points

for i = 1:numel(genotypeNames)
    heatmapcolor = heatmapCols{i};  % should be a valid colormap (Nx3 RGB)
    genotype = genotypeNames{i};

    % Data to plot
    fixedRNA = results.(genotype).RNAfraction;  % should be 60 x Nbins matrix
    figure;
    h = heatmap(xlabels, ylabels, fixedRNA);
    % Customize heatmap appearance
    h.Colormap = heatmapcolor;  
    h.XLabel = 'Position';
    h.YLabel = 'Minutes';
    % Optional: Make axes in correct order
    h.YDisplayLabels = ylabels;
    h.XDisplayLabels = xlabels;
    h.GridVisible = false;
    h.ColorLimits = [0 1];
    saveas(gcf, strcat(Path,sprintf('Presentations/PlacedImages/FixedRNAHeatmap_%s', genotype),genename,'.svg'));

    % Data to plot
    fixedRNABinary = double(fixedRNA > 0.1);
    figure;
    h = heatmap(xlabels, ylabels, fixedRNABinary);
    % Customize heatmap appearance
    h.Colormap = heatmapcolor;  
    h.XLabel = 'Position';
    h.YLabel = 'Minutes';
    % Optional: Make axes in correct order
    h.YDisplayLabels = ylabels;
    h.XDisplayLabels = xlabels;
    h.GridVisible = false;
    h.ColorLimits = [0 1];
    saveas(gcf, strcat(Path,sprintf('Presentations/PlacedImages/FixedRNAHeatmapBinary_%s', genotype),genename,'.svg'));

end

%%
groupSize = 3;  
nbins = length(binedges) - 1;
ngroups = 5;%ceil(nbins / groupSize);
useCumInt = true;

figure
for g = 1:nGenos
    genotype = genotypeNames{g};
    colorscheme = datasetCols{g};
    cmap = cbrewer2(colorscheme, 5);  % One color per group
    data = results.(genotype);

    if useCumInt
        XAxisMat = data.CumInt;
        xlabelstr = 'Integrated Med';
    else
        XAxisMat = data.Smooth;
        xlabelstr = 'Instantaneous Med';
    end

    Y = data.RNAfraction;
    X = [zeros(lagtime, nbins); XAxisMat(1:(end-lagtime), :)];

    subplot(nRows, nCols, g)
    hold on

    for groupIdx = 1:ngroups
        binStart = (groupIdx - 1) * groupSize + 1;
        binEnd = min(groupIdx * groupSize, nbins);

        xdata_group = [];
        ydata_group = [];

        for b = binStart:binEnd
            xdata = X(:, b);
            ydata = Y(:, b);
            valid = isfinite(xdata) & isfinite(ydata);
            xdata = xdata(valid);
            ydata = ydata(valid);

            scatter(xdata, ydata, 12, ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', cmap(ngroups-groupIdx+1, :), ...
                'MarkerFaceAlpha',0.25)
            xdata_group = [xdata_group; xdata];
            ydata_group = [ydata_group; ydata];
        end

        % Fit Hill function per group
        valid = isfinite(xdata_group) & isfinite(ydata_group);
        xfitdata = xdata_group(valid);
        yfitdata = ydata_group(valid);

        if numel(xfitdata) > 5  % Avoid crashing if very sparse
            startPoints = [1 1];
            s = fitoptions('Method','NonlinearLeastSquares','Robust','LAR', ...
                'Lower',[0 0],'Upper',[Inf Inf],'Startpoint',startPoints);
            HillEqn = fittype('x.^a1 ./ (x.^a1 + a2.^a1)', 'options', s);

            try
                [fitobj, gof] = fit(xfitdata, yfitdata, HillEqn);
                xplot = linspace(min(xfitdata), max(xfitdata), 500);
                yplot = feval(fitobj, xplot);
                plot(xplot, yplot, '-', 'Color', cmap(ngroups-groupIdx+1,:), 'LineWidth', 2)
            catch
                warning('Fit failed for group %d in %s', groupIdx, genotype);
            end
        end
    end

    title(genotype)
    xlabel(xlabelstr)
    ylabel('Fraction Nuclei On')
    ylim([0 1])
    xlim([0, 0.25])
    colormap(gca, cmap);
    hold off
end

if useCumInt
    filelabelstr = 'Int';
else
    filelabelstr = 'Inst';
end

saveas(gcf, strcat(Path,'Presentations/PlacedImages/HillPlotsGrouped',filelabelstr,genename,'.svg'));

groupSize = 3;  
nbins = length(binedges) - 1;
ngroups = 5;%ceil(nbins / groupSize);
useCumInt = false;

figure
for g = 1:nGenos
    genotype = genotypeNames{g};
    colorscheme = datasetCols{g};
    cmap = cbrewer2(colorscheme, 5);  % One color per group
    data = results.(genotype);

    if useCumInt
        XAxisMat = data.CumInt;
        xlabelstr = 'Integrated Med';
    else
        XAxisMat = data.Smooth;
        xlabelstr = 'Instantaneous Med';
    end

    Y = data.RNAfraction;
    X = [zeros(lagtime, nbins); XAxisMat(1:(end-lagtime), :)];

    subplot(nRows, nCols, g)
    hold on

    for groupIdx = 1:ngroups
        binStart = (groupIdx - 1) * groupSize + 1;
        binEnd = min(groupIdx * groupSize, nbins);

        xdata_group = [];
        ydata_group = [];

        for b = binStart:binEnd
            xdata = X(:, b);
            ydata = Y(:, b);
            valid = isfinite(xdata) & isfinite(ydata);
            xdata = xdata(valid);
            ydata = ydata(valid);

            scatter(xdata, ydata, 12, ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', cmap(ngroups-groupIdx+1, :), ...
                'MarkerFaceAlpha',0.25)
            xdata_group = [xdata_group; xdata];
            ydata_group = [ydata_group; ydata];
        end

        % Fit Hill function per group
        valid = isfinite(xdata_group) & isfinite(ydata_group);
        xfitdata = xdata_group(valid);
        yfitdata = ydata_group(valid);

        if numel(xfitdata) > 5  % Avoid crashing if very sparse
            startPoints = [1 1];
            s = fitoptions('Method','NonlinearLeastSquares','Robust','LAR', ...
                'Lower',[0 0],'Upper',[Inf Inf],'Startpoint',startPoints);
            HillEqn = fittype('x.^a1 ./ (x.^a1 + a2.^a1)', 'options', s);

            try
                [fitobj, gof] = fit(xfitdata, yfitdata, HillEqn);
                xplot = linspace(min(xfitdata), max(xfitdata), 500);
                yplot = feval(fitobj, xplot);
                plot(xplot, yplot, '-', 'Color', cmap(ngroups-groupIdx+1,:), 'LineWidth', 2)
            catch
                warning('Fit failed for group %d in %s', groupIdx, genotype);
            end
        end
    end

    title(genotype)
    xlabel(xlabelstr)
    ylabel('Fraction Nuclei On')
    ylim([0 1])
    xlim([0, 0.025])
    colormap(gca, cmap);
    hold off
end

if useCumInt
    filelabelstr = 'Int';
else
    filelabelstr = 'Inst';
end

saveas(gcf, strcat(Path,'Presentations/PlacedImages/HillPlotsGrouped',filelabelstr,genename,'.svg'));

%% Check for lag time

useCumInt = false;  % Use instantaneous signal
lagtimes = [0, 2, 5]; 
lagLabels = {'0 min', '2 min', '5 min'};
lagColors = cbrewer2('Blues',numel(lagtimes));  % Distinct colors for each lag

genotypeNames = fieldnames(results);
datasetCols2 = {myblue, mymagenta, myyellow, mygreen};

figure
hold on
legendEntries = {};

for g = 1%:numel(genotypeNames)
    genotype = genotypeNames{g};
    thisColor = datasetCols2{g};
    data = results.(genotype);

    % Select X-axis variable
    if useCumInt
        XAxisMat = data.CumInt;
        xlabelstr = 'Integrated Med';
    else
        XAxisMat = data.Smooth;
        xlabelstr = 'Instantaneous Med';
    end

    Y = data.RNAfraction;  % 60 x Nbins

    for l = 1:numel(lagtimes)
        lag = lagtimes(l);
        X = [zeros(lag, size(XAxisMat,2)); XAxisMat(1:(end-lag), :)];

        xdata = X(:);
        ydata = Y(:);

        valid = isfinite(xdata) & isfinite(ydata);
        xdata = xdata(valid);
        ydata = ydata(valid);

        % Fit Hill function
        startPoints = [1, 1];
        s = fitoptions('Method','NonlinearLeastSquares','Robust','LAR',...
            'Lower',[0 0],'Upper',[Inf Inf],'Startpoint',startPoints);
        HillEqn = fittype('x.^a1 ./ (x.^a1 + a2.^a1)', 'options', s);

        [fitobj, gof] = fit(xdata, ydata, HillEqn);

        % Plot scatter points
        scatter(xdata, ydata, 10, ...
            'MarkerEdgeColor', lagColors(l,:), ...
            'MarkerFaceColor', lagColors(l,:), ...
            'MarkerFaceAlpha', 0.2)

        % Plot fitted curve
        xfit = linspace(min(xdata), max(xdata), 1000);
        yfit = feval(fitobj, xfit);
        plot(xfit, yfit, '-', 'Color', lagColors(l,:), 'LineWidth', 2)

        % Optional: text with RMSE
        ypos = 0.2 - 0.05*l - 0.08*(g-1);
        text(0.6, ypos, ...
            sprintf('%s %s: RMSE=%.4f', genotype, lagLabels{l}, gof.rmse), ...
            'Color', lagColors(l,:), 'FontSize', 9, 'Units', 'normalized')
    end
end

xlabel(xlabelstr)
ylabel('Fraction Nuclei On')
ylim([0, 1])
xlim([0, max(xfit)])
fontsize(14, "points")
hold off

filelabelstr = useCumInt * "Int" + ~useCumInt * "Inst";
saveas(gcf, strcat(Path,'Presentations/PlacedImages/HillLagCompare_',filelabelstr,genename,'.svg'));
