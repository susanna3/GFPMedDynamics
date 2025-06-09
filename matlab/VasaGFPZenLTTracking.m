% Plan for zenLT tracking

%% Load Data
Path = 'C:/Users/Susanna Brantley/Brantley Lab Dropbox/Susanna Brantley/Duke/DiTaliaLab/ImageAnalysis/Datasets/';
%Path = 'D:/Dropbox/Duke/DiTaliaLab/ImageAnalysis/Datasets/';
%Path = '/Users/susannabrantley/Dropbox/Duke/DiTaliaLab/ImageAnalysis/Datasets/';

%ExpFolder = 'ZenLT_BcdGFP_HisRFP/';

s.matlab.fonts.editor.code.Name.TemporaryValue = 'Arial';
set(0,'DefaultAxesFontName','Arial');

% Load all embryo data and organize into a structured array

AllDataStructured.NC14 = { ...
    readmatrix(strcat(Path,"250120_zenLT_vasaGFP_HisRFP_1_nc13_cleaned.csv")); ...
    readmatrix(strcat(Path,"250210_zenLT_vasaGFP_HisRFP_1_nc13_cleaned.csv")); ...
    readmatrix(strcat(Path,"250210_zenLT_vasaGFP_HisRFP_3_nc13_cleaned.csv")) ...
};

AllDataStructured.NC13 = { ...
    readmatrix(strcat(Path,"250210_zenLT_vasaGFP_HisRFP_1_nc12_cleaned.csv")); ...
    readmatrix(strcat(Path,"250210_zenLT_vasaGFP_HisRFP_3_nc12_cleaned.csv")) ...
};

AllDataStructured.NC12 = { ...
    readmatrix(strcat(Path,"250210_zenLT_vasaGFP_HisRFP_1_nc11_cleaned.csv")); ...
    readmatrix(strcat(Path,"250210_zenLT_vasaGFP_HisRFP_3_nc11_cleaned.csv")) ...
};

% Example: Accessing the data
NC14_Data = AllDataStructured.NC14;  % Cell array of NC14 data
NC13_Data = AllDataStructured.NC13;  % Cell array of NC13 data
NC12_Data = AllDataStructured.NC12;  % Cell array of NC12 data


%% Track nuclei
% save nuc, cyto, and position

timeoffset = [9,9,9,4,4,0,0];
midline_precalc = [352,442,384,442,384,442,384];

stages = fieldnames(AllDataStructured);

% Loop over each NC stage
for s = 1:length(stages)
    stage_name = stages{s}; % Get the stage name (e.g., 'NC14', 'NC13', 'NC12')
    datasets = AllDataStructured.(stage_name); % Extract datasets for this stage
    
    fprintf('Processing %s (Number of datasets: %d)\n', stage_name, length(datasets));

    AllZenNucData = [];
    AllZenRatioData = [];
    AllZenXData = [];
    AllZenYData = [];

    for embryo = 1:length(datasets)
        ZenData = datasets{embryo}; % Extract matrix
        %Separate data by frame
        [~,~,X] = unique(ZenData(:,1));
        frame_array = accumarray(X,1:size(ZenData,1),[],@(r){ZenData(r,:)});
        Timepoints = length(frame_array);

        % ImageNo | ObjectNo | X | Y | Nuc | Cyto | ScaleFactor | NCRatio |
        % NucRescale

        X_index = 3;
        Y_index = 4;
        ObjectNum_index = 2;
        scale_factor = 7;
        Nuc_index = 9;
        Ratio_index = 8;
        TimeDistanceThreshold = 50;

        [NucleiOverTime,TrackedIDArray,TrackedNucZenArray, TrackedRatioArray, TrackedXArray, TrackedYArray] = tracknuclei_zenLT( ...
            frame_array,TimeDistanceThreshold, X_index,Y_index,ObjectNum_index,Nuc_index,Ratio_index,scale_factor);


        % Find the midline at last time point
        TrackedRatioArray(any(TrackedRatioArray == 0, 2), :) = [];
        TrackedYArray(any(TrackedYArray == 0, 2), :) = [];
        TrackedYArray(any(TrackedRatioArray > 3, 2), :) = [];
        TrackedRatioArray(any(TrackedRatioArray > 3, 2), :) = [];

        % Fit the data using a smoothing spline
        ft = fit(TrackedYArray(:,end), TrackedRatioArray(:,end), 'SmoothingSpline', 'SmoothingParam', 0.0000005);

        % Generate a fine range of x-values for smooth evaluation
        y_fine = linspace(min(TrackedYArray(:,end)), max(TrackedYArray(:,end)), 1000);
        y_fitted = feval(ft, y_fine);

        % Find the maximum value and its corresponding x
        [max_value, midline] = max(y_fitted)
        max_x = y_fine(midline);

        % Plot
        figure
        hold on
        scatter(TrackedYArray(:,end), TrackedRatioArray(:,end));
        plot(y_fine, y_fitted, 'r', 'LineWidth', 2);
        plot(max_x, max_value, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b'); % Mark max point
        ylim([0,2])
        hold off

        midline = midline_precalc(1,embryo);

        % Specify the midline distance for all nuclei
        Midline = abs(TrackedYArray(:,end)-midline)./4.4;
        TrackedRatioArray = [Midline,TrackedRatioArray];

        AllZenNucData = [AllZenNucData;TrackedNucZenArray];
        AllZenRatioData = [AllZenRatioData;TrackedRatioArray];
        AllZenXData = [AllZenXData;TrackedXArray];
        AllZenYData = [AllZenYData;TrackedYArray];

    end
    AllStagesData.(stage_name) = AllZenRatioData;
end

%% Bin all data by position

for stage = 1:length(stages)
    stage_name = stages{stage};
    RatioData = AllStagesData.(stage_name);

    num_timepoints = size(RatioData, 2);

    distbins = 0:3:60;
    binnums = discretize(RatioData(:,1),distbins);
    
    BinnedZen = [];

    for bin = 1:length(distbins)-1
        nucsinbin = RatioData(binnums(:,1)==bin,2:end);
        avgpert = mean(nucsinbin,1);
        BinnedZen = [BinnedZen;avgpert];
    end

    AllBinnedData.(stage_name) = BinnedZen;
end


    
%% Plot NC ratio versus distance from midline

% Define colormaps for each stage using cbrewer2
stage_colors = {'PuRd', 'Blues', 'Greens'}; % Different gradients for each stage
colormap_sizes = [50, 10, 5]; % Adjust based on number of rows

figure;
hold on;
xlabel('Distance (microns)');
ylabel('NC Zen-LT Ratio');
fontsize(14,"points");

for stage = 1:length(stages)
    stage_name = stages{stage};
    RatioData = AllBinnedData.(stage_name);

    % Define distance bins
    distbins = 0:3:60; % Adjust if needed
    distance = distbins(1:end-1) + diff(distbins)/2; % Bin centers

    % Generate colormap for the current stage
    colormap_stage = cbrewer2(stage_colors{stage}, colormap_sizes(stage));

    % Get timepoints
    num_timepoints = size(RatioData, 2);

    % Plot each timepoint as a separate fitted curve
    for t = 1:num_timepoints
        % Fit a smoothing spline
        fit_curve = fit(distance(:), RatioData(:,t), 'smoothingspline', 'SmoothingParam', 0.0001);

        % Generate fine distance points for smooth plotting
        distance_fine = linspace(min(distance), max(distance), 200);
        fit_values = feval(fit_curve, distance_fine);

        % Plot the fitted curve
        plot(distance_fine, fit_values, 'Color', colormap_stage(t,:), 'LineWidth', 1.5);
    end
end

legend(stages, 'Location', 'Best'); % Add legend for stages
hold off;


%%

% Load colormap using cbrewer2

for stage = 1:length(stages)
    stage_name = stages{stage};
    BinnedRatio = AllBinnedData.(stage_name);

    colormap_blue = cbrewer2('Blues', 20);  % 20 shades of blue for distances

    % Define time and distance axes
    num_rows = size(BinnedRatio, 1); % Number of distance bins (20)
    num_cols = size(BinnedRatio, 2); % Number of time points (50)

    time = linspace(0, (num_cols - 1), num_cols); % Time in minutes
    distance = linspace(0, (num_rows - 1) * 3, num_rows); % Distance in microns

    % ---- Plot fitted curves over time (Each row/distance as a separate curve) ----
    figure;
    hold on;
    fontsize(14,"points");
    for i = 1:num_rows
        % Fit a smoothing spline to each row (distance)
        fit_curve = fit(time(:), BinnedRatio(i,:)', 'smoothingspline', 'SmoothingParam', 0.0001);

        % Generate smooth time points for plotting
        time_fine = linspace(min(time), max(time), 200);
        fit_values = feval(fit_curve, time_fine);

        % Plot the fitted curve
        plot(time_fine, fit_values, 'Color', colormap_blue(21-i,:), 'LineWidth', 1.5);
    end
    xlabel('Time (minutes)');
    ylabel('Zen-LT NC Ratio');
    ylim([1,1.5])

    % Add colorbar to indicate distance
    colormap(colormap_blue);
    c = colorbar;
    c.Label.String = 'Distance (microns)';
    c.Ticks = linspace(0, 1, 3); % Adjust as needed
    c.TickLabels = round(linspace(max(distance), min(distance), 3)); % Show actual distance values

    hold off;
end


%% Plot MedGFP gradient and/or RNA gradient for ush versus zenLT

% First run SogPrediction_Live_NoRNAFit.m to get BinnedMed1 

BinnedRatio = AllBinnedData.NC14;
BinnedMedCompare = BinnedMed1(2:end,11:end);

% Define distance axis
num_rows = size(BinnedRatio, 1); % Number of distance bins (20)
distance = linspace(0, (num_rows - 1) * 3, num_rows); % Distance in microns

% Select columns to plot
cols_to_plot = [20, 35, 50]; 
num_cols = length(cols_to_plot);

% Get color maps
colors_ratio = cbrewer2('Blues', num_cols); % Different shades of PuRd for BinnedRatio
colors_med = cbrewer2('Blues', num_cols);  % Different shades of Blues for BinnedMed1

figure;
hold on;

% ---- Left Y-Axis for BinnedRatio (Fitted) ----
yyaxis left
for i = 1:num_cols
    col = cols_to_plot(i);
    
    % Fit smoothing spline
    fit_ratio = fit(distance(:), BinnedRatio(:, col), 'smoothingspline', 'SmoothingParam', 0.0001);
    
    % Generate fine points for smooth curve
    dist_fine = linspace(min(distance), max(distance), 200);
    fit_curve_ratio = feval(fit_ratio, dist_fine);
    
    % Plot fitted curve with PuRd colors
    plot(dist_fine, fit_curve_ratio, '-', 'Color', colors_ratio(i,:), 'LineWidth', 2, ...
        'DisplayName', ['BinnedRatio Fit Col ', num2str(col)]);
end
ylabel('BinnedRatio');
ylim([min(BinnedRatio(:)), max(BinnedRatio(:))]); % Adjust limits for clarity
ax = gca;
ax.YColor = 'k'; % Set y-axis color to black for visibility

% ---- Right Y-Axis for BinnedMedCompare (Fitted) ----
yyaxis right
for i = 1:num_cols
    col = cols_to_plot(i);
    
    % Fit smoothing spline
    fit_med = fit(distance(:), BinnedMedCompare(:, col), 'smoothingspline', 'SmoothingParam', 0.01);
    
    % Generate fine points for smooth curve
    fit_curve_med = feval(fit_med, dist_fine);
    
    % Plot fitted curve with Blues colors
    plot(dist_fine, fit_curve_med, '--', 'Color', colors_med(i,:), 'LineWidth', 2, ...
        'DisplayName', ['BinnedMedCompare Fit Col ', num2str(col)]);
end
ylabel('BinnedMedCompare');
ylim([min(BinnedMedCompare(:)), max(BinnedMedCompare(:))]); % Adjust limits for clarity
ax.YColor = 'k'; % Set y-axis color to black for visibility

% ---- Labels and Legend ----
xlabel('Distance (microns)');
title('Fitted Curves for BinnedRatio (PuRd) and BinnedMedCompare (Blues)');
legend('Location', 'best');
hold off;


