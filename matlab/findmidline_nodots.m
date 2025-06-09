%Find the midline and convert Y-coord to distance from midline
%Then corrects for background based on cells far from the midline

function [CorrectedNuclei, FullyTrackedNuclei] =  findmidline_nodots(FullyTrackedY,FullyTrackedNuclei, ...
    NumTimepoints,StartTime,MidlineTimepoint,PixelCutoffLow,PixelCutoffHigh,PixelToMicron,BackgroundCutoffLow,BackgroundCutoffHigh,...
    Y_index,MicronY_index,MedIntensity_index)
%%

MidlineOverTime = [zeros(StartTime,2)];

figure
hold on
for t = 1:NumTimepoints;
    distanceData = FullyTrackedY(:,t+1);
    medData = FullyTrackedNuclei(:,t+1);
    dataset = [distanceData medData];
    %dataset = sortrows(dataset,1);
    filterbyy = dataset(dataset(:,1)>PixelCutoffLow & dataset(:,1)<PixelCutoffHigh,:);
    yData = filterbyy(:,2);
    xData = filterbyy(:,1);
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.00005;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    yfitted = feval(fitresult,xData);
    plot(xData,yfitted)
    scatter(xData,yData)
    fittedmed = [yfitted xData];
    amplitude = max(yfitted);
    midlineinfo = fittedmed(fittedmed(:,1)==amplitude,:);
    MidlineOverTime = [MidlineOverTime; midlineinfo];
end
hold off

%% correct distance into microns from midline

MidlineCoord = MidlineOverTime(MidlineTimepoint,2);
YAtMaxGradient = FullyTrackedY(:,MidlineTimepoint-StartTime);
DistToMid = abs(YAtMaxGradient-MidlineCoord)./PixelToMicron;
FullyTrackedNuclei(:,1) = DistToMid;

%% Find background intensity based on average intensity near the midline in first frame
background_datapoint = FullyTrackedNuclei(:,(25-StartTime+1));
background_med = mean(mean(background_datapoint))

CorrectedNuclei = FullyTrackedNuclei(:,1);

for t = 1:NumTimepoints
    CorrectedMed = [];
    background_dataset = FullyTrackedNuclei(:,t+1);
    %background_med_nuclei = background_dataset(DistToMid(:,1)>BackgroundCutoffLow & DistToMid(:,1)<BackgroundCutoffHigh,:);
    %background_med = mean(background_med_nuclei);
    [num_row num_col] = size(background_dataset);
    for row = 1:num_row
        old_med = background_dataset(row,1);
        corrected_med = old_med - background_med;
        if corrected_med < 0
            corrected_med = 0;
            CorrectedMed = [CorrectedMed; corrected_med];
        else
            CorrectedMed = [CorrectedMed; corrected_med];
        end
    end
    CorrectedNuclei = [CorrectedNuclei CorrectedMed];
end

CorrectedNuclei = sortrows(CorrectedNuclei,1);
end