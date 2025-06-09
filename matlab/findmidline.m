%Find the midline and convert Y-coord to distance from midline
%Then corrects for background based on cells far from the midline

function [CorrectedNuclei, CorrectedDots, YAtMaxGradient, XAtMaxGradient] =  findmidline(FullyTrackedY,FullyTrackedX,FullyTrackedNuclei,FullyTrackedDotCount, ...
    NumTimepoints,PixelCutoffLow,PixelCutoffHigh,PixelToMicron,BackgroundCutoffLow,BackgroundCutoffHigh)

%% 
MidlineOverTime = [];

%figure
%hold on
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
    opts.SmoothingParam = 0.000005;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    yfitted = feval(fitresult,xData);
    %plot(xData,yfitted)
    fittedmed = [yfitted xData];
    amplitude = max(yfitted);
    midlineinfo = fittedmed(fittedmed(:,1)==amplitude,:);
    MidlineOverTime = [MidlineOverTime; midlineinfo];
end
%hold off

%% correct distance into microns from midline

MidlineCoord = MidlineOverTime(NumTimepoints-10,2);
YAtMaxGradient = FullyTrackedY(:,NumTimepoints-10);
XAtMaxGradient = FullyTrackedX(:,NumTimepoints-10);
DistToMid = abs(YAtMaxGradient-MidlineCoord)./PixelToMicron;
FullyTrackedNuclei(:,1) = DistToMid;

CorrectedNuclei = FullyTrackedNuclei(:,1);

%% correct for background based on average intensity of cells far from midline
for t = 1:NumTimepoints
    CorrectedMed = [];
    background_dataset = FullyTrackedNuclei(:,t+1);
    background_med_nuclei = background_dataset(DistToMid(:,1)>BackgroundCutoffLow & DistToMid(:,1)<BackgroundCutoffHigh,:);
    background_med = mean(background_med_nuclei);
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

FullyTrackedDotCount(:,1) = DistToMid;
CorrectedDots = sortrows(FullyTrackedDotCount,1);
CorrectedNuclei = sortrows(CorrectedNuclei,1);
CorrectedNuclei = CorrectedNuclei(CorrectedNuclei(:,1)<50,:);
CorrectedDots = CorrectedDots(CorrectedDots(:,1)<50,:);

YAtMaxGradient = (YAtMaxGradient-MidlineCoord)./PixelToMicron;
XAtMaxGradient = XAtMaxGradient./PixelToMicron;
%% 
end