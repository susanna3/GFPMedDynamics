
Path = 'C:/Users/Susanna Brantley/Brantley Lab Dropbox/Susanna Brantley/Duke/DiTaliaLab/ImageAnalysis/Datasets/';
%Path = 'D:/Dropbox/Duke/DiTaliaLab/ImageAnalysis/Datasets/';
%Path = '/Users/susannabrantley/Dropbox/Duke/DiTaliaLab/ImageAnalysis/Datasets/';

%Read in data file for Windows for reference embryo
D1_binned = readmatrix(strcat(Path,"BinnedRefEmbryo_D1binned_221228_008.csv"));

%File Labels:
%ImageNo | ObjectNo | Area | Int | X | Y | Min | Gen | Z 

%Load all embryo data
%Embryohnt1 = readmatrix("D:/Dropbox/Duke/DiTaliaLab/ImageAnalysis/Datasets/230119_003_wt_GFPMed_HisiRFP_cellprofiler_dotscleaned.csv");
Embryo1 = readmatrix(strcat(Path,"MedZenPaper/live_wt_gfpmed_ushms2_240502_1_cleaned.csv"));
Embryo2 = readmatrix(strcat(Path,"MedZenPaper/live_wt_gfpmed_ushms2_240429_3_cleaned.csv"));
Embryo3 = readmatrix(strcat(Path,"MedZenPaper/live_wt_gfpmed_ushms2_240429_2_cleaned.csv"));
Embryo4 = readmatrix(strcat(Path,"MedZenPaper/live_wt_gfpmed_ushms2_240429_1_cleaned.csv"));
Embryo5 = readmatrix(strcat(Path,"MedZenPaper/live_wt_gfpmed_ushms2_230104_001_cleaned.csv"));
Embryo6 = readmatrix(strcat(Path,"MedZenPaper/live_wt_gfpmed_ushms2_221228_021_cleaned.csv"));
Embryo7 = readmatrix(strcat(Path,"MedZenPaper/live_wt_gfpmed_ushms2_221228_008_009_cleaned.csv"));

EmbryoData = {Embryo1;Embryo2;Embryo3;Embryo4;Embryo5;Embryo6;Embryo7};

Bins = [11 15 15 10 9 12];
LowCutoff = [350 300 300 350 275 150];
HighCutoff = [725 750 800 750 700 850];
BackgroundLow = [45 55 55 40 55 55]; 
BackgroundHigh = [55 65 65 45 65 65];

%hntinfo
% Bins = [11];
% LowCutoff = [300];
% HighCutoff = [700];
% BackgroundLow = [45]; 
% BackgroundHigh = [55];

NumEmbryos = 7;

%Set shared parameters
NumTimepoints = 35;
NumZs = 11;
DistanceThreshold = 3.2;
TimeDistanceThreshold = 30;
PixelToMicron = 4.4;
StartTime = 25;
EndTime = 59;
Micron_cutoff = 75;

%%
WidthPerEmbryo = [];
AmplitudePerEmbryo = [];
AllTrackedNuclei = [];
AllTrackedDots = [];
FirstTranscript = {};

for embryo = 1:NumEmbryos
    NumBins = Bins(1,embryo);
    PixelCutoffLow = LowCutoff(1,embryo);
    PixelCutoffHigh = HighCutoff(1,embryo);
    BackgroundCutoffLow = BackgroundLow(1,embryo);
    BackgroundCutoffHigh = BackgroundHigh(1,embryo);
    MedData = EmbryoData{embryo};
    %Separate data by frame
    [~,~,X] = unique(MedData(:,1));
    frame_array = accumarray(X,1:size(MedData,1),[],@(r){MedData(r,:)});

    %set column indexes for tracking in z
    ObjectNum_index = 2;
    Area_index = 3;
    MedIntensity_index = 4;
    X_index = 5;
    Y_index = 6;
    DotCount_index = 7;
    DotInt_index = 8;

    NextNuc_index = 12;
    NextArea_index = 13;
    NextMed_index = 14;
    NextX_index = 15;
    NextY_index = 16;
    NextDotCount_index = 17;
    NextDotInt_index = 18;

    [StackedNucleusArray,StackedAreaArray,StackedMedArray,StackedXArray,StackedYArray,StackedDotCountArray,StackedDotIntArray] = stackbyz( ...
    frame_array,NumTimepoints,NumZs,DistanceThreshold, ...
    X_index,Y_index,ObjectNum_index,MedIntensity_index,DotCount_index,DotInt_index,Area_index, ...
    NextNuc_index,NextArea_index,NextMed_index,NextX_index,NextY_index,NextDotCount_index,NextDotInt_index);

    NucleiStacked = avgperstack( ...
    StackedNucleusArray,StackedAreaArray,StackedMedArray,StackedXArray,StackedYArray,StackedDotCountArray,StackedDotIntArray,...
    NumTimepoints);

    %set coordinates for tracking in time
    ObjectNum_index = 1;
    Area_index = 2;
    MedIntensity_index = 4;
    X_index = 5;
    Y_index = 6;
    DotCount_index = 7;
    DotInt_index = 8;

    NextNuc_index = 9;
    NextMed_index = 10;
    NextX_index = 11;
    NextY_index = 12;
    NextDotCount_index = 13;
    NextDotInt_index = 14;

    [NucleiOverTime,TrackedNucleusArray,TrackedMedArray,TrackedXArray,TrackedYArray,TrackedDotCountArray,TrackedDotIntArray] = tracknuclei( ...
    NucleiStacked,TimeDistanceThreshold,NumTimepoints,PixelCutoffLow,PixelCutoffHigh, ...
    X_index,Y_index,ObjectNum_index,MedIntensity_index,DotCount_index,DotInt_index, ...
    NextNuc_index,NextMed_index,NextX_index,NextY_index,NextDotCount_index,NextDotInt_index);

    [FullyTrackedNuclei,FullyTrackedX,FullyTrackedY,FullyTrackedDotCount,FullyTrackedDotInt] = tracktoend( ...
    EndTime,NumTimepoints, ...
    TrackedYArray,TrackedMedArray,TrackedXArray,TrackedDotCountArray,TrackedDotIntArray);

    [CorrectedNuclei,CorrectedDots] =  findmidline(FullyTrackedY,FullyTrackedNuclei,FullyTrackedDotCount, ...
    NumTimepoints,PixelCutoffLow,PixelCutoffHigh,PixelToMicron,BackgroundCutoffLow,BackgroundCutoffHigh);

    [D2_binned,FullyTrackedNucleiRescaled] = rescale(D1_binned,CorrectedNuclei,NumBins,NumTimepoints);

    [AmplitudeOverTime,WidthOverTime] = widthamp(NumTimepoints,FullyTrackedNucleiRescaled);
    WidthPerEmbryo = [WidthPerEmbryo WidthOverTime];
    AmplitudePerEmbryo = [AmplitudePerEmbryo AmplitudeOverTime];

    %AmplitudeOverDist = ampofdistbins(NumTimepoints,FullyTrackedNucleiRescaled,StartTime,EndTime,NumBins);

    first_transcript = firsttranscript(FullyTrackedNucleiRescaled, CorrectedDots, ...
    NumTimepoints,StartTime,EndTime);
    FirstTranscript = [FirstTranscript; first_transcript];
    NucleiInfoPerMovie = [first_transcript(:,1:3),FullyTrackedNucleiRescaled];
    AllTrackedNuclei = [AllTrackedNuclei;NucleiInfoPerMovie];
    AllTrackedDots = [AllTrackedDots;CorrectedDots];

    frame_array={};
end
%% Plot fitted curves
cmap=cbrewer2('Blues',55);
figure
hold on
title('Exponential Fit');
for t = 1:7:NumTimepoints;
    hold on
    xData = FullyTrackedNucleiRescaled(:,1);
    yData = FullyTrackedNucleiRescaled(:,t+1);
    xData = xData(yData ~= 0);    
    yData = yData(yData ~= 0);

    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.0005;
    [fitresult, gof] = fit( xData, yData, ft, opts );

    %[ft, stats] = fit(xData,yData,'exp1');
    yfitted = feval(fitresult,xData);
    c = cmap(t+10,:);
    sz = 10;
    shade = t+20;
    scatter(xData,yData,sz,c);
    plot(xData,yfitted,'-','color',cmap(shade,:));
end
hold off

%% Plot Width and Amplitude 
cmap=cbrewer2('Blues',15);
figure
hold on
title("Amplitude Over Time");
ylabel("Amplitude (max of fitted curve)");
xlabel("Time (Minutes)");
for i = 1:NumEmbryos;
     yData = AmplitudePerEmbryo(:,i);
     xData = (StartTime:EndTime);
     shade = 5+(i*2);
     plot(xData,yData,'-','color',cmap(shade,:));
end
hold off

figure
hold on
title("Width Over Time");
ylabel("Width = Avg Width at 40,50,&60%");
xlabel("Time (Minutes)");
for i = 1:2:NumEmbryos*2;
    yData = WidthPerEmbryo(:,i);
    xData = (StartTime:EndTime);
    shade = 5+i;
    plot(xData,yData,'-','color',cmap(shade,:));
end
hold off 

figure
hold on
title("Width Over Time");
ylabel("Width = Integral/Amplitude");
xlabel("Time (Minutes)");
for i = 2:2:NumEmbryos*2;
    yData = WidthPerEmbryo(:,i);
    xData = (StartTime:EndTime);
    shade = 5+i;
    plot(xData,yData,'-','color',cmap(shade,:));
    ylim([0 75])
end
hold off 
%% Check X cutoff

cmap=cbrewer2('Blues',56);
figure
hold on
for mat = 1:10:NumTimepoints
    %figure
    xData = NucleiStacked{mat}(:,6)./3.2;
    %xData = xData./PixelToMicron;
    yData = NucleiStacked{mat}(:,4);
    c = cmap(mat+10,:);
    sz = 10;
    scatter(xData,yData,sz,c);
    ylim([0 .20]);
    %xlim([200 800]);
end
hold off

%% Plot Med threshold and on time

DistFromMid_transcribingnuc = [];
MedLevel_transcribingnuc = [];
TimeOn_allnuc = [];
MedIntegral_transcribingnuc = [];

for e = 1:NumEmbryos
    first_transcript = FirstTranscript{e};
    MS2_ontime = (first_transcript(:,1))+StartTime-1;
    Transcribing_Nuclei = first_transcript(first_transcript(:,1)<31,:);
    Med_level_at_on = Transcribing_Nuclei(:,2);
    Med_integral_at_on = Transcribing_Nuclei(:,3);
    dist_from_mid = Transcribing_Nuclei(:,5);

    DistFromMid_transcribingnuc = [DistFromMid_transcribingnuc;dist_from_mid];
    MedLevel_transcribingnuc = [MedLevel_transcribingnuc; Med_level_at_on];
    TimeOn_allnuc = [TimeOn_allnuc; MS2_ontime];
    MedIntegral_transcribingnuc = [MedIntegral_transcribingnuc; Med_integral_at_on];
end

figure
hold on
%ylim([0 50])
title("Med Threshold at Promoter Firing in Transcribing Nuclei")
histogram(MedLevel_transcribingnuc,'BinEdges', [0 0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1])
hold off

%     figure
%     hold on
%     title("Med Integral at Promoter Firing in Transcribing Nuclei")
%     histogram(Med_integral_at_on,'BinLimits',[0.00001,0.01],'BinWidth',0.0005)
%     hold off

% figure
% hold on
% ylim([0 100])
% title("Time To First Promoter Firing")
% histogram(TimeOn_allnuc)
% hold off
% 
% figure
% hold on
% title("Dist of Transcribing Nuclei")
% histogram(DistFromMid_transcribingnuc);
% hold off

figure
hold on
fontsize(gcf,14,'points')
scatter(DistFromMid_transcribingnuc,MedLevel_transcribingnuc)
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'SmoothingSpline');
opts.SmoothingParam = 0.0005;
[fitresult, gof] = fit( DistFromMid_transcribingnuc, MedLevel_transcribingnuc, ft );

plot((0:75),feval(fitresult,(0:75)),'Color','magenta','LineWidth',2)
%title("Med@ON v. Distance")
ylim([0 0.05])
xlabel("Microns")
ylabel("Med level at promoter firing")

figure
hold on
fontsize(gcf,14,'points')
scatter(DistFromMid_transcribingnuc,TimeOn_allnuc(TimeOn_allnuc(:,1)<55,1))
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'SmoothingSpline');
opts.SmoothingParam = 0.0005;
[fitresult, gof] = fit( DistFromMid_transcribingnuc, TimeOn_allnuc(TimeOn_allnuc(:,1)<55,1), ft );

plot((0:60),feval(fitresult,(0:60)),'Color','magenta','LineWidth',2)
%title("Time On v. Distance")
xlabel("Microns")
ylabel("Minutes")

figure
hold on
fontsize(gcf,14,'points')
scatter(DistFromMid_transcribingnuc,MedIntegral_transcribingnuc)
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'SmoothingSpline');
opts.SmoothingParam = 0.0005;
[fitresult, gof] = fit( DistFromMid_transcribingnuc, MedIntegral_transcribingnuc, ft );

plot((0:75),feval(fitresult,(0:75)),'Color','magenta','LineWidth',2)
%title("Integral v. Distance")
ylim([0 0.5])
xlabel("Microns")
ylabel("Integrated Med")

%% Plot relationship between Hill function and Distance

%BinnedZenWTData = load('ZenMadModelData/BinnedWTRNAZen.mat');
%BinnedZenWTData = BinnedZenWTData.BinnedRNA;
BinnedZenWTData = BinnedRNA;

% Bin MS2 nuclei so can assign zen value

inithillco = 3;
initKd = 0.000175;

NucleiInfo = [DistFromMid_transcribingnuc,TimeOn_allnuc(TimeOn_allnuc(:,1)<55,1)];
DistBins = discretize(NucleiInfo(:,1),[0:binsize:60,200]);
TimeBins = discretize(NucleiInfo(:,2),[0:5:60]);
NucleiInfo = [NucleiInfo,DistBins,TimeBins];

ZenConc = [];
IntZenConc = [];

for n = 1:length(NucleiInfo)
    distbin = NucleiInfo(n,3);
    timebin = NucleiInfo(n,4);
    zenconc = BinnedZenWTData(timebin,distbin);
    Q = BinnedZenWTData(1:timebin,distbin);
    Q(isnan(Q))=0;
    intzen = sum(Q);
    ZenConc = [ZenConc;zenconc];
    IntZenConc = [IntZenConc;intzen];
end

% Still need to think about this....

DistanceValues = DistFromMid_transcribingnuc(IntZenConc(:,1)>0,:);
MedIntValues = MedIntegral_transcribingnuc(IntZenConc(:,1)>0,:);
IntZenConc = IntZenConc(IntZenConc(:,1)>0,:);

MedModel = (MedIntValues.^inithillco)./(initKd + (MedIntValues.^inithillco));

figure
hold on
fontsize(gcf,14,'points')
scatter(DistanceValues,MedModel)
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'SmoothingSpline');
opts.SmoothingParam = 0.0005;
[fitresult, gof] = fit( DistanceValues, MedModel, ft );

plot((0:75),feval(fitresult,(0:75)),'Color','magenta','LineWidth',2)
%title("Integral v. Distance")
%ylim([0 0.000001])
xlabel("Microns")
ylabel("Zen RNA")

%% Determine if threshold is the same across time

% FirstTranscript_e1 = FirstTranscript{2};
% FirstTranscript_e2 = FirstTranscript{3};
% TranscribingNuclei = [FirstTranscript_e1;FirstTranscript_e2];

BinnedOnTime = discretize(Transcribing_Nuclei(:,1)+24, [25 30 35 40 45 50 55]);
MedBinnedByTime = [BinnedOnTime Transcribing_Nuclei(:,1) Med_level_at_on];
[~,~,X] = unique(MedBinnedByTime(:,1));

[num_total_nuclei width_total_nuc] = size(MedBinnedByTime);

BinnedMedValues = {};
for bin = 2:4
    figure
    hold on
    ylim([0 70])
    binned_info = [];
    for nuc = 1:num_total_nuclei
        check_bin = MedBinnedByTime(nuc,1);
        if check_bin == bin
            nuc_threshold = MedBinnedByTime(nuc,3);
            binned_info = [binned_info; nuc_threshold];
        else
            continue
        end
    end
    h = histogram(binned_info, 'BinEdges', [0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1]);
    hold off
    BinnedMedValues = [BinnedMedValues; binned_info];
end

Min30_35 = BinnedMedValues{1};
Min35_40 = BinnedMedValues{2};
Min40_45 = BinnedMedValues{3};

[h p k2stat] = kstest2(Min30_35,Min35_40)
[h2 p2 k2stat2] = kstest2(Min35_40,Min40_45)
[h3 p3 k2stat3] = kstest2(Min30_35,Min40_45)

%% Probabilty of promoter firing

first_transcript = [];
[num_tracked_nuclei timepoints_tracked] = size(FullyTrackedNucleiRescaled);

NucForProb = sortrows(FullyTrackedNucleiRescaled,1);
DotsForProb = sortrows(CorrectedDots,1);

NucPerTime = [];
TranscribingPerTime = [];
%BinsPerTime = [];

for t = 1:NumTimepoints
    nuclei = NucForProb(:,t+1);
    dots = DotsForProb(:,t+1);
    [numnuc w] = size(nuclei);
    medbins = discretize(nuclei,(0:0.001:0.15),'IncludedEdge','left');
    %max(medbins)
    nuclei = [medbins nuclei];
    nucon = [];

    nucperbin = zeros([1 150]);
    transcribingperbin = zeros([1 150]);
    for b = 1:max(medbins)
        totalnuc = 0;
        totaltranscribing = 0;
        for n = 1:numnuc
            bin = nuclei(n,1);
            if bin == b
                dotstatus = dots(n,1);
                if dotstatus==0
                    totalnuc = totalnuc+1
                    med = NucForProb(n,t+1)
                    dot = DotsForProb(n,t+1)
                elseif dotstatus>0
                    totalnuc= totalnuc+1;
                    totaltranscribing = totaltranscribing+1;
                    nucon = [nucon;n];
                    med = NucForProb(n,t+1);
                    dot = DotsForProb(n,t+1);
                end
            end
        end
        nucperbin(1,b) = totalnuc;
        transcribingperbin(1,b) = totaltranscribing;
    end
    NucPerTime = [NucPerTime;nucperbin];
    TranscribingPerTime = [TranscribingPerTime;transcribingperbin];
    [remove_nuc wid] = size(nucon);
    for n = 1:remove_nuc
        nuc_id = nucon(n,1);
        NucForProb(n,:) = [];
        DotsForProb(n,:) = [];
    end
end
                
%NucPerTime tells you how many nuclei could turn on a dot that haven't yet
%TranscibingPerTime tells you how many of the nuclei in NucPerTime do turn
%on a dot at that time point

sumnuc = sum(NucPerTime,1);
sumtranscribing = sum(TranscribingPerTime,1);
nucdenominator = [];
for s = 1:150
    denom = sum(sumtranscribing(s:end));
    nucdenominator = [nucdenominator denom];
end

ProbabilityMat = sumtranscribing(:,1:150)./nucdenominator;

figure
hold on
title("Probability of Promoter Firing")
ylabel("Probability On")
xlabel("Med Threshold")
ylim([0 1])
xData = (0.001:0.001:0.15);
yData = ProbabilityMat(1:150);
plot(xData,yData);
hold off

%% Rewrite MS2 data for comparison to fixed (distance binned)

TimeOn = AllTrackedNuclei(:,1);
Distance = AllTrackedNuclei(:,4);
DistanceBins = discretize(Distance(:,1),[0:3:60,200],"IncludedEdge",'right');
MedLevels = AllTrackedNuclei(:,5:end);

numbins=21;

DistBinnedMed = [];
DistBinnedProb = zeros(NumTimepoints,numbins);

for bin = 1:numbins
    medperbin = MedLevels(DistanceBins(:,1)==bin,:);
    avgmedperbin = mean(medperbin,1);
    DistBinnedMed = [DistBinnedMed,avgmedperbin'];
    totalnuc = length(medperbin);
    nucinbin = TimeOn(DistanceBins(:,1)==bin,:);

    for t = 1:35
        transcribingnuc = nucinbin(nucinbin(:,1)<=t,:);
        nucon = length(transcribingnuc);
        percenton = nucon/totalnuc;
        DistBinnedProb(t,bin) = percenton;
    end

end

figure
hold on
title(genename)
xlabel('log(MedGFP)')
ylabel('log(P/1-P)')
ylim([-4,6]);
xlim([-4,1]);
scatter(log(DistBinnedMed.*6),log(DistBinnedProb./(1-DistBinnedProb)),5,'MarkerFaceColor','blue','MarkerEdgeColor','blue')
hold off

