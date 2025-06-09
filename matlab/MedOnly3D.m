%% Set Image Colors

mymagenta = '#D81B60';
myblue = '#1E88E5';
mygreen = '#004D40';
myyellow = '#FFC107';

%FigurePath = strcat(Path,'Presentations/PlacedImages/');

set(0,'DefaultAxesFontName','Arial');

% load heatmap settings
mymagentamap = load('mygreenmap.mat').mygreenmap;
mybluemap = load('mybluemap.mat').mybluemap;
mygreenmap = load('mymagentamap.mat').mymagentamap;
myyellowmap = load('myyellowmap.mat').myyellowmap;

s.matlab.fonts.editor.code.Name.TemporaryValue = 'Arial';


%% Read in data
%Path = 'C:/Users/Susanna Brantley/Brantley Lab Dropbox/Susanna Brantley/Duke/DiTaliaLab/ImageAnalysis/Datasets/';
Path = 'D:/Brantley Lab Dropbox/Susanna Brantley/Duke/DiTaliaLab/ImageAnalysis/Datasets/';
%Path = '/Users/susannabrantley/Dropbox/Duke/DiTaliaLab/ImageAnalysis/Datasets/GFPMed_HisiRFP_sogmut/';

%Read in data file for Windows for reference embryo
%CellularizationFront = readmatrix(strcat(Path,"CellularizationFrontReference_live.csv"));
RefWT = readmatrix(strcat(Path,"GFPMed_HisiRFP/230727_wt_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));
RefHet = readmatrix(strcat(Path,"GFPMed_HisiRFP_sogmut/230727_het2_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));
D1_binned_wt = readmatrix(strcat(Path,'BinnedRefEmbryo_D1binned_wt_230730.csv'));
D1_binned_sog = readmatrix(strcat(Path,'BinnedRefEmbryo_D1binned_sog_230730_2.csv'));

%File Labels:
%ImageNo | ObjectNo | Area | Int | X | Y | Min | Gen | Z 

%Load all embryo data
%EmbryoMadMed = readmatrix("D:/Dropbox/Duke/DiTaliaLab/ImageAnalysis/Datasets/211107_wt_GFPMed_pMad_cellprofiler_cleaned.csv");

WTEmbryo1 = readmatrix(strcat(Path,"GFPMed_HisiRFP/230727_wt_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));
WTEmbryo2 = readmatrix(strcat(Path,"GFPMed_HisiRFP/230730_wt_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));
WTEmbryo3 = readmatrix(strcat(Path,"GFPMed_HisiRFP/221228_008_009_wt_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));
%WTEmbryo4 = readmatrix(strcat(Path,"GFPMed_HisiRFP/221228_021_wt_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));
WTEmbryo5 = readmatrix(strcat(Path,"GFPMed_HisiRFP/230104_001_wt_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned"));
%WTEmbryo6 = readmatrix(strcat(Path,"GFPMed_HisiRFP/230119_003_wt_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));
WTEmbryo7 = readmatrix(strcat(Path,"GFPMed_HisiRFP/230217_002_wt_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));
%WTEmbryo8 = readmatrix(strcat(Path,"GFPMed_HisiRFP/230218_003_wt_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));
WTEmbryo9 = readmatrix(strcat(Path,"240429_1_wt_GFPMed_ushMS2_cleaned.csv"));
WTEmbryo10 = readmatrix(strcat(Path,"240429_2_wt_GFPMed_ushMS2_cleaned.csv"));
WTEmbryo11 = readmatrix(strcat(Path,"240429_3_wt_GFPMed_ushMS2_cleaned.csv"));
WTEmbryo12 = readmatrix(strcat(Path,"240502_1_wt_GFPMed_ushMS2_cleaned.csv"));

HetEmbryo3 = readmatrix(strcat(Path,"GFPMed_HisiRFP_sogmut/230726_het_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));
HetEmbryo1 = readmatrix(strcat(Path,"GFPMed_HisiRFP_sogmut/230727_1_sog_GFPMed_HisiRFP_ByZ_cleaned.csv"));
HetEmbryo4 = readmatrix(strcat(Path,"GFPMed_HisiRFP_sogmut/230727_het_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));
HetEmbryo5 = readmatrix(strcat(Path,"GFPMed_HisiRFP_sogmut/230727_het2_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));

SogEmbryo6 = readmatrix(strcat(Path,"GFPMed_HisiRFP_sogmut/230729_sog_GFPMed_HisiRFP_ByZ_cellprofiler_cleaned.csv"));
%SogEmbryo7 = readmatrix(strcat(Path,"GFPMed_HisiRFP_sogmut/230727_1_sog_GFPMed_HisiRFP_ByZ_cleaned.csv"));
SogEmbryo9 = readmatrix(strcat(Path,"GFPMed_HisiRFP_sogmut/230730_1_sog_GFPMed_HisiRFP_ByZ_cleaned.csv"));
SogEmbryo10 = readmatrix(strcat(Path,"GFPMed_HisiRFP_sogmut/230730_2_sog_GFPMed_HisiRFP_ByZ_cleaned.csv"));

%ZenEmbryo1 = readmatrix(strcat(Path,"GFPMed_HisiRFP/250424_wt_GFPMed_HisiRFP_ncratio_cleaned.csv"));
%ZenEmbryo1 = readmatrix(strcat(Path,"GFPMed_HisiRFP/250501_wt_GFPMed_HisiRFP_ncratio_cleaned.csv"));
ZenEmbryo1 = readmatrix(strcat(Path,"250424_zen_GFPMed_HisiRFP_cleaned.csv"));
ZenEmbryo2 = readmatrix(strcat(Path,"250501_zen_GFPMed_HisiRFP_cleaned.csv"));

WTEmbryoData = {WTEmbryo1;WTEmbryo2;WTEmbryo3;WTEmbryo5;WTEmbryo7;WTEmbryo9(:,1:8);WTEmbryo10(:,1:8);WTEmbryo11(:,1:8);WTEmbryo12(:,1:8)};
HetEmbryoData = {HetEmbryo4;HetEmbryo5;HetEmbryo3};
SogEmbryoData ={SogEmbryo6;SogEmbryo9;SogEmbryo10};
ZenEmbryoData = {ZenEmbryo1,ZenEmbryo2};

%Set shared parameters
DistanceThreshold = 3.2;
TimeDistanceThreshold = 10;
PixelToMicron = 4.4;
EndTime = 60;
%PixelCutoffLow = 200; %300 for 008 | 300 for 021 | 350 for 001 | 350 for 002 | 275 for 003
%PixelCutoffHigh = 750; %800 for 008 | 750 for 021 | 725 for 001 | 750 for 002 | 700 for 003
Micron_cutoff = 75;
%NumBins = 9; %15 for 008 | 10 for 002 | 9 for 003

 %%
NumEmbryos = length(WTEmbryoData);
Bins = [15 15 15 15 15 15 15 15 15];
LowCutoff = [300 200 300 300 350 325 300 350 350];
HighCutoff = [700 600 900 850 750 750 700 725 750];
BackgroundLow = [55 55 55 55 55 55 55 55 55]; 
BackgroundHigh = [65 65 65 65 65 65 65 65 65];
StartTimes = [10,10,24,24,24,2,9,5,10];
ZSizes = [11,11,11,11,11,12,13,15,7];

% NumEmbryos = length(HetEmbryoData);
% Bins = [15 15 15];
% LowCutoff = [200 200 300];
% HighCutoff = [750 750 850];
% BackgroundLow = [55 55 55]; 
% BackgroundHigh = [65 65 65];
% StartTimes = [11,11,11];
% ZSizes = [11,11,11];

WidPerEmbryo = [];
AmplitudePerEmbryo = [];
AllTrackedNuclei = [];
WTScaleFactors = [];

for embryo = 1:NumEmbryos
    NumBins = Bins(1,embryo);
    NumZs = ZSizes(1,embryo);
    StartTime = StartTimes(1,embryo);
    NumTimepoints = EndTime-StartTime;
    PixelCutoffLow = LowCutoff(1,embryo);
    PixelCutoffHigh = HighCutoff(1,embryo);
    BackgroundCutoffLow = BackgroundLow(1,embryo);
    BackgroundCutoffHigh = BackgroundHigh(1,embryo);
    %NumTimepoints = NumTimepoints(1,embryo);
    MedData = WTEmbryoData{embryo};
    %Separate data by frame
    [~,~,X] = unique(MedData(:,1));
    frame_array = accumarray(X,1:size(MedData,1),[],@(r){MedData(r,:)});

    %set column indexes for tracking in z
    ObjectNum_index = 2;
    Area_index = 3;
    MedIntensity_index = 4;
    X_index = 5;
    Y_index = 6;

    NextNuc_index = 8;
    NextArea_index = 9;
    NextMed_index = 10;
    NextX_index = 11;
    NextY_index = 12;

    [StackedNucleusArray,StackedAreaArray,StackedMedArray,StackedXArray,StackedYArray] = stackbyz_nodots( ...
        frame_array,NumTimepoints,NumZs,DistanceThreshold, ...
        X_index,Y_index,ObjectNum_index,MedIntensity_index,Area_index, ...
        NextNuc_index,NextArea_index,NextMed_index,NextX_index,NextY_index);

    NucleiStacked = avgperstack_nodots( ...
        StackedNucleusArray,StackedAreaArray,StackedMedArray,StackedXArray,StackedYArray,...
        NumTimepoints);

    %set coordinates for tracking in time
    ObjectNum_index = 1;
    Area_index = 2;
    MedIntensity_index = 4;
    X_index = 5;
    Y_index = 6;

    NextNuc_index = 7;
    NextMed_index = 8;
    NextX_index = 9;
    NextY_index = 10;

    [NucleiOverTime,TrackedNucleusArray,TrackedMedArray,TrackedXArray,TrackedYArray] = tracknuclei_nodots( ...
        NucleiStacked,TimeDistanceThreshold,NumTimepoints,PixelCutoffLow,PixelCutoffHigh, ...
        X_index,Y_index,ObjectNum_index,MedIntensity_index, ...
        NextNuc_index,NextMed_index,NextX_index,NextY_index);

    [FullyTrackedNuclei,FullyTrackedX,FullyTrackedY] = tracktoend_nodots( ...
        EndTime,NumTimepoints, ...
        TrackedYArray,TrackedMedArray,TrackedXArray);

    MidlineTimepoint = 55;
    MicronY_index = 7;

    [CorrectedNuclei, FullyTrackedNucleiMidline] =  findmidline_nodots(FullyTrackedY,FullyTrackedNuclei, ...
    NumTimepoints,StartTime,MidlineTimepoint,PixelCutoffLow,PixelCutoffHigh,PixelToMicron,BackgroundCutoffLow,BackgroundCutoffHigh,...
    Y_index,MicronY_index,MedIntensity_index);
    
    FullyTrackedNucleiRescaled = CorrectedNuclei;
    %[scale_factor,FullyTrackedNucleiRescaled] = rescale(D1_binned_wt,NumTimepoints,CorrectedNuclei);
    %WTScaleFactors = [WTScaleFactors,scale_factor];
    Positions = FullyTrackedNucleiRescaled(:,1);
    MedLevels = FullyTrackedNucleiRescaled(:,2:end);
    PlaceHolders = zeros(length(FullyTrackedNucleiRescaled),StartTime);
    AllTrackedNuclei = [AllTrackedNuclei;Positions,PlaceHolders,MedLevels];

    [AmplitudeOverTime,WidthOverTime] = widthamp(NumTimepoints,FullyTrackedNucleiRescaled,11);
    WidthOverTime = [zeros(StartTime,2);WidthOverTime];
    AmplitudeOverTime = [zeros(StartTime,1);AmplitudeOverTime];
    WidPerEmbryo = [WidPerEmbryo WidthOverTime];
    AmplitudePerEmbryo = [AmplitudePerEmbryo AmplitudeOverTime];

    %AmplitudeOverDist = ampofdistbins(NumTimepoints,FullyTrackedNucleiRescaled,StartTime,EndTime,NumBins);

    %frame_array={};

    cmap=cbrewer2('Blues',70);
    figure
    hold on
    for mat = 1:5:NumTimepoints
        %figure
        xData = NucleiStacked{mat}(:,6);
        %xData = xData./PixelToMicron;
        yData = NucleiStacked{mat}(:,4);
        c = cmap(mat+10,:);
        sz = 10;
        scatter(xData,yData,sz,c);
        %ylim([0 .05]);
        %xlim([200 800]);
    end
    hold off

end


%% Analyze Sog data

Bins = [15 15 15 15 15];
LowCutoff = [150 175 175 175 175];
HighCutoff = [900 750 750 750 750];
BackgroundLow = [55 55 55 55 55]; 
BackgroundHigh = [65 65 65 65 65];
StartTimes = [11,11,11,11,11];
ZSizes = [11,11,11,11,11];

NumEmbryos = length(HetEmbryoData);
%StartTime = 10;

WidPerEmbryoSog = [];
AmplitudePerEmbryoSog = [];
AllTrackedNucleiSog = [];

%scale_factor = mean(WTScaleFactors);

for embryo = 1:NumEmbryos
    NumBins = Bins(1,embryo);
    StartTime = StartTimes(1,embryo);
    NumZs = ZSizes(1,embryo);
    NumTimepoints = EndTime-StartTime;
    PixelCutoffLow = LowCutoff(1,embryo);
    PixelCutoffHigh = HighCutoff(1,embryo);
    BackgroundCutoffLow = BackgroundLow(1,embryo);
    BackgroundCutoffHigh = BackgroundHigh(1,embryo);
    %NumTimepoints = NumTimepoints(1,embryo);
    MedData = HetEmbryoData{embryo};
    %Separate data by frame
    [~,~,X] = unique(MedData(:,1));
    frame_array = accumarray(X,1:size(MedData,1),[],@(r){MedData(r,:)});

    %set column indexes for tracking in z
    ObjectNum_index = 2;
    Area_index = 3;
    MedIntensity_index = 4;
    X_index = 5;
    Y_index = 6;

    NextNuc_index = 9;
    NextArea_index = 10;
    NextMed_index = 11;
    NextX_index = 12;
    NextY_index = 13;

    [StackedNucleusArray,StackedAreaArray,StackedMedArray,StackedXArray,StackedYArray] = stackbyz_nodots( ...
        frame_array,NumTimepoints,NumZs,DistanceThreshold, ...
        X_index,Y_index,ObjectNum_index,MedIntensity_index,Area_index, ...
        NextNuc_index,NextArea_index,NextMed_index,NextX_index,NextY_index);

    NucleiStacked = avgperstack_nodots( ...
        StackedNucleusArray,StackedAreaArray,StackedMedArray,StackedXArray,StackedYArray,...
        NumTimepoints);

    %set coordinates for tracking in time
    ObjectNum_index = 1;
    Area_index = 2;
    MedIntensity_index = 4;
    X_index = 5;
    Y_index = 6;

    NextNuc_index = 7;
    NextMed_index = 8;
    NextX_index = 9;
    NextY_index = 10;

    [NucleiOverTime,TrackedNucleusArray,TrackedMedArray,TrackedXArray,TrackedYArray] = tracknuclei_nodots( ...
        NucleiStacked,TimeDistanceThreshold,NumTimepoints,PixelCutoffLow,PixelCutoffHigh, ...
        X_index,Y_index,ObjectNum_index,MedIntensity_index, ...
        NextNuc_index,NextMed_index,NextX_index,NextY_index);

    [FullyTrackedNuclei,FullyTrackedX,FullyTrackedY] = tracktoend_nodots( ...
        EndTime,NumTimepoints, ...
        TrackedYArray,TrackedMedArray,TrackedXArray);

    MidlineTimepoint = 55;
    MicronY_index = 7;

    [CorrectedNuclei, FullyTrackedNucleiMidline] =  findmidline_nodots(FullyTrackedY,FullyTrackedNuclei, ...
    NumTimepoints,StartTime,MidlineTimepoint,PixelCutoffLow,PixelCutoffHigh,PixelToMicron,BackgroundCutoffLow,BackgroundCutoffHigh,...
    Y_index,MicronY_index,MedIntensity_index);

    [scale_factor,FullyTrackedNucleiRescaled] = rescale(D1_binned_sog,NumTimepoints,CorrectedNuclei);
    Positions = FullyTrackedNucleiRescaled(:,1);
    MedLevels = FullyTrackedNucleiRescaled(:,2:end);
    PlaceHolders = zeros(length(FullyTrackedNucleiRescaled),StartTime);
    AllTrackedNucleiSog = [AllTrackedNucleiSog;Positions,PlaceHolders,MedLevels];

    [AmplitudeOverTime,WidthOverTime] = widthamp(NumTimepoints,FullyTrackedNucleiRescaled,StartTime);
    WidthOverTime = [zeros(StartTime,2);WidthOverTime];
    AmplitudeOverTime = [zeros(StartTime,1);AmplitudeOverTime];
    WidPerEmbryoSog = [WidPerEmbryoSog WidthOverTime];
    AmplitudePerEmbryoSog = [AmplitudePerEmbryoSog AmplitudeOverTime];

    %AmplitudeOverDist = ampofdistbins(NumTimepoints,FullyTrackedNucleiRescaled,StartTime,EndTime,NumBins);

    %frame_array={};

    cmap=cbrewer2('Blues',70);
    figure
    hold on
    for mat = 1:5:NumTimepoints
        %figure
        xData = NucleiStacked{mat}(:,6);
        %xData = xData./PixelToMicron;
        yData = NucleiStacked{mat}(:,4);
        c = cmap(mat+10,:);
        sz = 10;
        scatter(xData,yData,sz,c);
        %plot(xData,yData,'Color',c)
        %ylim([0 .05]);
        %xlim([200 800]);
    end
    hold off

end

%% Zen data

NumEmbryos = length(ZenEmbryoData);
Bins = [15 15 15];
LowCutoff = [200 325 300];
HighCutoff = [750 750 850];
BackgroundLow = [55 55 55]; 
BackgroundHigh = [65 65 65];
ScaleFactors = [2.7,8.88];
StartTimes = [16,11];
ZSizes = [7,12];

ZenWidPerEmbryo = [];
ZenAmplitudePerEmbryo = [];
ZenAllTrackedNuclei = [];
WTScaleFactors = [];

for embryo = 1:NumEmbryos
    NumBins = Bins(1,embryo);
    NumZs = ZSizes(1,embryo);
    StartTime = StartTimes(1,embryo);
    ScaleFactor = ScaleFactors(1,embryo);
    NumTimepoints = EndTime-StartTime;
    PixelCutoffLow = LowCutoff(1,embryo);
    PixelCutoffHigh = HighCutoff(1,embryo);
    BackgroundCutoffLow = BackgroundLow(1,embryo);
    BackgroundCutoffHigh = BackgroundHigh(1,embryo);
    %NumTimepoints = NumTimepoints(1,embryo);
    MedData = ZenEmbryoData{embryo};
    %Separate data by frame
    [~,~,X] = unique(MedData(:,1));
    frame_array = accumarray(X,1:size(MedData,1),[],@(r){MedData(r,:)});

    %set column indexes for tracking in z
    ObjectNum_index = 2;
    Area_index = 3;
    MedIntensity_index = 4;
    X_index = 5;
    Y_index = 6;

    NextNuc_index = 8;
    NextArea_index = 9;
    NextMed_index = 10;
    NextX_index = 11;
    NextY_index = 12;

    [StackedNucleusArray,StackedAreaArray,StackedMedArray,StackedXArray,StackedYArray] = stackbyz_nodots( ...
        frame_array,NumTimepoints,NumZs,DistanceThreshold, ...
        X_index,Y_index,ObjectNum_index,MedIntensity_index,Area_index, ...
        NextNuc_index,NextArea_index,NextMed_index,NextX_index,NextY_index);

    NucleiStacked = avgperstack_nodots( ...
        StackedNucleusArray,StackedAreaArray,StackedMedArray,StackedXArray,StackedYArray,...
        NumTimepoints);

    %set coordinates for tracking in time
    ObjectNum_index = 1;
    Area_index = 2;
    MedIntensity_index = 4;
    X_index = 5;
    Y_index = 6;

    NextNuc_index = 7;
    NextMed_index = 8;
    NextX_index = 9;
    NextY_index = 10;

    [NucleiOverTime,TrackedNucleusArray,TrackedMedArray,TrackedXArray,TrackedYArray] = tracknuclei_nodots( ...
        NucleiStacked,TimeDistanceThreshold,NumTimepoints,PixelCutoffLow,PixelCutoffHigh, ...
        X_index,Y_index,ObjectNum_index,MedIntensity_index, ...
        NextNuc_index,NextMed_index,NextX_index,NextY_index);

    [FullyTrackedNuclei,FullyTrackedX,FullyTrackedY] = tracktoend_nodots( ...
        EndTime,NumTimepoints, ...
        TrackedYArray,TrackedMedArray,TrackedXArray);

    MidlineTimepoint = 55;
    MicronY_index = 7;

    [CorrectedNuclei, FullyTrackedNucleiMidline] =  findmidline_nodots(FullyTrackedY,FullyTrackedNuclei, ...
    NumTimepoints,StartTime,MidlineTimepoint,PixelCutoffLow,PixelCutoffHigh,PixelToMicron,BackgroundCutoffLow,BackgroundCutoffHigh,...
    Y_index,MicronY_index,MedIntensity_index);
    
    FullyTrackedNucleiRescaled = CorrectedNuclei;
    %[scale_factor,FullyTrackedNucleiRescaled] = rescale(D1_binned_wt,NumTimepoints,CorrectedNuclei);
    %WTScaleFactors = [WTScaleFactors,scale_factor];
    Positions = FullyTrackedNucleiRescaled(:,1);
    MedLevels = FullyTrackedNucleiRescaled(:,2:end).*ScaleFactor;
    PlaceHolders = zeros(length(FullyTrackedNucleiRescaled),StartTime);
    ZenAllTrackedNuclei = [ZenAllTrackedNuclei;Positions,PlaceHolders,MedLevels];

    [ZenAmplitudeOverTime,ZenWidthOverTime] = widthamp(NumTimepoints,FullyTrackedNucleiRescaled,11);
    ZenWidthOverTime = [zeros(StartTime,2);ZenWidthOverTime];
    ZenAmplitudeOverTime = [zeros(StartTime,1);ZenAmplitudeOverTime];
    ZenWidPerEmbryo = [ZenWidPerEmbryo ZenWidthOverTime];
    ZenAmplitudePerEmbryo = [ZenAmplitudePerEmbryo ZenAmplitudeOverTime];

    %AmplitudeOverDist = ampofdistbins(NumTimepoints,FullyTrackedNucleiRescaled,StartTime,EndTime,NumBins);

    %frame_array={};

    cmap=cbrewer2('Blues',70);
    figure
    hold on
    for mat = 1:5:NumTimepoints
        %figure
        xData = NucleiStacked{mat}(:,6);
        %xData = xData./PixelToMicron;
        yData = NucleiStacked{mat}(:,4);
        c = cmap(mat+10,:);
        sz = 10;
        scatter(xData,yData,sz,c);
        %ylim([0 .05]);
        %xlim([200 800]);
    end
    hold off

end


%% Plot avg width and amplitude
% Should convert these plots to fits

AvgWTAmp = [];
AvgWTAmpDev = [];
for x = 1:40
    ampatt = nonzeros(AmplitudePerEmbryo(x,:));
    avgampatt = mean(ampatt);
    stdatt = std(ampatt)./sqrt(length(ampatt));
    AvgWTAmp = [AvgWTAmp;avgampatt];
    AvgWTAmpDev = [AvgWTAmpDev;stdatt];
end

AvgSogAmp = [];
AvgSogAmpDev = [];
for x = 1:40
    ampatt = nonzeros(AmplitudePerEmbryoSog(x,:));
    avgampatt = mean(ampatt);
    stdatt = std(ampatt)./sqrt(length(ampatt));
    AvgSogAmp = [AvgSogAmp;avgampatt];
    AvgSogAmpDev = [AvgSogAmpDev;stdatt];
end

AvgWTAmpDev = [AvgWTAmpDev,AvgWTAmp-AvgWTAmpDev,AvgWTAmp+AvgWTAmpDev];
AvgSogAmpDev = [AvgSogAmpDev,AvgSogAmp-AvgSogAmpDev,AvgSogAmp+AvgSogAmpDev];

figure
hold on
fontsize(14,"points")
select = ~isnan(AvgWTAmp);
x = [21:60]';
x = x(select);
WTAmp = AvgWTAmp(select);
WTDevUp = AvgWTAmpDev(:,2);
WTDevUp = WTDevUp(select);
WTDevDown = AvgWTAmpDev(:,3);
WTDevDown = WTDevDown(select);
fill([x', flip(x')], [WTDevUp', flip(WTDevDown')], [0.12, 0.53, 0.9],'FaceAlpha',0.2, 'EdgeColor','none')
plot(x,WTAmp,'LineWidth',2,'Color',myblue)

select = ~isnan(AvgSogAmp);
x = [21:60]';
x = x(select);
SogAmp = AvgSogAmp(select);
SogDevUp = AvgSogAmpDev(:,2);
SogDevUp = SogDevUp(select);
SogDevDown = AvgSogAmpDev(:,3);
SogDevDown = SogDevDown(select);
fill([x', flip(x')], [SogDevUp', flip(SogDevDown')], [0.85, 0.11, 0.38],'FaceAlpha',0.2, 'EdgeColor','none')
plot(x,SogAmp,'LineWidth',2,'Color',mymagenta)
%xlim([20,60])
xlabel("Minutes")
ylabel("Amplitude")
hold off

WidthPerEmbryo = WidPerEmbryo(:,2:2:6);
WidthPerEmbryoSog = WidPerEmbryoSog(:,2:2:6);

AvgWTWidth = [];
AvgWTWidthDev = [];
for x = 1:40
    att = nonzeros(WidthPerEmbryo(x,:));
    avgatt = mean(att);
    stdatt = std(att)./sqrt(length(att));
    AvgWTWidth = [AvgWTWidth;avgatt];
    AvgWTWidthDev = [AvgWTWidthDev;stdatt];
end

AvgSogWidth = [];
AvgSogWidthDev = [];
for x = 1:40
    att = nonzeros(WidthPerEmbryoSog(x,:));
    avgatt = mean(att);
    stdatt = std(att)./sqrt(length(att));
    AvgSogWidth = [AvgSogWidth;avgatt];
    AvgSogWidthDev = [AvgSogWidthDev;stdatt];
end

AvgWTWidthDev = [AvgWTWidthDev,AvgWTWidth-AvgWTWidthDev,AvgWTWidth+AvgWTWidthDev];
AvgSogWidthDev = [AvgSogWidthDev,AvgSogWidth-AvgSogWidthDev,AvgSogWidth+AvgSogWidthDev];

figure
hold on
fontsize(14,"points")
select = ~isnan(AvgWTWidth);
x = [21:60]';
x = x(select);
WTWidth = AvgWTWidth(select);
WTDevUp = AvgWTWidthDev(:,2);
WTDevUp = WTDevUp(select);
WTDevDown = AvgWTWidthDev(:,3);
WTDevDown = WTDevDown(select);
fill([x', flip(x')], [WTDevUp', flip(WTDevDown')], [0.12, 0.53, 0.9],'FaceAlpha',0.2, 'EdgeColor','none')
plot(x,WTWidth,'LineWidth',2,'Color',myblue)

select = ~isnan(AvgSogWidth);
x = [21:60]';
x = x(select);
SogWidth = AvgSogWidth(select);
SogDevUp = AvgSogWidthDev(:,2);
SogDevUp = SogDevUp(select);
SogDevDown = AvgSogWidthDev(:,3);
SogDevDown = SogDevDown(select);
fill([x', flip(x')], [SogDevUp', flip(SogDevDown')], [0.85, 0.11, 0.38],'FaceAlpha',0.2, 'EdgeColor','none')
plot(x,SogWidth,'LineWidth',2,'Color',mymagenta)
%scatter(x,SogWidth)
%xlim([20,60])
xlabel("Minutes")
ylabel("Width (microns)")
hold off



%% Bin by distance and time and plot fitted curves

binedges = 0:5:75;
distbinnednuclei = discretize(AllTrackedNuclei(:,1),binedges,'IncludedEdge','left');
figure
hold on
fontsize(14,'points')
cmap = cbrewer2('Blues',length(binedges)+10);
for bin = 1:length(binedges)-1
    nucinbin = AllTrackedNuclei(distbinnednuclei(:,1)==bin,2:end);
    medovertime = mean(nucinbin,1);
    semovertime = std(nucinbin,1)./sqrt(length(nucinbin));
    c = cmap(length(binedges)+10-bin,:);
    x = 1:60;
    fill([x, flip(x)], [medovertime+semovertime, flip(medovertime-semovertime)], c,'FaceAlpha',0.5, 'EdgeColor','none')
    plot(x,medovertime,'LineWidth',2,'Color',c)
end
xlim([20,60])
ylim([0,0.05])
xlabel('Time (Minutes)')
ylabel('Nuclear GFP-Med')
hold off

binedges = 0:5:75;
distbinnednuclei = discretize(ZenAllTrackedNuclei(:,1),binedges,'IncludedEdge','left');
figure
hold on
fontsize(14,'points')
cmap = cbrewer2('Greens',length(binedges)+10);
for bin = 1:length(binedges)-1
    nucinbin = ZenAllTrackedNuclei(distbinnednuclei(:,1)==bin,2:end);
    medovertime = mean(nucinbin,1);
    semovertime = std(nucinbin,1)./sqrt(length(nucinbin));
    c = cmap(length(binedges)+10-bin,:);
    x = 1:60;
    fill([x, flip(x)], [medovertime+semovertime, flip(medovertime-semovertime)], c,'FaceAlpha',0.5, 'EdgeColor','none')
    plot(x,medovertime,'LineWidth',2,'Color',c)
end
xlim([20,60])
%ylim([0,0.05])
xlabel('Time (Minutes)')
ylabel('Nuclear GFP-Med')
hold off

%%
binedges = 0:5:75;
distbinnednuclei = discretize(AllTrackedNucleiSog(:,1),binedges,'IncludedEdge','left');
MedOverTime = [];
figure
hold on
fontsize(14,'points')
cmap = cbrewer2('Oranges',length(binedges)+10);
for bin = 1:length(binedges)-1
    nucinbin = AllTrackedNucleiSog(distbinnednuclei(:,1)==bin,2:end);
    medovertime = mean(nucinbin,1);
    MedOverTime = [MedOverTime,medovertime'];
    semovertime = std(nucinbin,1)./sqrt(length(nucinbin));
    c = cmap(length(binedges)+10-bin,:);
    x = 1:60;
    fill([x, flip(x)], [medovertime+semovertime, flip(medovertime-semovertime)], c,'FaceAlpha',0.5, 'EdgeColor','none')
    plot(x,medovertime,'LineWidth',2,'Color',c)
end
xlim([20,60])
ylim([0,0.04])
xlabel('Time (Minutes)')
ylabel('Nuclear GFP-Med')
hold off

figure
heatmap(MedOverTime)

%% Find the derivative over time

binedges = 0:5:75;
distbinnednuclei = discretize(AllTrackedNuclei(:,1),binedges,'IncludedEdge','left');
figure
hold on
fontsize(14,'points')
cmap = cbrewer2('Blues',length(binedges)+10);
for bin = 1:length(binedges)-1
    nucinbin = AllTrackedNuclei(distbinnednuclei(:,1)==bin,2:end);
    medovertime = mean(nucinbin,1);
    ft = fit((1:60)',medovertime','SmoothingSpline','SmoothingParam',0.005);
    d1 = differentiate(ft,20:60);
    c = cmap(length(binedges)+10-bin,:);
    x = 20:60;
    plot(x,d1,'LineWidth',2,'Color',c,'LineStyle','--')
end
xlim([20,60])
xlabel('Time (Minutes)')
ylabel('dMed/dt')
hold off

binedges = 0:5:75;
distbinnednuclei = discretize(AllTrackedNucleiSog(:,1),binedges,'IncludedEdge','left');
figure
hold on
fontsize(14,'points')
cmap = cbrewer2('PuRd',length(binedges)+10);
for bin = 1:length(binedges)-1
    nucinbin = AllTrackedNucleiSog(distbinnednuclei(:,1)==bin,2:end);
    medovertime = mean(nucinbin,1);
    ft = fit((1:60)',medovertime','SmoothingSpline','SmoothingParam',0.005);
    d1 = differentiate(ft,20:60);
    c = cmap(length(binedges)+10-bin,:);
    x = 20:60;
    plot(x,d1,'LineWidth',2,'Color',c,'LineStyle','--')
end
xlim([20,60])
xlabel('Time (Minutes)')
ylabel('dMed/dt')
hold off
%% Find x vlaue per threshold

AllTrackedNuclei = sortrows(AllTrackedNuclei,1,'ascend');
widthamp(NumTimepoints,AllTrackedNuclei,11);

AllFittedMedSog = AllTrackedNuclei(:,1);

cmap=cbrewer2('Blues',60);
figure
hold on
for t = 1:NumTimepoints;
    xData = AllTrackedNuclei(:,1);
    yData = AllTrackedNuclei(:,t+1);

    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.0000005;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    yfitted = feval(fitresult,xData);
    AllFittedMed = [AllFittedMed yfitted];
    shade = 10+t;
    plot(xData,yfitted,'-','color',cmap(shade,:))
end
hold off

%%
AllFittedMed = sortrows(AllFittedMed,1,'descend');
MaxMed = [];

for row = 1:length(AllFittedMed);
    max_at_dist = max(AllFittedMed(row,2:end));
    MaxMed =[MaxMed;max_at_dist];
end

MaxThresholdDist = [];
MaxThresholdTime = [];
for threshold = 0.01:0.005:0.065
    threshold_dist = [threshold];
    threshold_time = [threshold];
    [l w] = size(MaxMed);
    for p = 1:l
        med_level = MaxMed(p,1);
        x_dist = AllFittedMed(p,1);
        if med_level>=threshold
            threshold_dist = [threshold_dist x_dist];
            threshold_time = [threshold_time t];
            break
        else
            continue
        end
    end
    dist_at_threshold = threshold_dist(1,2);
    time_at_threshold = threshold_time(1,2);
    MaxThresholdDist = [MaxThresholdDist; dist_at_threshold];
    MaxThresholdTime = [MaxThresholdTime; time_at_threshold];
end

OverTimeThresholdDist = [];
OverTimeThresholdTime = [];

for time = 5:5:NumTimepoints
    ThresholdDist = [];
    ThresholdTime = [];
    for threshold = 0.01:0.005:0.065
        threshold_dist = [threshold];
        threshold_time = [threshold];
        for t = time:NumTimepoints
            fitted_data = AllFittedMed(:,t+1);
            [l w] = size(fitted_data);
            for p = 1:l
                med_level = fitted_data(p,1);
                x_dist = AllFittedMed(p,1);
                if med_level>=threshold
                    threshold_dist = [threshold_dist x_dist];
                    threshold_time = [threshold_time t];
                    break
                else
                    continue
                end
            end
        end
        dist_at_threshold = threshold_dist(1,2);
        time_at_threshold = threshold_time(1,2);
        ThresholdDist = [ThresholdDist; dist_at_threshold];
        ThresholdTime = [ThresholdTime; time_at_threshold];
    end
    OverTimeThresholdDist = [OverTimeThresholdDist ThresholdDist];
    OverTimeThresholdTime = [OverTimeThresholdTime ThresholdTime];
end

OverTimeThresholdDist = [OverTimeThresholdDist MaxThresholdDist];
OverTimeThresholdTime = [OverTimeThresholdTime MaxThresholdTime];

cmap=cbrewer2('Blues',10);
figure
title("Max Distance v. Threshold")
xlabel("Threshold")
ylabel("Max Distance (um)")
hold on
for condition = 1:7;
    xData = 0.01:0.005:0.065;
    yData = OverTimeThresholdDist(:,condition);
    %c = cmap(OverTimeThresholdTime(:,condition));
    sz = 10;
    shade = 2+condition;
    plot(xData,yData,'-','color',cmap(shade,:));
    scatter(xData,yData,sz);
    ylim([0 40]);
    %xlim([0 0.07]);
end
hold off

%% Plot fitted curves
cmap=cbrewer2('Blues',55);
    figure
    hold on
for t = 1:5:NumTimepoints;
    xData = FullyTrackedNuclei(:,1);
    yData = FullyTrackedNuclei(:,t+1);
    xData = xData(yData ~= 0);    
    yData = yData(yData ~= 0);

    ft = fittype( 'exp2' );
    %opts = fitoptions( 'Method', 'SmoothingSpline');
    %opts.SmoothingParam = 0.005;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    yfitted = feval(fitresult,xData);
    c = cmap(t+10,:);
    sz = 10;
    shade = t+10;
    scatter(xData,yData,sz,c);
    plot(xData,yfitted,'-','color',cmap(shade,:));

end
hold off
%% Plot Width and Amplitude 
cmap=cbrewer2('Blues',20);
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
%xlim([25 60])
%ylim([0 50])
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
end
%xlim([25 60]);
%ylim([0 50])
hold off 
%% Check X cutoff

cmap=cbrewer2('Blues',56);
figure
hold on
for mat = 1:5:NumTimepoints
    %figure
    xData = NucleiStacked{mat}(:,6);
    %xData = xData./PixelToMicron;
    yData = NucleiStacked{mat}(:,4);
    c = cmap(mat+10,:);
    sz = 10;
    scatter(xData,yData,sz,c);
    ylim([0 .07]);
    %xlim([200 800]);
end
hold off

%% Fit curve to avgMed at each timepoint to find midline
%plot rescaled curve

AmplitudeOverTime = [];

figure
hold on
title('Med v. Distance (color = time)');
xlabel('Distance (um)');
ylabel('Avg Med Concentration')
for t = 1:NumTimepoints;
    distanceData = CorrectedNuclei{t}(:,MicronY_index);
    medData = CorrectedNuclei{t}(:,CorrectedMed_index);
    dataset = [distanceData medData];
    dataset = sortrows(dataset,1);
    filterbyy = dataset(dataset(:,1)<Micron_cutoff,:);
    yData = filterbyy(:,2);
    xData = filterbyy(:,1);
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.0005;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    yfitted = feval(fitresult,xData);
    fittedmed = [yfitted xData];
    amplitude = max(yfitted);
    AmplitudeOverTime = [AmplitudeOverTime; amplitude];
    shade = t+10;
    plot(xData,yfitted,'-','color',cmap(shade,:));
end
ylim([0 0.15]);
xlim([0 75]);
hold off

% Time = (EndTime+1-NumTimepoints:EndTime);
% Amplitude = AmplitudeOverTime(:,1);
% 
% figure
% plot(Time,Amplitude)
% title('Amplitude Over Time');
% xlabel('Time (min)');
% ylabel('Amplitude')
% ylim([0 0.15]);

%% Track nuclei over time

%plot NucleiOverTime 
figure
hold on
for t = 15;
    nuclei = NucleiOverTime{t};
    nuclei = sortrows(nuclei,MicronY_index);
    xData = nuclei(:,MicronY_index);
    yData = nuclei(:,CorrectedMed_index);
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.00005;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    yfitted = feval(fitresult,xData);
    fittedmed = yfitted;
    plot(xData,yData);
end
hold off

%% Order nuclei in a time course


%% Find nuclei tracked through entire movie

cmap=cbrewer2('Blues',NumTimepoints+20);
figure
hold on
title('Med v. Distance (color = time)');
xlabel('Distance (um)');
ylabel('Avg Med Concentration')
for t = 1:NumTimepoints;
    xData = FullyTrackedNuclei(:,1);
    yData = FullyTrackedNuclei(:,t+1);
    xData = xData(yData ~= 0);    
    yData = yData(yData ~= 0);
    %plot(xData, yData)
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.0005;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    yfitted = feval(fitresult,xData);
    fittedmed = yfitted;
    AmplitudeOverTime = [AmplitudeOverTime; max(yfitted)];
    shade = t+10;
    plot(xData,fittedmed,'-','color',cmap(shade,:));
    xlim([0 75]);
    ylim([0 0.15]);
end
hold off

%% Plot rescaled curves

AmplitudeOverTime = [];
WidthOverTime = [];

cmap=cbrewer2('Blues',NumTimepoints+20);
figure
hold on
title('Med v. Distance (color = time)');
xlabel('Distance (um)');
ylabel('Avg Med Concentration')
for t = 1:NumTimepoints;
    xData = FullyTrackedNucleiRescaled(:,1);
    yData = FullyTrackedNucleiRescaled(:,t+1);
    xData = xData(yData ~= 0);    
    yData = yData(yData ~= 0);

    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.0005;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    yfitted = feval(fitresult,xData);
    AmplitudeOverTime = [AmplitudeOverTime; max(yfitted)];
    halfMax = (min(yfitted) + max(yfitted)) / 2;
    max40 = 0.4*(min(yfitted) + max(yfitted));
    max60 = 0.6*(min(yfitted) + max(yfitted));
    x1 = find(yfitted >= halfMax, 1, 'first');
    x2 = find(yfitted >= halfMax, 1, 'last');
    x3 = find(yfitted >= max40, 1, 'first');
    x4 = find(yfitted >= max40, 1, 'last');
    x5 = find(yfitted >= max60, 1, 'first');
    x6 = find(yfitted >= max60, 1, 'last');
    width50 = xData(x2,:)-xData(x1,:);
    width40 = xData(x4,:)-xData(x3,:);
    width60 = xData(x6,:)-xData(x5,:);
    avgwidth = mean([width50 width40 width60]);
    width = trapz(xData,yData)./max(yfitted);
    WidthOverTime = [WidthOverTime; avgwidth width];

    shade = t+10;
    plot(xData,yfitted,'-','color',cmap(shade,:));
    xlim([0 75]);
    ylim([0 0.15]);
end
hold off

Time = (StartTime:EndTime);
Amplitude = AmplitudeOverTime(:,1);

figure
plot(Time,Amplitude)
title('Amplitude Over Time');
xlabel('Time (min)');
ylabel('Amplitude')
ylim([0 0.15]);

Time = (StartTime:EndTime);
Width = WidthOverTime(:,1);

figure
hold on
for width = 1:2;
    Width = WidthOverTime(:,width);
    plot(Time,Width)
end
ylim([10 50]);
title('Width Over Time');
xlabel('Time (min)');
ylabel('Width (um)');
hold off


%% Plot rescaled by amplitude

GaussianFitFeatures = [];

cmap=cbrewer2('Blues',NumTimepoints+20);
figure
hold on
title('Med Rescale By Amplitude (color = time)');
xlabel('Distance (um)');
ylabel('Avg Med Concentration')
for t = 1:NumTimepoints;
    xData = FullyTrackedNucleiRescaled(:,1);
    yData = FullyTrackedNucleiRescaled(:,t+1);
    xData = xData(yData ~= 0);    
    yData = yData(yData ~= 0);

    ft = fit(xData, yData, 'gauss1');
    yfitted = feval(fitresult,xData);
    integral = trapz(xData,yData);
    sigma = integral/(ft.a1*sqrt(2*pi));
    amp = ft(1);

    GaussianFitFeatures = [GaussianFitFeatures ; sigma amp];

    shade = t+10;
    plot(xData,yfitted./amp,'-','color',cmap(shade,:));
    %xlim([0 75]);
    %ylim([0 1.1]);
end
hold off

Time = (StartTime:EndTime);
RescaledWidth = GaussianFitFeatures(:,1);

figure
plot(Time,RescaledWidth)
%ylim([15 50]);
title('Rescaled Width Over Time');
xlabel('Time (min)');
ylabel('Width (um)');

%% Determine when each nucleus first turns on MS2

first_transcript = [];
[num_tracked_nuclei timepoints_tracked] = size(FullyTrackedNucleiRescaled);

for xx = 1:num_tracked_nuclei %number of nuclei that have been tracked
    counter = 4;
    dist_from_mid = FullyTrackedNucleiRescaled(xx,1);
    %nuc_index = FullyTrackedNucleiRescaled(xx,2);
    for yy = 5:NumTimepoints
        Med_level = FullyTrackedNucleiRescaled(xx,yy+1);
        if FullyTrackedDots(xx,yy) == 0
            counter = counter + 1;
            if counter == NumTimepoints
                checkMed = sum(FullyTrackedNucleiRescaled(xx,2:counter));
                if checkMed > 0
                    integrateMed = cumtrapz(FullyTrackedNucleiRescaled(xx,2:counter));
                    totalMed = integrateMed(end);
                    %collect_mat = [dist_index counter Med_level totalMed];
                    collect_mat = [counter Med_level totalMed xx];
                    first_transcript = [first_transcript; collect_mat];
                else
                    totalMed = 0;
                    collect_mat = [counter Med_level totalMed xx];
                    first_transcript = [first_transcript; collect_mat];
                end
                %break
            else
                continue
            end
        else
            counter = counter + 1;
            checkMed = sum(FullyTrackedNucleiRescaled(xx,2:counter));
            if checkMed > 0
                integrateMed = cumtrapz(FullyTrackedNucleiRescaled(xx,2:counter));
                totalMed = integrateMed(end);
                collect_mat = [counter Med_level totalMed xx];
                first_transcript = [first_transcript; collect_mat];
            else
                totalMed = 0;
                collect_mat = [counter Med_level totalMed xx];
                first_transcript = [first_transcript; collect_mat];
            end
            break
        end
    end
end

%first_transcript = sortrows(first_transcript);

%Remove nuclei that don't turn on dots
%first_transcript = first_transcript(1:67,:);

MS2_ontime = (first_transcript(:,1))+StartTime-1;
Transcribing_Nuclei = first_transcript(first_transcript(:,1)<31,:);
Med_level_at_on = Transcribing_Nuclei(:,2);
%Med_integral_at_on = first_transcript(:,3);

figure
hold on
title("Med Threshold at Promoter Firing in Transcribing Nuclei")
histogram(Med_level_at_on,'BinEdges', [0 0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.055 0.06 0.065 0.07 0.075 0.08 0.085 0.09 0.095 0.1])
hold off
% figure
% histogram(Med_integral_at_on,'BinLimits',[0.00001,0.01],'BinWidth',0.0005)
figure
hold on
title("Time To First Promoter Firing")
histogram(MS2_ontime)
hold off

%% Plot Med levels of nuclei that turn on MS2 and those that don't

cmap=cbrewer2('Blues',9);

[tracked_nuc tracked_timepoints] = size(FullyTrackedNucleiRescaled);

figure
hold on
title("Med v. Time (color = transcribing or not)")
for nuc = 1:tracked_nuc;
    xData = [StartTime:EndTime];
    yData = FullyTrackedNucleiRescaled(nuc,2:end);
    if first_transcript(nuc,1) == NumTimepoints
        binshade = 7;
        plot(xData,yData,'-','color',cmap(binshade,:));
        ylim([0 0.15]);
    else
        binshade = 4;
        plot(xData,yData,'-','color',cmap(binshade,:));
        ylim([0 0.15]);
    end
end
hold off

%% Determine if threshold is the same across time

BinnedOnTime = discretize(Transcribing_Nuclei(:,1)+24, [25 30 35 40 45 50 55]);
MedBinnedByTime = [BinnedOnTime Transcribing_Nuclei(:,1) Med_level_at_on];
[~,~,X] = unique(MedBinnedByTime(:,1));

[num_total_nuclei width_total_nuc] = size(MedBinnedByTime);

BinnedMedValues = {};
for bin = 2:4
    figure
    hold on
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

%% pMad v Med comparison

Minutes = [13 21 38 40.5 41 42.5 43.5 46 49.9 50.1 51 57.5];
Y_index = 8;
MedIntensity_index = 4;
MadIntensity_index = 6;
MidlineOverTime = [];

figure
hold on
for t = 1:NumTimepoints
    timepoint_nuclei = NucleiStacked{t};
    timepoint_nuclei = sortrows(timepoint_nuclei,Y_index);
    xData = timepoint_nuclei(:,Y_index);
    yData = timepoint_nuclei(:,MadIntensity_index);
    
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.000005;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    yfitted = feval(fitresult,xData);
    plot(xData,yfitted)
    fittedmed = [yfitted xData];
    amplitude = max(yfitted);
    midlineinfo = fittedmed(fittedmed(:,1)==amplitude,:);
    MidlineOverTime = [MidlineOverTime; midlineinfo];
    distfrommid = abs(xData-midlineinfo);
    NucleiStacked{t} = [timepoint_nuclei distfrommid./4.4];
end
hold off

Dist_index = 10;

cmap=cbrewer2('Blues',9);
figure
hold on
for t = 1:NumTimepoints
    timepoint_nuclei = NucleiStacked{t};
    timepoint_nuclei = sortrows(timepoint_nuclei,Dist_index);
    xData = timepoint_nuclei(:,Dist_index);
    yData = timepoint_nuclei(:,MadIntensity_index);

    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.00000005;
    [fitresult, gof] = fit( xData, yData, ft, opts );
    yfitted = feval(fitresult,xData);
    plot(xData,yfitted)
end
hold off

cmap=cbrewer2('Blues',20);
figure
hold on
xlabel('GFPMed Concentration');
ylabel('pMad Concentration');
for t = 1:NumTimepoints
    timepoint_nuclei = NucleiStacked{t};
    timepoint_nuclei = sortrows(timepoint_nuclei,Dist_index);
    xData = timepoint_nuclei(:,MedIntensity_index);
    yData = timepoint_nuclei(:,MadIntensity_index);
    c = cmap(t+10,:);
    sz = 10;
    shade = t+5;
    scatter(xData,yData,sz,c);
    %binshade = 5+t;
    %plot(xData,yData,'-','color',cmap(binshade,:));
end
hold off
