Path = 'C:/Users/Susanna Brantley/Brantley Lab Dropbox/Susanna Brantley/Duke/DiTaliaLab/ImageAnalysis/Datasets/';
%Path = 'D:/Dropbox/Duke/DiTaliaLab/ImageAnalysis/Datasets/';
%Path = '/Users/susannabrantley/Dropbox/Duke/DiTaliaLab/ImageAnalysis/Datasets/';

ExpFolder = 'ZenLT_BcdGFP_HisRFP/';

s.matlab.fonts.editor.code.Name.TemporaryValue = 'Arial';
set(0,'DefaultAxesFontName','Arial');

%Read in data file for Windows for reference embryo
%D1_binned = readmatrix(strcat(Path,"BinnedRefEmbryo_D1binned_230217_002.csv"));

%File Labels:
%ImageNo | ObjectNo | Area | Int | Mean | X | Y | Gen | Z 

%Load all embryo data
Nuclei1 = readmatrix(strcat(Path,ExpFolder,"240621_nc14_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Nuclei.csv"));
Nuclei2 = readmatrix(strcat(Path,ExpFolder,"240628_nc14_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Nuclei.csv"));
Nuclei3 = readmatrix(strcat(Path,ExpFolder,"240628_2_nc14_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Nuclei.csv"));

Cyto1 = readmatrix(strcat(Path,ExpFolder,"240621_nc14_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Cytoplasm.csv"));
Cyto2 = readmatrix(strcat(Path,ExpFolder,"240628_nc14_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Cytoplasm.csv"));
Cyto3 = readmatrix(strcat(Path,ExpFolder,"240628_2_nc14_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Cytoplasm.csv"));

NC13Nuc1 = readmatrix(strcat(Path,ExpFolder,"240628_nc13_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Nuclei.csv"));
NC13Cyto1 = readmatrix(strcat(Path,ExpFolder,"240628_nc13_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Cytoplasm.csv"));
NC13Nuc2 = readmatrix(strcat(Path,ExpFolder,"240628_2_nc13_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Nuclei.csv"));
NC13Cyto2 = readmatrix(strcat(Path,ExpFolder,"240628_2_nc13_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Cytoplasm.csv"));

NC12Nuc1 = readmatrix(strcat(Path,ExpFolder,"240628_nc12_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Nuclei.csv"));
NC12Cyto1 = readmatrix(strcat(Path,ExpFolder,"240628_nc12_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Cytoplasm.csv"));

%NC11Nuc1 = readmatrix(strcat(Path,ExpFolder,"240628_nc11_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Nuclei.csv"));
%NC11Cyto1 = readmatrix(strcat(Path,ExpFolder,"240628_nc11_zenLT_bcdGFPhet_HisiRFP_cellprofiler_Cytoplasm.csv"));

NC14Data = {[Nuclei1,Cyto1(:,3)],[Nuclei2,Cyto2(:,3)],[Nuclei3,Cyto3(:,3)]};
NC13Data = {[NC13Nuc1,NC13Cyto1(:,3)],[NC13Nuc2,NC13Cyto2(:,3)]};
NC12Data = {[NC12Nuc1,NC12Cyto1(:,3)]};
%NC11Data = {[NC11Nuc1,NC11Cyto1(:,3)]};

FullData = {[NC12Nuc1,NC12Cyto1(:,3)],[NC13Nuc1,NC13Cyto1(:,3)],[Nuclei2,Cyto2(:,3)]};

NumEmbryos = length(FullData);

%Set shared parameters
PixelCutoffLow = 0;
PixelCutoffHigh = 1024;
BackgroundCutoffLow = 75;
BackgroundCutoffHigh = 300;
DistanceThreshold = 3.2;
TimeDistanceThreshold = 50;
PixelToMicron = 4.4;
Micron_cutoff = 75;
EndTime = 60;
NumBins = 15;

%% Identify dataset-sepcific parameters

ZSlices = [];
Timepoints = [];
Start = [];

for e = 1:NumEmbryos
    eData = FullData{e};
    numz = 1;%max(eData(:,8));
    numt = max(eData(:,1))/numz;
    start = EndTime-numt+1;
    ZSlices = [ZSlices,numz];
    Timepoints = [Timepoints,numt];
    Start = [Start,start];
end

%% Let's do some quick visualization of the data

figure
hold on
cmap1=cbrewer2('Greens',60);
cmap2=cbrewer2('Blues',60);
cmap3=cbrewer2('Reds',60);
cmaps = {['Greens'],['Reds'],['Blues']};

Timing = [0];

%ylim([0.25,0.6])
%xlim([0,200])
for i = 1:NumEmbryos
    WidthOverTime = [];
    Data = FullData{i};
    %Ratio = Data(:,3)./Data(:,6);
    %Data = [Data,Ratio];
    
    %Separate data by frame
    [~,~,X] = unique(Data(:,1));
    frame_array = accumarray(X,1:size(Data,1),[],@(r){Data(r,:)});
    firstframe = frame_array{1};
    %background = mean([firstframe(:,3);firstframe(:,6)]);
    background = 0;
    timepoints = length(frame_array)+sum(Timing);
    Timing = [Timing;timepoints];
    
    for x = 1:length(frame_array)
        color = cmaps{i};
        cmap=cbrewer2(color,length(frame_array));
        c = cmap(x,:);
        data = frame_array{x};
        position = data(:,4)./4.4;
        zen = data(:,3);
        correctedzen = zen-background;
        correctedzen = max(correctedzen,0);
        cytozen = data(:,6);
        correctedcyto = cytozen-background;
        correctedcyto = max(correctedcyto,0);
        ratio = correctedzen./correctedcyto;
        rposition = position(ratio(:,1)~=Inf,1);
        ratio = ratio(ratio(:,1)~=Inf,1);
        rposition = rposition(ratio(:,1)<1.5,1);
        ratio = ratio(ratio(:,1)<1.5,1);
        ratio = ratio(rposition(:,1)>50,1);
        rposition = rposition(rposition(:,1)>50,1);
        ratio = ratio(rposition(:,1)<175,1);
        rposition = rposition(rposition(:,1)<175,1);
        %ratio = data(:,7);
        [fitresult, gof] = fit(rposition,ratio,'smoothingspline','SmoothingParam',0.000005 );
        yfitted = feval(fitresult,1:250);
        fitforwidth = feval(fitresult,80:160);
        amplitude = max(feval(fitresult,80:160));
        minimum = min(feval(fitresult,80:160));
        width = integrate(fitresult,160,80)./amplitude;
        halfMax = minimum+((amplitude-minimum)/2);
        x1 = find(fitforwidth >= halfMax, 1, 'first');
        x2 = find(fitforwidth >= halfMax, 1, 'last');
        %width = x2-x1;
        %scatter(1:250,yfitted,1,'MarkerEdgeColor',c);
        plot((50:175),yfitted(50:175,:),'Color',c,'LineWidth',1.5)
        WidthOverTime = [WidthOverTime;width];
    end
    timecorrection = Timing(i,1);
    plot((1:length(WidthOverTime))+timecorrection,WidthOverTime(:,1),'Color',cmap(length(frame_array)-2,:),'LineWidth',2)
end
ylim([1,1.5]);
%xlim([50,175])
fontsize(14,"points")
xlabel('Position (Microns)')
ylabel('Nuclear/Cytoplasmic ZenLT')
hold off

%% 

figure
hold on
cmap1=cbrewer2('Greens',60);
cmap2=cbrewer2('Blues',60);
cmap3=cbrewer2('Reds',60);
cmaps = {['Greens'],['Reds'],['Blues']};

Timing = [0];

%ylim([0.25,0.6])
%xlim([0,200])
for i = 1:NumEmbryos
    WidthOverTime = [];
    Data = FullData{i};
    %Ratio = Data(:,3)./Data(:,6);
    %Data = [Data,Ratio];
    
    %Separate data by frame
    [~,~,X] = unique(Data(:,1));
    frame_array = accumarray(X,1:size(Data,1),[],@(r){Data(r,:)});
    firstframe = frame_array{1};
    %background = mean([firstframe(:,3);firstframe(:,6)]);
    background = 0;
    timepoints = length(frame_array)+sum(Timing);
    Timing = [Timing;timepoints];
    
    for x = 1:length(frame_array)
        color = cmaps{i};
        cmap=cbrewer2(color,length(frame_array));
        c = cmap(x,:);
        data = frame_array{x};
        position = data(:,4)./4.4;
        zen = data(:,3);
        correctedzen = zen-background;
        correctedzen = max(correctedzen,0);
        cytozen = data(:,6);
        correctedcyto = cytozen-background;
        correctedcyto = max(correctedcyto,0);
        ratio = correctedzen./correctedcyto;
        rposition = position(ratio(:,1)~=Inf,1);
        ratio = ratio(ratio(:,1)~=Inf,1);
        rposition = rposition(ratio(:,1)<1.5,1);
        ratio = ratio(ratio(:,1)<1.5,1);
        ratio = ratio(rposition(:,1)>50,1);
        rposition = rposition(rposition(:,1)>50,1);
        ratio = ratio(rposition(:,1)<175,1);
        rposition = rposition(rposition(:,1)<175,1);
        %ratio = data(:,7);
        [fitresult, gof] = fit(rposition,ratio,'smoothingspline','SmoothingParam',0.000005 );
        yfitted = feval(fitresult,1:250);
        fitforwidth = feval(fitresult,80:160);
        amplitude = max(feval(fitresult,80:160));
        minimum = min(feval(fitresult,80:160));
        width = integrate(fitresult,160,80)./amplitude;
        halfMax = minimum+((amplitude-minimum)/2);
        x1 = find(fitforwidth >= halfMax, 1, 'first');
        x2 = find(fitforwidth >= halfMax, 1, 'last');
        %width = x2-x1;
        %scatter(1:250,yfitted,1,'MarkerEdgeColor',c);
        plot((50:200),yfitted(50:200,:),'Color',c,'LineWidth',1.5)
        WidthOverTime = [WidthOverTime;width];
    end
    timecorrection = Timing(i,1);
    %plot((1:length(WidthOverTime))+timecorrection,WidthOverTime(:,1),'Color',cmap(length(frame_array)-2,:),'LineWidth',2)
end
ylim([1,1.5]);
%xlim([50,175])
fontsize(14,"points")
xlabel('Position (Microns)')
ylabel('Nuclear/Cytoplasmic ZenLT')
hold off

%% Track nuclei over time, find the midline, and measure Med after background correction
% Will also calculate width and amp over time
% Normalize all experimental data to one another using a reference matrix

WidthPerEmbryo = zeros(NumEmbryos,60);
AmplitudePerEmbryo = zeros(NumEmbryos,60);
AllTrackedNuclei = [];

for embryo = 1:NumEmbryos
    MedData = NC14Data{embryo};
    %Separate data by frame
    [~,~,X] = unique(MedData(:,1));
    frame_array = accumarray(X,1:size(MedData,1),[],@(r){MedData(r,:)});

    ObjectNum_index = 2;
    Area_index = 3;
    MedIntensity_index = 4;
    X_index = 5;
    Y_index = 6;
    NumTimepoints = Timepoints(1,embryo);
    StartTime = Start(1,embryo);
    NumZs = ZSlices(1,embryo);

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

    MidlineTimepoint = 31;

    [CorrectedNuclei] =  findmidline_nodots(FullyTrackedY,FullyTrackedNuclei, ...
    NumTimepoints,MidlineTimepoint,PixelCutoffLow,PixelCutoffHigh,PixelToMicron,BackgroundCutoffLow,BackgroundCutoffHigh);

    [D2_binned,FullyTrackedNucleiRescaled] = rescale(D1_binned,CorrectedNuclei);

    [AmplitudeOverTime,WidthOverTime] = widthamp(NumTimepoints,FullyTrackedNucleiRescaled,StartTime);

    %Need to account for changes in length of movies here...
    %WidthPerEmbryo(embryo,:) = WidthOverTime;
    %AmplitudePerEmbryo(embryo,:) = AmplitudeOverTime;
    %AllTrackedNuclei = [AllTrackedNuclei;FullyTrackedNucleiRescaled];

end

