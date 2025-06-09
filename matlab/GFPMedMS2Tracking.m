Path = 'C:/Users/Susanna Brantley/Brantley Lab Dropbox/Susanna Brantley/Duke/DiTaliaLab/';
%Path = 'D:/Dropbox/Duke/DiTaliaLab/';
%Path = '/Users/susannabrantley/Dropbox/Duke/DiTaliaLab/';

DataPath = 'ImageAnalysis/Datasets/';

%Read in data file for Windows for reference embryo
D1_binned = readmatrix(strcat(Path,DataPath,"BinnedRefEmbryo_D1binned_240429_2.csv"));

%File Labels:
%ImageNo | ObjectNo | Area | Int | Mean | X | Y | Gen | Z 

%Load all embryo data
Embryo1 = readmatrix(strcat(Path,DataPath,"240429_1_wt_GFPMed_ushMS2_cleaned.csv"));
Embryo2 = readmatrix(strcat(Path,DataPath,"240429_2_wt_GFPMed_ushMS2_cleaned.csv"));
Embryo3 = readmatrix(strcat(Path,DataPath,"240429_3_wt_GFPMed_ushMS2_cleaned.csv"));
Embryo4 = readmatrix(strcat(Path,DataPath,"240502_1_wt_GFPMed_ushMS2_cleaned.csv"));

SogEmbryo1 = readmatrix(strcat(Path,DataPath,"230729_1_sog_GFPMed_nodots_cleaned.csv"));
SogEmbryo2 = readmatrix(strcat(Path,DataPath,"230730_1_sog_GFPMed_nodots_cleaned.csv"));
SogEmbryo3 = readmatrix(strcat(Path,DataPath,"230730_2_sog_GFPMed_nodots_cleaned.csv"));

EmbryoData = {Embryo1;Embryo2;Embryo3;Embryo4};

NumEmbryos = length(EmbryoData);

%Set shared parameters
PixelCutoffLow = 0; %300
PixelCutoffHigh = 1024; %750
BackgroundCutoffLow = 50;
BackgroundCutoffHigh = 55;
DistanceThreshold = 3.2;
TimeDistanceThreshold = 30;
PixelToMicron = 4.4;
Micron_cutoff = 75;
EndTime = 60;
NumBins = 15;

%% Set Image Colors

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

%% Identify dataset-sepcific parameters

ZSlices = [];
Timepoints = [];
Start = [];

for e = 1:NumEmbryos
    eData = EmbryoData{e};
    numz = max(eData(:,8));
    numt = max(eData(:,1))/numz;
    start = EndTime-numt+1;
    ZSlices = [ZSlices,numz];
    Timepoints = [Timepoints,numt];
    Start = [Start,start];
end

%% Track nuclei over time, find the midline, and measure Med after background correction
% Will also calculate width and amp over time
% Normalize all experimental data to one another using a reference matrix

WidthPerEmbryo = [];
AmplitudePerEmbryo = [];
AllTrackedMed = [];
AllTrackedDots = [];

for embryo = 1:NumEmbryos
    MedData = EmbryoData{embryo};
    %Separate data by frame
    [~,~,X] = unique(MedData(:,1));
    frame_array = accumarray(X,1:size(MedData,1),[],@(r){MedData(r,:)});

    NumTimepoints = Timepoints(1,embryo);
    StartTime = Start(1,embryo);
    NumZs = ZSlices(1,embryo);

    %Set coordinates for tracking in Z
    ObjectNum_index = 2;
    Area_index = 3;
    MedIntensity_index = 4;
    X_index = 5;
    Y_index = 6;
    DotCount_index = 7;
    DotInt_index = 10;

    NextNuc_index = 11;
    NextArea_index = 12;
    NextMed_index = 13;
    NextX_index = 14;
    NextY_index = 15;
    NextDotCount_index = 16;
    NextDotInt_index = 17;

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
    
    [CorrectedNuclei,CorrectedDots,YAtMaxGradient,XAtMaxGradient] =  findmidline(FullyTrackedY,FullyTrackedX,FullyTrackedNuclei,FullyTrackedDotCount, ...
    NumTimepoints,PixelCutoffLow,PixelCutoffHigh,PixelToMicron,BackgroundCutoffLow,BackgroundCutoffHigh);

    [D2_binned,FullyTrackedNucleiRescaled] = rescale(D1_binned,NumTimepoints,CorrectedNuclei);

    %WidthPerEmbryo = [WidthPerEmbryo,WidthOverTime(NumTimepoints-39:NumTimepoints,2)];
    %AmplitudePerEmbryo = [AmplitudePerEmbryo,AmplitudeOverTime(NumTimepoints-39:NumTimepoints,1)];
    AllTrackedMed20_60 = [FullyTrackedNucleiRescaled(:,1),FullyTrackedNucleiRescaled(:,NumTimepoints-39:NumTimepoints)];
    AllTrackedMed = [AllTrackedMed;AllTrackedMed20_60];
    AllTrackedDots20_60 = [CorrectedDots(:,1),CorrectedDots(:,NumTimepoints-39:NumTimepoints)];
    AllTrackedDots = [AllTrackedDots;AllTrackedDots20_60];
    
    % figure
    % hold on
    % scatter(XAtMaxGradient,YAtMaxGradient);
    % ylabel('Microns from Midline');
    % xlabel('Microns along AP Axis');
    % hold off

end

%% Plot width and amp over time with standard deviation

avgWidth = mean(WidthPerEmbryo,2);
avgAmp = mean(AmplitudePerEmbryo,2);
stdWidth = std(WidthPerEmbryo,0,2);
stdAmp = std(AmplitudePerEmbryo,0,2);

x = 21:60;
y = avgWidth';
sd = stdWidth';
f=figure();
hold on
plot(x,y)
fontsize(gcf,14,'points')
patch([x fliplr(x)], [y-sd  fliplr(y+sd)],[0.85, 0.11, 0.38],'FaceAlpha',0.2, 'EdgeColor','none' )
plot(x,y,'LineWidth',2,'Color',mymagenta,'LineStyle','--')
ylabel("Width at Half Max (um)")
xlabel("Time (min)")
ylim([0,50])
saveas(f,strcat(FigurePath,'widthovertime'),'svg');
hold off

x = 21:60;
y = avgAmp';
sd = stdAmp';
f=figure();
hold on
fontsize(gcf,14,'points')
plot(x,y)
patch([x fliplr(x)], [y-sd  fliplr(y+sd)],[0.12, 0.53, 0.9],'FaceAlpha',0.2, 'EdgeColor','none' )
plot(x,y,'LineWidth',2,'Color',myblue,'LineStyle','--')
ylabel("Amplitude")
xlabel("Time (min)")
ylim([0,0.03])
saveas(f,strcat(FigurePath,'amplitudeovertime'),'svg');
hold off

%% Bin Med by distance and plot over time

distedges = 0:5:45;
binnedbydist = discretize(AllTrackedMed(:,1),distedges);
BinnedByDist = [];

for distbin = 1:length(distedges)-1
    datainbin = [];
    for row = 1:length(AllTrackedMed)
        bin = binnedbydist(row,1);
        if bin == distbin
            data = AllTrackedMed(row,2:end);
            datainbin = [datainbin;data];
        else
            continue
        end
    end
    meandata = mean(datainbin,1);
    BinnedByDist = [BinnedByDist;meandata];
end

cmap=cbrewer2('Blues',15);
f=figure();
hold on
fontsize(gcf,14,'points')
for bin = 1:9
    xdata = (21:60)';
    ydata = BinnedByDist(bin,:)';
    ft = fit(xdata, ydata,'SmoothingSpline','SmoothingParam',0.05);
    yfitted = feval(ft,xdata);
    c = cmap(15-bin,:);
    scatter(xdata,ydata,10,'MarkerEdgeColor',c,'MarkerFaceColor',c)
    plot(xdata,yfitted,'color',c,'LineWidth',1)
end
xlabel('Time (min)')
ylabel('Nuclear GFP-Med')
saveas(f,strcat(FigurePath,'binnedmedovertime'),'svg');
hold off

%% Plot the derivative over time at each position bin to see the plateau
cmap=cbrewer2('Blues',15);
f=figure();
hold on
fontsize(gcf,14,'points')
for bin = 1:9
    xdata = (21:60)';
    ydata = BinnedByDist(bin,:)';
    ft = fit(xdata, ydata,'SmoothingSpline','SmoothingParam',0.05);
    rate = differentiate(ft,xdata);
    c = cmap(15-bin,:);
    scatter(xdata,rate,10,'MarkerEdgeColor',c,'MarkerFaceColor',c)
    plot(xdata,rate,'color',c,'LineWidth',2)
end
xlabel('Time (min)')
ylabel('d/dt(GFPMed)')
saveas(f,strcat(FigurePath,'ratemedovertime'),'svg');
hold off

%% Plot dots over time 

% this function defines when a nucleus first starts to express an MS2 dot)
% 40 is the number of timepoints, 21 is the start time, and 60 is the end
% time
first_transcript = firsttranscript(AllTrackedMed, AllTrackedDots, 40,21,60,0);

Transcribing_Nuclei = first_transcript(first_transcript(:,1)<40,:);
MS2_ontime = Transcribing_Nuclei(:,1)+20;
Med_level_at_on = Transcribing_Nuclei(:,2);
Med_integral_at_on = Transcribing_Nuclei(:,3);
Distance = Transcribing_Nuclei(:,5);

% the data are very noisy, so let's bin by distance and average
distedges = 0:3:51;
binnedbydist = discretize(Distance(:,1),distedges,'IncludedEdge','left');
BinnedByDistTimeOn = zeros(1,17);
BinnedByDistMedAtOn = zeros(1,17);
BinnedByDistIntMedAtOn = zeros(1,17); 

BinnedByDistTimeOnStd = zeros(1,17);
BinnedByDistMedAtOnStd = zeros(1,17);
BinnedByDistIntMedAtOnStd = zeros(1,17); 

for bin = 1:17
    ontimeinbin =[];
    medinbin = [];
    intmedinbin = [];
    for n = 1:length(binnedbydist)
        checkbin = binnedbydist(n,1);
        if checkbin == bin
            ontimeinbin = [ontimeinbin;MS2_ontime(n,1)];
            medinbin = [medinbin;Med_level_at_on(n,1)];
            intmedinbin = [intmedinbin;Med_integral_at_on(n,1)];
        else
            continue
        end
    end
    BinnedByDistTimeOn(1,bin)=mean(ontimeinbin);
    BinnedByDistMedAtOn(1,bin)=mean(medinbin);
    BinnedByDistIntMedAtOn(1,bin)=mean(intmedinbin);

    BinnedByDistTimeOnStd(1,bin)=std(ontimeinbin);
    BinnedByDistMedAtOnStd(1,bin)=std(medinbin);
    BinnedByDistIntMedAtOnStd(1,bin)=std(intmedinbin);
end

x = 3:3:51;
y = BinnedByDistTimeOn;
sd = BinnedByDistTimeOnStd;
f=figure();
hold on
fontsize(gcf,14,'points')
patch([x fliplr(x)], [y-sd  fliplr(y+sd)],[0.85, 0.11, 0.38],'FaceAlpha',0.2, 'EdgeColor','none' )
scatter(x,y,40,'MarkerEdgeColor',mymagenta,'MarkerFaceColor',mymagenta)
ylabel("Time On (min)")
xlabel("Distance from Midline (um)")
xlim([0,50])
ylim([20,60])
P = polyfit(x,y,1);
yfit = P(1)*x+P(2);
plot(x,yfit,'r--','Color',mymagenta,'LineWidth',1);
SStot = sum((y-mean(y)).^2); 
SSres = sum((y-yfit).^2); 
Rsq = 1-SSres/SStot;
text(6,57.5,['R2 = ', num2str(Rsq)],'FontSize',14);
saveas(f,strcat(FigurePath,'ushms2ontime'),'svg');
hold off

x = 3:3:51;
y = BinnedByDistMedAtOn;
sd = BinnedByDistMedAtOnStd;
f=figure();
hold on
fontsize(gcf,14,'points')
patch([x fliplr(x)], [y-sd  fliplr(y+sd)],[0.12, 0.53, 0.9],'FaceAlpha',0.2, 'EdgeColor','none' )
scatter(x,y,40,'MarkerEdgeColor',myblue,'MarkerFaceColor',myblue)
ylabel("Nuclear Medea")
xlabel("Distance from Midline (um)")
xlim([0,50])
ylim([0,0.015])
P = polyfit(x,y,1);
yfit = P(1)*x+P(2);
plot(x,yfit,'r--','Color',myblue,'LineWidth',1);
SStot = sum((y-mean(y)).^2); 
SSres = sum((y-yfit).^2); 
Rsq = 1-SSres/SStot;
text(6,0.014,['R2 = ', num2str(Rsq)],'FontSize',14);
saveas(f,strcat(FigurePath,'ushms2onmed'),'svg');
hold off

x = 3:3:51;
y = BinnedByDistIntMedAtOn;
sd = BinnedByDistIntMedAtOnStd;
f=figure();
hold on
%plot(x,y)
fontsize(gcf,14,'points')
patch([x fliplr(x)], [y-sd  fliplr(y+sd)],[0.12, 0.53, 0.9],'FaceAlpha',0.2, 'EdgeColor','none' )
scatter(x,y,40,'MarkerEdgeColor',myblue,'MarkerFaceColor',myblue)
ylabel("Integrated Nuclear Medea")
xlabel("Distance from Midline (um)")
xlim([0,50])
ylim([0,0.15])
[P,s] = polyfit(x,y,1);
yfit = P(1)*x+P(2);
plot(x,yfit,'r--','Color',myblue,'LineWidth',1);
Rsq = 1 - s.normr^2 / norm(y-mean(y))^2;
% SStot = sum((y-mean(y)).^2); 
% SSres = sum((y-yfit).^2); 
% Rsq = 1-SSres/SStot;
text(6,0.14,['R2 = ', num2str(Rsq)],'FontSize',14);
saveas(f,strcat(FigurePath,'ushms2onintmed'),'svg');
hold off

%% Plot integrated Med over time

first_transcript = firsttranscript(AllTrackedMed, AllTrackedDots, 40,21,60,2);

TimeOn = first_transcript(:,1);
Position = first_transcript(:,5);
TrackedMed = AllTrackedMed(:,2:end);
IntTrackedMed = cumtrapz(1:40,TrackedMed,2);

distedges = [0:5:50];
binnedbydist = discretize(AllTrackedMed(:,1),distedges);
BinnedByDist = [];

for distbin = 1:length(distedges)-1
    datainbin = [];
    for row = 1:length(AllTrackedMed)
        bin = binnedbydist(row,1);
        if bin == distbin
            data = IntTrackedMed(row,:);
            datainbin = [datainbin;data];
        else
            continue
        end
    end
    meandata = mean(datainbin,1);
    BinnedByDist = [BinnedByDist;meandata];
end

cmap=cbrewer2('Blues',15);
f=figure();
hold on
fontsize(gcf,14,'points')
ylim([0,0.4]);
for bin = 1:10
    xdata = [21:60]';
    ydata = BinnedByDist(bin,:)';
    ft = fit(xdata, ydata,'SmoothingSpline','SmoothingParam',0.05);
    yfitted = feval(ft,xdata);
    c = cmap(15-bin,:);
    scatter(xdata,ydata,10,'MarkerEdgeColor',c,'MarkerFaceColor',c)
    plot(xdata,yfitted,'color',c,'LineWidth',1)
end
xlabel('Time (min)')
ylabel('Integrated Nuclear GFP-Med')
saveas(f,strcat(FigurePath,'binnedmedintovertime'),'svg');
hold off

%% Experimental thing for Stefano
distedges = 0:5:45;
binnedbydist = discretize(AllTrackedMed(:,1),distedges);
BinnedByDist = [];

for distbin = 1:length(distedges)-1
    datainbin = [];
    for row = 1:length(AllTrackedMed)
        bin = binnedbydist(row,1);
        if bin == distbin
            data = AllTrackedMed(row,2:end);
            datainbin = [datainbin;data];
        else
            continue
        end
    end
    meandata = mean(datainbin,1);
    BinnedByDist = [BinnedByDist;meandata];
end

figure
hold on
Slopes = [];
for bin = 1:8
    xdata = (35:50)';
    ydata = BinnedByDist(bin,15:30)';
    P = polyfit(xdata, ydata,1)
    plot(xdata,ydata)
    slope = P(1);
    Slopes = [Slopes,slope];
end
hold off

first_transcript = firsttranscript(AllTrackedMed, AllTrackedDots, 40,21,60,2);

Transcribing_Nuclei = first_transcript(first_transcript(:,1)<40,:);
MS2_ontime = Transcribing_Nuclei(:,1)+20;
Med_level_at_on = Transcribing_Nuclei(:,2);
Med_integral_at_on = Transcribing_Nuclei(:,3);
Distance = Transcribing_Nuclei(:,5);

% the data are very noisy, so let's bin by distance and average
distedges = 0:5:40;
binnedbydist = discretize(Distance(:,1),distedges);
BinnedByDistTimeOn = zeros(1,8);
BinnedByDistMedAtOn = zeros(1,8);
BinnedByDistIntMedAtOn = zeros(1,8); 

BinnedByDistTimeOnStd = zeros(1,8);
BinnedByDistMedAtOnStd = zeros(1,8);
BinnedByDistIntMedAtOnStd = zeros(1,8); 

for bin = 1:8
    ontimeinbin =[];
    medinbin = [];
    intmedinbin = [];
    for n = 1:length(binnedbydist)
        checkbin = binnedbydist(n,1);
        if checkbin == bin
            ontimeinbin = [ontimeinbin;MS2_ontime(n,1)];
            medinbin = [medinbin;Med_level_at_on(n,1)];
            intmedinbin = [intmedinbin;Med_integral_at_on(n,1)];
        else
            continue
        end
    end
    BinnedByDistTimeOn(1,bin)=mean(ontimeinbin);
    BinnedByDistMedAtOn(1,bin)=mean(medinbin);
    BinnedByDistIntMedAtOn(1,bin)=mean(intmedinbin);

    BinnedByDistTimeOnStd(1,bin)=std(ontimeinbin);
    BinnedByDistMedAtOnStd(1,bin)=std(medinbin);
    BinnedByDistIntMedAtOnStd(1,bin)=std(intmedinbin);
end 

x = Slopes;
y = BinnedByDistIntMedAtOn;
sd = BinnedByDistIntMedAtOnStd;
f=figure();
hold on
%plot(x,y)
fontsize(gcf,14,'points')
patch([x fliplr(x)], [y-sd  fliplr(y+sd)],[0.12, 0.53, 0.9],'FaceAlpha',0.2, 'EdgeColor','none' )
scatter(x,y,40,'MarkerEdgeColor',myblue,'MarkerFaceColor',myblue)
ylabel("Integrated Nuclear Medea")
xlabel("Slope")
%xlim([5,40])
%ylim([0,0.15])
[P,s] = polyfit(sqrt(x),y,1);
yfit = P(1)*sqrt(x)+P(2);
plot(x,yfit,'r--','Color',myblue,'LineWidth',1);
Rsq = 1 - s.normr^2 / norm(y-mean(y))^2;
% SStot = sum((y-mean(y)).^2); 
% SSres = sum((y-yfit).^2); 
% Rsq = 1-SSres/SStot;
text(6,0.14,['R2 = ', num2str(Rsq)],'FontSize',14);
%saveas(f,strcat(FigurePath,'ushms2onintmed'),'svg');
hold off

x = Slopes;
y = BinnedByDistMedAtOn;
sd = BinnedByDistMedAtOnStd;
f=figure();
hold on
%plot(x,y)
fontsize(gcf,14,'points')
patch([x fliplr(x)], [y-sd  fliplr(y+sd)],[0.12, 0.53, 0.9],'FaceAlpha',0.2, 'EdgeColor','none' )
scatter(x,y,40,'MarkerEdgeColor',myblue,'MarkerFaceColor',myblue)
ylabel("Instantaneous Nuclear Medea")
xlabel("Slope")
%xlim([5,40])
%ylim([0,0.15])
[P,s] = polyfit(sqrt(x),y,1);
yfit = P(1)*sqrt(x)+P(2);
plot(x,yfit,'r--','Color',myblue,'LineWidth',1);
Rsq = 1 - s.normr^2 / norm(y-mean(y))^2;
% SStot = sum((y-mean(y)).^2); 
% SSres = sum((y-yfit).^2); 
% Rsq = 1-SSres/SStot;
text(6,0.14,['R2 = ', num2str(Rsq)],'FontSize',14);
%saveas(f,strcat(FigurePath,'ushms2onintmed'),'svg');
hold off