%% Set Image Colors

Path = 'C:/Users/Susanna Brantley/Brantley Lab Dropbox/Susanna Brantley/Duke/DiTaliaLab/';
%Path = 'D:/Dropbox/Duke/DiTaliaLab/';
%Path = '/Users/susannabrantley/Dropbox/Duke/DiTaliaLab/';

DataPath = 'ImageAnalysis/Datasets/';

mymagenta = '#D81B60';
myblue = '#1E88E5';
mygreen = '#004D40';
myyellow = '#FFC107';

FigurePath = strcat(Path,'Presentations/PlacedImages/');

set(0,'DefaultAxesFontName','Arial');

% load heatmap settings
mymagentamap = load('mygreenmap.mat').mygreenmap;
mybluemap = load('mybluemap.mat').mybluemap;
mygreenmap = load('mymagentamap.mat').mymagentamap;
myyellowmap = load('myyellowmap.mat').myyellowmap;

s.matlab.fonts.editor.code.Name.TemporaryValue = 'Arial';

%% Load Data

WTData = load("AllDataWTLive.mat", 'AllTrackedNuclei');
WTData = WTData.AllTrackedNuclei;

HetData = load("AllDataSogHetLive.mat", 'AllTrackedNucleiSog');
HetData = HetData.AllTrackedNucleiSog;

SogData = load("AllDataSogLive.mat", 'AllTrackedNucleiSog');
SogData = SogData.AllTrackedNucleiSog;

AllWTDataUsh = load("AllDataWTtup.mat", 'AllWTNucleiData');
AllWTTimingUsh = load("AllTimingWTtup.mat", 'AllWTNucleiTiming');
AllWTDataUsh = AllWTDataUsh.AllWTNucleiData;
AllWTTimingUsh = AllWTTimingUsh.AllWTNucleiTiming;

AllMutDataUsh = load("AllDataSogC15.mat", 'AllMutantNucleiData');
AllMutTimingUsh = load("AllTimingSogC15.mat", 'AllMutantNucleiTiming');
AllMutDataUsh = AllMutDataUsh.AllMutantNucleiData;
AllMutTimingUsh = AllMutTimingUsh.AllMutantNucleiTiming;

AllHetDataUsh = load("AllDataSogHettup.mat", 'AllMutantNucleiData');
AllHetTimingUsh = load("AllTimingSogHettup.mat", 'AllMutantNucleiTiming');
AllHetDataUsh = AllHetDataUsh.AllMutantNucleiData;
AllHetTimingUsh = AllHetTimingUsh.AllMutantNucleiTiming;

AllZenDataUsh = load("AllDataZenC15.mat", 'AllMutantNucleiData');
AllZenTimingUsh = load("AllTimingZenC15.mat", 'AllMutantNucleiTiming');
AllZenDataUsh = AllZenDataUsh.AllMutantNucleiData;
AllZenTimingUsh = AllZenTimingUsh.AllMutantNucleiTiming;

%% Bin Med Data

% I will change the structure of the data so that nuclei are organized by
% position and timing information is included in a new column
binedges = 0:5:75;
WTPositionBins = discretize(WTData(:,1),binedges,'IncludedEdge','left');
AllWTData = [WTPositionBins,WTData];

BinnedMed1 = zeros(size(binedges,2)-1,60);
BinnedUpMed1 = zeros(size(binedges,2)-1,60);
BinnedLowMed1 = zeros(size(binedges,2)-1,60);

for binno = 1:max(WTPositionBins)
    nucinbin = AllWTData(AllWTData(:,1)==binno,3:end);
    for t = 1:60
        nucatt = nucinbin(:,t);
        avgmed = mean(nucatt);
        devmed = std(nucatt)./sqrt(length(nucatt));
        BinnedMed1(binno,t) = avgmed;
        BinnedUpMed1(binno,t) = avgmed+devmed;
        BinnedLowMed1(binno,t) = avgmed-devmed;
    end
end

%binedges = 0:5:100;
SogPositionBins = discretize(HetData(:,1),binedges,'IncludedEdge','left');
AllSogData = [SogPositionBins,HetData];

BinnedMedSog1 = zeros(size(binedges,2)-1,60);
BinnedUpMedSog1 = zeros(size(binedges,2)-1,60);
BinnedLowMedSog1 = zeros(size(binedges,2)-1,60);

for binno = 1:max(SogPositionBins)
    nucinbin = AllSogData(AllSogData(:,1)==binno,3:end);
    for t = 1:60
        nucatt = nucinbin(:,t);
        avgmed = mean(nucatt);
        devmed = std(nucatt)./sqrt(length(nucatt));
        BinnedMedSog1(binno,t) = avgmed;
        BinnedUpMedSog1(binno,t) = avgmed+devmed;
        BinnedLowMedSog1(binno,t) = avgmed-devmed;
    end
end

BinnedMed = [];
BinnedMedRate = [];
cmap = cbrewer2('Blues',25);
figure
hold on
for x = 1:size(binedges,2)-1
    xdata = 1:60;
    ydata = BinnedMed1(x,:);
    ft = fit(xdata',ydata','SmoothingSpline','SmoothingParam',0.01);
    fittedmed = feval(ft,1:60);
    rateofmed = differentiate(ft,1:60);
    BinnedMed = [BinnedMed,fittedmed];
    BinnedMedRate = [BinnedMedRate,rateofmed];
    c=cmap(25-x,:);
    scatter(xdata,ydata,'MarkerEdgeColor',c)
    plot(xdata,fittedmed,'Color',c);
end
ylim([0,0.04])
%xlim([20,60])
hold off

BinnedMedSog = [];
BinnedMedRateSog = [];
cmap = cbrewer2('Blues',25);
figure
hold on
for x = 1:size(binedges,2)-1
    xdata = 1:60;
    ydata = BinnedMedSog1(x,:);
    ft = fit(xdata',ydata','SmoothingSpline','SmoothingParam',0.01);
    rateofmed = differentiate(ft,1:60);
    fittedmed = feval(ft,1:60);
    BinnedMedSog = [BinnedMedSog,fittedmed];
    BinnedMedRateSog = [BinnedMedRateSog,rateofmed];
    c=cmap(25-x,:);
    scatter(xdata,ydata,'MarkerEdgeColor',c)
    plot(xdata,fittedmed,'Color',c);
end
ylim([0,0.04])
%xlim([20,60])
hold off

% BinnedMed = BinnedMed-mean(BinnedMed(:,1:20),2);
% BinnedMedSog = BinnedMedSog-mean(BinnedMedSog(:,1:20),2);
% 
BinnedMed = max(BinnedMed,0);
BinnedMedSog = max(BinnedMedSog,0);

BinnedMed = BinnedMed;
BinnedMedSog = BinnedMedSog;

%% Plot rate of Med 

cmap = cbrewer2('Blues',25);
figure
hold on
fontsize(14,"points")
for x = 1:size(binedges,2)-2
    xdata = 1:60;
    ydata = BinnedMedRate(:,x);
    c = cmap(25-x,:);
    plot(xdata',ydata,'Color',c,'LineWidth',2,'LineStyle','--')
end
xlabel('Time (Minutes)')
ylabel('dMed/dt')
ylim([0,0.0025])
xlim([20,60])
hold off

cmap = cbrewer2('Oranges',25);
figure
hold on
fontsize(14,"points")
for x = 1:size(binedges,2)-2
    xdata = 1:60;
    ydata = BinnedMedRateSog(:,x);
    c = cmap(25-x,:);
    plot(xdata',ydata,'Color',c,'LineWidth',2,'LineStyle','--')
end
xlabel('Time (Minutes)')
ylabel('dMed/dt')
ylim([0,0.0025])
xlim([20,60])
hold off


%% Find width and amp 

% WidthOverTime = [];
% AmplitudeOverTime = [];
% SogWidthOverTime = [];
% SogAmplitudeOverTime = [];
% ZenWidthOverTime = [];
% ZenAmplitudeOverTime = [];
% figure
% hold on
% for t = 20:60
%     timepointdata = BinnedMed(t,:);
%     amplitude = max(timepointdata);
%     ft = fit([0:3:60]',timepointdata','SmoothingSpline','SmoothingParam',0.0005);
%     %plot(1:100,feval(ft,1:100))
%     yfitted = feval(ft,1:60);
%     halfMax = amplitude/2;
%     width = trapz(timepointdata)/amplitude;
%     %width = find(yfitted >= halfMax, 1, 'last');
%     WidthOverTime = [WidthOverTime,width*2];
%     AmplitudeOverTime = [AmplitudeOverTime,amplitude];
% 
%     sogtimepointdata = BinnedMedSog(t,1:20);
%     sogamplitude = max(sogtimepointdata(:,1));
%     sogft = fit([0:3:57]',sogtimepointdata','SmoothingSpline','SmoothingParam',0.0005);
%     yfitted = feval(sogft,1:60);
%     scatter(0:3:57,sogtimepointdata)
%     plot(1:100,feval(sogft,1:100))
%     halfMax = sogamplitude/2;
%     sogwidth = trapz(sogtimepointdata)/sogamplitude;
%     %sogwidth = find(yfitted >= halfMax, 1, 'last');
%     SogWidthOverTime = [SogWidthOverTime,sogwidth*2];
%     SogAmplitudeOverTime = [SogAmplitudeOverTime,sogamplitude];
% end
% hold off
% 
% figure
% hold on
% fontsize(14,'points')
% plot(20:60,WidthOverTime,'Color',myblue,'LineWidth',2)
% plot(20:60,SogWidthOverTime,'Color',mymagenta,'LineWidth',2)
% xlabel('Time (Minutes)')
% ylabel('Width')
% hold off
% 
% figure
% hold on
% fontsize(14,'points')
% plot(20:60,AmplitudeOverTime,'Color',myblue,'LineWidth',2)
% plot(20:60,SogAmplitudeOverTime,'Color',mymagenta,'LineWidth',2)
% xlabel('Time (Minutes)')
% ylabel('Amplitude')
% hold off


%% Calculate integrated Med

cmap = cbrewer2('Blues',25);

MeanCumInt = [];
inttime = 60; %optinttime;
offtime = 0;
figure
hold on
for x = 1:size(binedges,2)-1
    avgdata = BinnedMed(:,x);
    avgdata = avgdata;
    int = [zeros(offtime,1);cumtrapz(avgdata(offtime+1:inttime,1))];
    for t = (inttime+1):60
        intwindow = t-inttime;
        timeint = max(cumtrapz(avgdata(intwindow:t,1)));
        int = [int;timeint];
    end
    c=cmap(25-x,:);
    scatter(1:60,int,'MarkerEdgeColor',c)
    MeanCumInt = [MeanCumInt,int];
end
%xlim([20,60])
ylim([0,1])
hold off



SogMeanCumInt = [];
%inttime = 10; %optinttime;
%offtime = 0;
figure
hold on
for x = 1:size(binedges,2)-1
    avgdata = BinnedMedSog(:,x);
    avgdata = avgdata;
    int = [zeros(offtime,1);cumtrapz(avgdata(offtime+1:inttime,1))];
    for t = (inttime+1):60
        intwindow = t-inttime;
        timeint = max(cumtrapz(avgdata(intwindow:t,1)));
        int = [int;timeint];
    end
    c=cmap(25-x,:);
    scatter(1:60,int,'MarkerEdgeColor',c)
    SogMeanCumInt = [SogMeanCumInt,int];
end
%xlim([20,60])
ylim([0,1])
hold off

%% Bin gene expression

NucByPositionUsh = cell(1,21);
counter = 0;
for d = binedges
    counter = counter+1;
    PositionData = [];
    for i = 1:length(AllWTDataUsh)
        age = AllWTTimingUsh(i,1);
        data = AllWTDataUsh{i};
        for n = 1:length(data)
            position = data(n,1);
            if position>=d &&  position<d+(binedges(1,2)-binedges(1,1))
                PositionData = [PositionData; age, data(n,:)];
            else
                continue
            end
        end
    end
    NucByPositionUsh{counter} = PositionData;
end

MutNucByPositionUsh = cell(1,21);
counter = 0;
for d = binedges
    counter = counter+1;
    PositionData = [];
    for i = 1:length(AllHetDataUsh)
        age = AllHetTimingUsh(i,1);
        data = AllHetDataUsh{i};
        for n = 1:length(data)
            position = data(n,1);
            if position>=d &&  position<d+(binedges(1,2)-binedges(1,1))
                PositionData = [PositionData; age, data(n,:)];
            else
                continue
            end
        end
    end
    MutNucByPositionUsh{counter} = PositionData;
end

cmap=cbrewer2('Blues',25);

%timebinedges = [0,20,25:5:60];

FractionOn = [];
HeatmapData = [];
RateOn = [];
figure
hold on
fontsize(14,"points")
for x = 1:size(binedges,2)-1%21
    nucatpos = NucByPositionUsh{x};
    [~,~,X] = unique(nucatpos(:,1));
    frame_array = accumarray(X,1:size(nucatpos,1),[],@(r){nucatpos(r,:)});
    collectedtimes = [0];
    collectedfracs =[0];
    for i = 1:length(frame_array)
        timebinned = frame_array{i};
        numon = sum(timebinned(:,4)>0);
        totalnuc = length(timebinned);
        fractionon = numon/totalnuc;
        time = timebinned(1,1);
        collectedtimes = [collectedtimes;time];
        collectedfracs = [collectedfracs;fractionon];
    end
    %collectedfracs = collectedfracs(collectedtimes(:,1)>20,:);
    %collectedtimes = collectedtimes(collectedtimes(:,1)>20,:);
    collectedtimes = [0;collectedtimes];
    collectedfracs = [0;collectedfracs];
    [fitresult, gof] = fit( collectedtimes,collectedfracs,'smoothingspline','SmoothingParam', 0.005 );
    yfitted = feval(fitresult,1:60);
    HeatmapData = [HeatmapData,yfitted];
    [logitCoef,dev] = glmfit(collectedtimes,collectedfracs,'binomial','logit');
    logitFit = glmval(logitCoef,1:60,'logit');
    %dev;
    %logitFit(logitFit(:,1)>0.9999,1)=0;
    %rateon = differentiate(logitCoef,1:60);
    FractionOn = [FractionOn,logitFit];
    %RateOn = [RateOn,rateon];
    c = cmap(25-x,:);
    scatter(collectedtimes,collectedfracs,'MarkerEdgeColor',c,'MarkerFaceColor',c)
    plot(1:60,logitFit,'LineWidth',2,'Color',c)
    ylim([0 1])
end
xlabel('Time (minutes)')
ylabel('Fraction Nuclei On')
xlim([20,60])
hold off

MutFractionOn = [];
MutHeatmapData = [];
figure
hold on
for x = 1:size(binedges,2)-1%21
    nucatpos = MutNucByPositionUsh{x};
    [~,~,X] = unique(nucatpos(:,1));
    frame_array = accumarray(X,1:size(nucatpos,1),[],@(r){nucatpos(r,:)});
    collectedtimes = [0];
    collectedfracs =[0];
    for i = 1:length(frame_array)
        timebinned = frame_array{i};
        numon = sum(timebinned(:,4)>0);
        totalnuc = length(timebinned);
        fractionon = numon/totalnuc;
        time = timebinned(1,1);
        collectedtimes = [collectedtimes;time];
        collectedfracs = [collectedfracs;fractionon];
    end
    %collectedfracs = collectedfracs(collectedtimes(:,1)>10,:);
    %collectedtimes = collectedtimes(collectedtimes(:,1)>10,:);
    collectedtimes = [0;collectedtimes];
    collectedfracs = [0;collectedfracs];
    [fitresult, gof] = fit( collectedtimes,collectedfracs,'smoothingspline','SmoothingParam', 0.005 );
    yfitted = feval(fitresult,1:60);
    MutHeatmapData = [MutHeatmapData,yfitted];
    [logitCoef,dev] = glmfit(collectedtimes,collectedfracs,'binomial','logit');
    logitFit = glmval(logitCoef,1:60,'logit');
    %logitFit(logitFit(:,1)>0.9999,1)=0;
    MutFractionOn = [MutFractionOn,logitFit];
    c = cmap(25-x,:);
    scatter(collectedtimes,collectedfracs,'MarkerEdgeColor',c,'MarkerFaceColor',c)
    plot(1:60,logitFit,'Color',c)
end
hold off

figure
heatmap(FractionOn)
caxis([0, 1]);
figure
heatmap(MutFractionOn)
caxis([0, 1]);


%%
% FractionOn = BinnedMS2Data_final(:,:);

%%
% BinnedMed = BinnedMed';
% BinnedMedSog = BinnedMedSog';
% BinnedMedZen = BinnedMedZen';

%% Plot probability v. med levels

lagtime = 0;
ka = 1;

ProbForLogistic = [];
SogProbForLogistic = [];
ZenProbForLogistic = [];

MadForLogistic = [];
SogMadForLogistic = [];
ZenMadForLogistic = [];

figure 
hold on
fontsize(14,"points")

for x = 1:size(binedges,2)-1
    % probability = FractionOn(:,x);
    % ProbForLogistic =[ProbForLogistic;probability];
    % instmad = [zeros(lagtime,1);BinnedMed(1:(end-lagtime),x)];
    % MadForLogistic = [MadForLogistic;instmad];
    % cmap = cbrewer2('Blues',30);
    % c=cmap(30-x,:);
    % scatter(instmad,probability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

    sogprobability = MutFractionOn(:,x);
    SogProbForLogistic =[SogProbForLogistic;sogprobability];
    soginstmad = [zeros(lagtime,1);BinnedMedSog(1:(end-lagtime),x)];
    SogMadForLogistic = [SogMadForLogistic;soginstmad];
    cmap = cbrewer2('PuRd',30);
    c = cmap(25-x,:);
    scatter(soginstmad,sogprobability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

    % zenprobability = ZenFractionOn(:,x);
    % ZenProbForLogistic =[ZenProbForLogistic;zenprobability];
    % zeninstmad = [zeros(lagtime,1);BinnedMedZen(1:(end-lagtime),x)];
    % ZenMadForLogistic = [ZenMadForLogistic;zeninstmad];
    % cmap = cbrewer2('Greens',30);
    % c = cmap(25-x,:);
    % scatter(zeninstmad,zenprobability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

end
ylim([0,1])
xlabel('Instantaneous Mad')
ylabel('Probability')

startPoints = [1 1];
s = fitoptions('Method','NonlinearLeastSquares',... %
    'Lower',[0 0],...
    'Upper',[inf inf],...
    'Startpoint',startPoints);

HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);

% [ffun,gofrr] = fit(MadForLogistic,ProbForLogistic,HillEqn);
% ffun
% yfit=feval(ffun,0:0.0001:10); %Fitted function
% plot(0:0.0001:10,yfit,'Color',myblue,'LineWidth',2);

[ffun,gofrr] = fit(SogMadForLogistic,SogProbForLogistic,HillEqn);
ffun
gofrr
yfit=feval(ffun,0:0.0001:10); %Fitted function
plot(0:0.0001:10,yfit,'Color',mymagenta,'LineWidth',2);

% [ffun,gofrr] = fit(ZenMadForLogistic,ZenProbForLogistic,HillEqn);
% ffun
% yfit=feval(ffun,0:0.0001:10); %Fitted function
% plot(0:0.0001:10,yfit,'Color',mygreen,'LineWidth',2);

xlim([0,0.05])
hold off

IntForLogistic = [];
SogIntForLogistic = [];
ZenIntForLogistic = [];

figure 
hold on
fontsize(14,"points")
for x = 1:size(binedges,2)-1
    % probability = FractionOn(1:60,x);
    % instmad = [zeros(lagtime,1);MeanCumInt(1:(end-lagtime),x)];
    % IntForLogistic = [IntForLogistic;instmad];
    % cmap = cbrewer2('Blues',30);
    % c=cmap(30-x,:);
    % scatter(instmad./ka,probability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

    sogprobability = MutFractionOn(1:60,x);
    soginstmad = [zeros(lagtime,1);SogMeanCumInt(1:(end-lagtime),x)];
    SogIntForLogistic = [SogIntForLogistic;soginstmad];
    cmap = cbrewer2('PuRd',30);
    c = cmap(30-x,:);
    scatter(soginstmad,sogprobability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

    % zenprobability = ZenFractionOn(1:60,x);
    % zeninstmad = [zeros(lagtime,1);ZenMeanCumInt(1:(end-lagtime),x)];
    % ZenIntForLogistic = [ZenIntForLogistic;zeninstmad];
    % cmap = cbrewer2('Greens',30);
    % c = cmap(30-x,:);
    % scatter(zeninstmad,zenprobability,'MarkerEdgeColor',c,'MarkerFaceColor',c)
end


startPoints = [1 1];
s = fitoptions('Method','NonlinearLeastSquares',... %
    'Lower',[0 0],...
    'Upper',[inf inf],...
    'Startpoint',startPoints);

HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);

% [ffun,gofrr] = fit(IntForLogistic,ProbForLogistic,HillEqn);
% ffun
% yfit=feval(ffun,0:0.0001:10); %Fitted function
% plot(0:0.0001:10,yfit,'Color',myblue,'LineWidth',2);

[ffun,gofrr] = fit(SogIntForLogistic,SogProbForLogistic,HillEqn);
ffun
gofrr
yfit=feval(ffun,0:0.0001:10); %Fitted function
plot(0:0.0001:10,yfit,'Color',mymagenta,'LineWidth',2);

% [ffun,gofrr] = fit(ZenIntForLogistic,ZenProbForLogistic,HillEqn);
% ffun
% yfit=feval(ffun,0:0.0001:10); %Fitted function
% plot(0:0.0001:10,yfit,'Color',mygreen,'LineWidth',2);

xlim([0,0.5])
ylim([0,1])
xlabel('Integrated Mad')
ylabel('Probability')
hold off

%% Let's only plot later timepoints since we cannot stage early ones well

lagtime = 0;
ka = 1;

ProbForLogistic = [];
SogProbForLogistic = [];

MadForLogistic = [];
SogMadForLogistic = [];

figure 
hold on
fontsize(14,"points")

for x = 1:21
    % probability = FractionOn(:,x);
    % ProbForLogistic =[ProbForLogistic;probability];
    % instmad = [zeros(lagtime,1);BinnedMed(1:(end-lagtime),x)];
    % MadForLogistic = [MadForLogistic;instmad];
    % cmap = cbrewer2('Blues',30);
    % c=cmap(30-x,:);
    % scatter(instmad,probability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

    sogprobability = MutFractionOn(41:end,x);
    SogProbForLogistic =[SogProbForLogistic;sogprobability];
    soginstmad = BinnedMedSog(41:end,x);
    SogMadForLogistic = [SogMadForLogistic;soginstmad];
    cmap = cbrewer2('PuRd',30);
    c = cmap(25-x,:);
    scatter(soginstmad,sogprobability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

end
ylim([0,1])
xlabel('Instantaneous Mad')
ylabel('Probability')

startPoints = [1 1];
s = fitoptions('Method','NonlinearLeastSquares',... %
    'Lower',[0 0],...
    'Upper',[inf inf],...
    'Startpoint',startPoints);

HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);

% [ffun,gofrr] = fit(MadForLogistic,ProbForLogistic,HillEqn);
% ffun
% yfit=feval(ffun,0:0.0001:10); %Fitted function
% plot(0:0.0001:10,yfit,'Color',myblue,'LineWidth',2);

[ffun,gofrr] = fit(SogMadForLogistic,SogProbForLogistic,HillEqn);
ffun
gofrr
yfit=feval(ffun,0:0.0001:10); %Fitted function
plot(0:0.0001:10,yfit,'Color',mymagenta,'LineWidth',2);

xlim([0,0.05])
hold off

IntForLogistic = [];
SogIntForLogistic = [];
ZenIntForLogistic = [];

figure 
hold on
fontsize(14,"points")
for x = 1:21
    % probability = FractionOn(1:60,x);
    % instmad = [zeros(lagtime,1);MeanCumInt(1:(end-lagtime),x)];
    % IntForLogistic = [IntForLogistic;instmad];
    % cmap = cbrewer2('Blues',30);
    % c=cmap(30-x,:);
    % scatter(instmad./ka,probability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

    sogprobability = MutFractionOn(41:end,x);
    soginstmad = SogMeanCumInt(41:end,x);
    SogIntForLogistic = [SogIntForLogistic;soginstmad];
    cmap = cbrewer2('PuRd',30);
    c = cmap(30-x,:);
    scatter(soginstmad,sogprobability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

end


startPoints = [1 1];
s = fitoptions('Method','NonlinearLeastSquares',... %
    'Lower',[0 0],...
    'Upper',[inf inf],...
    'Startpoint',startPoints);

HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);

% [ffun,gofrr] = fit(IntForLogistic,ProbForLogistic,HillEqn);
% ffun
% yfit=feval(ffun,0:0.0001:10); %Fitted function
% plot(0:0.0001:10,yfit,'Color',myblue,'LineWidth',2);

[ffun,gofrr] = fit(SogIntForLogistic,SogProbForLogistic,HillEqn);
ffun
gofrr
yfit=feval(ffun,0:0.0001:10); %Fitted function
plot(0:0.0001:10,yfit,'Color',mymagenta,'LineWidth',2);

xlim([0,0.5])
ylim([0,1])
xlabel('Integrated Mad')
ylabel('Probability')
hold off

%%
startPoints = [1 1];
s = fitoptions('Method','NonlinearLeastSquares',... %
    'Lower',[0    0],...
    'Upper',[inf  inf],...
    'Startpoint',startPoints);
xdata = IntForLogistic;
ydata = ProbForLogistic;
% HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
% [ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
% %fitmax = round(max(xdata),2,'significant');
% yfit=feval(ffun,0:0.001:5); %Fitted function
%plot(0:0.001:5,yfit,'Color',myblue,'LineWidth',2);
% txt = [num2str(gofrr.rmse)];
% text(0.4,0.8,txt,'FontSize',14,'Color',myblue)
% txt = [num2str(ffun.a2)];
% text(0.6,0.8,txt,'FontSize',14,'Color',myblue)

xdata = SogIntForLogistic;
ydata = SogProbForLogistic;
% HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
% [ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
% %fitmax = round(max(xdata),2,'significant');
% yfit=feval(ffun,0:0.001:5); %Fitted function
%plot(0:0.001:5,yfit,'Color',mymagenta,'LineWidth',2);
% txt = [num2str(gofrr.rmse)];
% text(0.4,0.6,txt,'FontSize',14,'Color',mymagenta)
% txt = [num2str(ffun.a2)];
% text(0.6,0.6,txt,'FontSize',14,'Color',mymagenta)
% 
% txt = ['Int Time: ' num2str(inttime)];
% text(0,0,txt,'FontSize',14,'Color','black')
%xlim([0,0.05])
% ylabel('Probability')
% xlabel('Integrated Mad')
% hold off

figure
hold on
fontsize(14,"points")
startPoints = [1 1];
s = fitoptions('Method','NonlinearLeastSquares','Robust','LAR'... %
    'Lower',[0 0],...
    'Upper',[inf inf],...
    'Startpoint',startPoints);
xdata = [IntForLogistic]'./ka;
ydata = [ProbForLogistic]';

HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
% %Logistic = fittype('1/(1+exp(-(x-a1)/a2))');
[ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
ffun
yfit=feval(ffun,0:0.01:8); %Fitted function
%scatter(SogIntForLogistic,SogProbForLogistic,'MarkerEdgeColor',mymagenta,'MarkerFaceColor','none');
scatter(IntForLogistic./ka,ProbForLogistic,'MarkerEdgeColor',myblue,'MarkerFaceColor',myblue)
plot(0:0.01:8,yfit,'Color','black','LineWidth',2);
txt = ['Error: ' num2str(gofrr.rmse)];
text(0.05,0.1,txt,'FontSize',14)
xlim([0,0.4])
xlabel('Integrated Mad')
ylabel('Probability')
ylim([0,1])
hold off

figure
hold on
fontsize(14,"points")
startPoints = [1 1];
s = fitoptions('Method','NonlinearLeastSquares',... %
   'Lower',[0 0],...
   'Upper',[inf inf],...
   'Startpoint',startPoints);
xdata = [MadForLogistic]';
ydata = [ProbForLogistic]';
HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
%Logistic = fittype('1/(1+exp(-(x-a1)/a2))');
[ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
ffun
yfit=feval(ffun,0:0.001:0.04); %Fitted function
%scatter(SogMadForLogistic,SogProbForLogistic,'MarkerEdgeColor',mymagenta,'MarkerFaceColor','none');
scatter(MadForLogistic,ProbForLogistic,'MarkerEdgeColor',myblue,'MarkerFaceColor',myblue)
plot(0:0.001:0.04,yfit,'Color','black','LineWidth',2);
txt = ['Error: ' num2str(gofrr.rmse)];
text(0.005,0.1,txt,'FontSize',14)
%gofrr
xlim([0,0.04])
xlabel('Instantaneous Mad')
ylabel('Probability')
ylim([0,1])
hold off

%% What if it isn't the integral but the rate?

% ChangeMetric = BinnedMedRate;%./BinnedMed;
% SogChangeMetric = BinnedMedRateSog;%./BinnedMedSog;
% 
% lagtime = 0;
% 
% ProbForRateLogistic = [];
% SogProbForRateLogistic = [];
% 
% RateForLogistic = [];
% SogRateForLogistic = [];
% 
% figure
% hold on
% for x = 1:15
%     probability = FractionOn(21:60,x);
%     rate = [zeros(1,lagtime),ChangeMetric(21:(end-lagtime),x)];
%     probability = probability(rate(:,1)~=Inf,:);
%     rate = rate(rate(:,1)~=Inf,:);
%     probability = probability(rate(:,1)~=-Inf,:);
%     rate = rate(rate(:,1)~=-Inf,:);
%     ProbForRateLogistic =[ProbForRateLogistic;probability];
%     RateForLogistic = [RateForLogistic;rate];
% 
%     cmap = cbrewer2('Blues',30);
%     c=cmap(30-x,:);
%     scatter(rate,probability,'MarkerEdgeColor',c)
% 
%     sogprobability = MutFractionOn(21:60,x);
%     sograte = [zeros(1,lagtime),SogChangeMetric(21:(end-lagtime),x)];
%     sogprobability = sogprobability(sograte(:,1)~=Inf,:);
%     sograte = sograte(sograte(:,1)~=Inf,:);
%     sogprobability = sogprobability(sograte(:,1)~=-Inf,:);
%     sograte = sograte(sograte(:,1)~=-Inf,:);
%     SogProbForRateLogistic =[SogProbForRateLogistic;sogprobability];
%     SogRateForLogistic = [SogRateForLogistic;sograte];
%     cmap = cbrewer2('PuRd',30);
%     c = cmap(30-x,:);
%     scatter(sograte,sogprobability,'MarkerEdgeColor',c)
% 
% end
% ylabel('Probability')
% xlabel('dc/dt / c')
% hold off

% figure
% hold on
% xdata = RateForLogistic;
% ydata = ProbForRateLogistic;
% % HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
% % [ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
% % yfit=feval(ffun,0:0.001:5); %Fitted function
% scatter(xdata,ydata,'MarkerEdgeColor',myblue)
% 
% xdata = SogRateForLogistic;
% ydata = SogProbForRateLogistic;
% % HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
% % [ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
% % yfit=feval(ffun,0:0.001:5); %Fitted function
% scatter(xdata,ydata,'MarkerEdgeColor',mymagenta)
% hold off

