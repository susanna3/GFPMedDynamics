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


%% Fixed pMad Data

AllWTData = load("AllDataWT.mat", 'AllWTNucleiData');
AllWTTiming = load("AllTimingWT.mat", 'AllWTNucleiTiming');
AllWTData = AllWTData.AllWTNucleiData;
AllWTTiming = AllWTTiming.AllWTNucleiTiming;

AllMutData = load("AllDataSog.mat", 'AllMutantNucleiData');
AllMutTiming = load("AllTimingSog.mat", 'AllMutantNucleiTiming');
AllMutData = AllMutData.AllMutantNucleiData;
AllMutTiming = AllMutTiming.AllMutantNucleiTiming;

AllZenData = load("AllDataZen.mat", 'AllMutantNucleiData');
AllZenTiming = load("AllTimingZen.mat", 'AllMutantNucleiTiming');
AllZenData = AllZenData.AllMutantNucleiData;
AllZenTiming = AllZenTiming.AllMutantNucleiTiming;

cmap = cbrewer2('Blues',101);

NucByPosition = cell(1,21);
BinnedMed1 = zeros(21,60);
RateOfMed = zeros(21,60);
BinnedUpMed1= zeros(21,60);
BinnedLowMed1 = zeros(21,60);
counter = 0;
figure
hold on
ylim([0,0.6])
for d = 0:3:60
    counter = counter+1;
    PositionData = [];
    [datalen datawid] = size(AllWTData);
    for i = 1:datalen
        age = AllWTTiming(i,1);
        data = AllWTData{i};
        for n = 1:length(data)
            position = data(n,1);
            if position>=d &&  position<d+3
                PositionData = [PositionData; age, data(n,:)];
            else
                continue
            end
        end
    end
    time = PositionData(:,1);
    pmad = PositionData(:,3);
    timebins = discretize(time,[0,15,30:10:60],'IncludedEdge','right');
    Times = [0];
    MadLevels = [0];
    for t = 1:5
        binnedmed = pmad(timebins(:,1)==t,:);
        avgbinnedmed = mean(binnedmed);
        if avgbinnedmed > 0
            Times = [Times;t*10];
            MadLevels = [MadLevels;avgbinnedmed];
        else
            continue
        end
    end
    [ft,gof] = fit(time,pmad,'SmoothingSpline','SmoothingParam',0.00005);
    ci = gof.rmse^2;
    yfitted = feval(ft,1:60);
    rateofmed = differentiate(ft,1:60);
    yfup = yfitted+ci;
    yfdown = yfitted-ci;
    BinnedMed1(counter,:) = yfitted;
    BinnedUpMed1(counter,:) = yfup;
    BinnedLowMed1(counter,:) = yfdown;
    RateOfMed(counter,:) = rateofmed;
    c = cmap(101-d,:);
    plot(1:60,yfitted,'Color',c)
    scatter(time,pmad,1,'MarkerEdgeColor',c,'MarkerFaceColor','none')
    plot(1:60,yfup,'Color','r','LineStyle','--')
    plot(1:60,yfdown,'Color','r','LineStyle','--')
    NucByPosition{counter} = PositionData;
end
hold off

% Now let's do the same for sog
MutNucByPosition = cell(1,21);
BinnedMedSog1 = zeros(21,60);
BinnedUpMedSog1= zeros(21,60);
BinnedLowMedSog1 = zeros(21,60);
SogRateOfMed = zeros(21,60);
counter = 0;
figure
hold on
ylim([0,0.6])
for d = 0:3:60
    counter = counter+1;
    PositionData = [];
    [datalen datawid] = size(AllMutData);
    for i = 1:datalen
        age = AllMutTiming(i,1);
        data = AllMutData{i};
        for n = 1:length(data)
            position = data(n,1);
            if position>=d &&  position<d+3
                PositionData = [PositionData; age, data(n,:)];
            else
                continue
            end
        end
    end
    time = PositionData(:,1);
    pmad = PositionData(:,3);

    timebins = discretize(time,[0,15,30:10:60],'IncludedEdge','right');
    Times = [0];
    MadLevels = [0];
    for t = 1:5
        binnedmed = pmad(timebins(:,1)==t,:);
        avgbinnedmed = mean(binnedmed);
        if avgbinnedmed > 0
            Times = [Times;t*10];
            MadLevels = [MadLevels;avgbinnedmed];
        else
            continue
        end
    end
    [ft,gof] = fit(time,pmad,'SmoothingSpline','SmoothingParam',0.0001);
    ci = gof.rmse^2;
    yfitted = feval(ft,1:60);
    rateofmed = differentiate(ft,1:60);
    yfup = yfitted+ci;
    yfdown = yfitted-ci;
    BinnedMedSog1(counter,:) = yfitted;
    BinnedUpMedSog1(counter,:) = yfup;
    BinnedLowMedSog1(counter,:) = yfdown;
    SogRateOfMed(counter,:) = rateofmed;
    c = cmap(101-d,:);
    plot(1:60,yfitted,'Color',c)
    scatter(time,pmad,0.25,'MarkerEdgeColor',c,'MarkerFaceColor','none')
    plot(1:60,yfup,'Color','r','LineStyle','--')
    plot(1:60,yfdown,'Color','r','LineStyle','--')
    MutNucByPosition{counter} = PositionData;
end
hold off

% Now let's do the same for zen
ZenNucByPosition = cell(1,21);
BinnedMedZen1 = zeros(21,60);
BinnedUpMedZen1= zeros(21,60);
BinnedLowMedZen1 = zeros(21,60);
ZenRateOfMed = zeros(21,60);
counter = 0;
figure
hold on
ylim([0,0.6])
for d = 0:3:60
    counter = counter+1;
    PositionData = [];
    [datalen datawid] = size(AllZenData);
    for i = 1:datalen
        age = AllZenTiming(i,1);
        data = AllZenData{i};
        for n = 1:length(data)
            position = data(n,1);
            if position>=d &&  position<d+3
                PositionData = [PositionData; age, data(n,:)];
            else
                continue
            end
        end
    end
    time = PositionData(:,1);
    pmad = PositionData(:,3);
    timebins = discretize(time,[0,15,30:10:60],'IncludedEdge','right');
    Times = [0];
    MadLevels = [0];
    for t = 1:5
        binnedmed = pmad(timebins(:,1)==t,:);
        avgbinnedmed = mean(binnedmed);
        if avgbinnedmed > 0
            Times = [Times;t*10];
            MadLevels = [MadLevels;avgbinnedmed];
        else
            continue
        end
    end
    [ft,gof] = fit(time,pmad,'SmoothingSpline','SmoothingParam',0.0001);
    ci = gof.rmse^2;
    yfitted = feval(ft,1:60);
    rateofmed = differentiate(ft,1:60);
    yfup = yfitted+ci;
    yfdown = yfitted-ci;
    BinnedMedZen1(counter,:) = yfitted;
    BinnedUpMedZen1(counter,:) = yfup;
    BinnedLowMedZen1(counter,:) = yfdown;
    ZenRateOfMed(counter,:) = rateofmed;
    c = cmap(101-d,:);
    plot(1:60,yfitted,'Color',c)
    scatter(time,pmad,0.25,'MarkerEdgeColor',c,'MarkerFaceColor','none')
    plot(1:60,yfup,'Color','r','LineStyle','--')
    plot(1:60,yfdown,'Color','r','LineStyle','--')
    ZenNucByPosition{counter} = PositionData;
end
hold off

%% Remove background from first 20min

BinnedMed = max(BinnedMed1,0);
BinnedMedSog = max(BinnedMedSog1,0);
BinnedMedZen = max(BinnedMedZen1,0);

BinnedMed = BinnedMed1-mean(BinnedMed1(:,1:20),2);
BinnedMedSog = BinnedMedSog1-mean(BinnedMedSog1(:,1:20),2);
BinnedMedZen = BinnedMedZen1-mean(BinnedMedZen1(:,1:20),2);

BinnedMed = max(BinnedMed,0);
BinnedMedSog = max(BinnedMedSog,0);
BinnedMedZen = max(BinnedMedZen,0);

%% Rescale to live data value

scale_factor = 0.1082;

BinnedMedZen = BinnedMedZen*scale_factor;

%% Get width and amplitude

WidthOverTime = [];
AmplitudeOverTime = [];
SogWidthOverTime = [];
SogAmplitudeOverTime = [];
ZenWidthOverTime = [];
ZenAmplitudeOverTime = [];
figure
hold on
for t = 10:60
    timepointdata = BinnedMed(:,t);
    amplitude = max(timepointdata);
    ft = fit([0:5:100]',timepointdata,'SmoothingSpline','SmoothingParam',0.005);
    plot(1:100,feval(ft,1:100))
    yfitted = feval(ft,1:100);
    halfMax = amplitude/2;
    %width = trapz(timepointdata)/amplitude;
    width = find(yfitted >= halfMax, 1, 'last');
    WidthOverTime = [WidthOverTime,width*2];
    AmplitudeOverTime = [AmplitudeOverTime,amplitude];

    sogtimepointdata = BinnedMedSog(:,t);
    sogamplitude = max(sogtimepointdata(:,1));
    sogft = fit([0:5:100]',sogtimepointdata,'SmoothingSpline','SmoothingParam',0.005);
    yfitted = feval(sogft,1:100);
    halfMax = sogamplitude/2;
    %width = trapz(timepointdata)/amplitude;
    sogwidth = find(yfitted >= halfMax, 1, 'last');
    SogWidthOverTime = [SogWidthOverTime,sogwidth*2];
    SogAmplitudeOverTime = [SogAmplitudeOverTime,sogamplitude];

    zentimepointdata = BinnedMedZen(:,t);
    zenamplitude = max(zentimepointdata(1:10,1));
    zenft = fit([0:5:100]',zentimepointdata,'SmoothingSpline','SmoothingParam',0.005);
    yfitted = feval(zenft,1:100);
    halfMax = zenamplitude/2;
    %width = trapz(timepointdata)/amplitude;
    zenwidth = find(yfitted >= halfMax, 1, 'last');
    ZenWidthOverTime = [ZenWidthOverTime,zenwidth*2];
    ZenAmplitudeOverTime = [ZenAmplitudeOverTime,zenamplitude];
end
hold off

figure
hold on
fontsize(14,'points')
plot(10:60,WidthOverTime,'Color',myblue,'LineWidth',2)
plot(10:60,SogWidthOverTime,'Color',mymagenta,'LineWidth',2)
plot(10:60,ZenWidthOverTime,'Color',mygreen,'LineWidth',2)
xlabel('Time (Minutes)')
ylabel('Width')
hold off

figure
hold on
fontsize(14,'points')
plot(10:60,AmplitudeOverTime,'Color',myblue,'LineWidth',2)
plot(10:60,SogAmplitudeOverTime,'Color',mymagenta,'LineWidth',2)
plot(10:60,ZenAmplitudeOverTime,'Color',mygreen,'LineWidth',2)
xlabel('Time (Minutes)')
ylabel('Amplitude')
hold off



%% Plot curves and derivative

cmap = cbrewer2('Blues',25);
figure
hold on
fontsize(14,"points")
for x = 1:21
    xdata = 1:60;
    ydata = BinnedMed(x,:);
    c = cmap(25-x,:);
    %ydataup = BinnedUpMed1(x,:);
    %ydatadown = BinnedLowMed1(x,:);
    %fill([xdata, flip(xdata)], [ydataup, flip(ydatadown)], c,'FaceAlpha',0.5, 'EdgeColor','none')
    plot(xdata,ydata,'Color',c,'LineWidth',1)
end
xlabel('Time (Minutes)')
ylabel('Nuclear pMad')
ylim([0,0.75])
hold off

cmap = cbrewer2('PuRd',30);
figure
hold on
fontsize(14,"points")
for x = 1:21
    xdata = 1:60;
    ydata = BinnedMedSog(x,:);
    c = cmap(30-x,:);
    %ydataup = BinnedUpMedSog1(x,:);
    %ydatadown = BinnedLowMedSog1(x,:);
    %fill([xdata, flip(xdata)], [ydataup, flip(ydatadown)], c,'FaceAlpha',0.5, 'EdgeColor','none')
    plot(xdata,ydata,'Color', c,'LineWidth',1)
end
xlabel('Time (Minutes)')
ylabel('Nuclear pMad')
ylim([0,0.75])
hold off

cmap = cbrewer2('Greens',22);
figure
hold on
fontsize(14,"points")
for x = 1:21
    xdata = 1:60;
    ydata = BinnedMedZen(x,:);
    c = cmap(22-x,:);
    %ydataup = BinnedUpMedZen(x,:);
    %ydatadown = BinnedLowMedZen1(x,:);
    %fill([xdata, flip(xdata)], [ydataup, flip(ydatadown)], c,'FaceAlpha',0.5, 'EdgeColor','none')
    plot(xdata,ydata,'Color',c,'LineWidth',1)
end
xlabel('Time (Minutes)')
ylabel('Nuclear pMad')
ylim([0,0.75])
hold off

cmap = cbrewer2('Blues',25);
figure
hold on
fontsize(14,"points")
for x = 1:21
    xdata = 1:60;
    ydata = RateOfMed(x,:);
    c = cmap(25-x,:);
    plot(xdata,ydata,'Color',c,'LineWidth',1.5,'LineStyle','--')
end
xlabel('Time (Minutes)')
ylabel('dMad/dt')
ylim([0,0.03])
hold off

cmap = cbrewer2('Greens',25);
figure
hold on
fontsize(14,"points")
for x = 1:21
    xdata = 1:60;
    ydata = SogRateOfMed(x,:);
    c = cmap(25-x,:);
    plot(xdata,ydata,'Color',c,'LineWidth',1.5,'LineStyle','--')
end
xlabel('Time (Minutes)')
ylabel('dMad/dt')
ylim([0,0.03])
hold off

cmap = cbrewer2('PuRd',25);
figure
hold on
fontsize(14,"points")
for x = 1:21
    xdata = 1:60;
    ydata = ZenRateOfMed(x,:);
    c = cmap(25-x,:);
    plot(xdata,ydata,'Color',c,'LineWidth',1.5,'LineStyle','--')
end
xlabel('Time (Minutes)')
ylabel('dMad/dt')
%ylim([0,0.03])
hold off

%% Finally, I need to restructure the data for my gene of interest as well

% Load in data

AllWTDataUsh = load("AllDataWTUsh.mat", 'AllWTNucleiData');
AllWTTimingUsh = load("AllTimingWTUsh.mat", 'AllWTNucleiTiming');
AllWTDataUsh = AllWTDataUsh.AllWTNucleiData;
AllWTTimingUsh = AllWTTimingUsh.AllWTNucleiTiming;

AllMutDataUsh = load("AllDataSogUsh.mat", 'AllMutantNucleiData');
AllMutTimingUsh = load("AllTimingSogUsh.mat", 'AllMutantNucleiTiming');
AllMutDataUsh = AllMutDataUsh.AllMutantNucleiData;
AllMutTimingUsh = AllMutTimingUsh.AllMutantNucleiTiming;

AllZenDataUsh = load("AllDataZenUsh.mat", 'AllMutantNucleiData');
AllZenTimingUsh = load("AllTimingZenUsh.mat", 'AllMutantNucleiTiming');
AllZenDataUsh = AllZenDataUsh.AllMutantNucleiData;
AllZenTimingUsh = AllZenTimingUsh.AllMutantNucleiTiming;

NucByPositionUsh = cell(1,21);
counter = 0;
for d = 0:3:60
    counter = counter+1;
    PositionData = [];
    for i = 1:length(AllWTDataUsh)
        age = AllWTTimingUsh(i,1);
        data = AllWTDataUsh{i};
        for n = 1:length(data)
            position = data(n,1);
            if position>=d &&  position<d+3
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
for d = 0:3:60
    counter = counter+1;
    PositionData = [];
    for i = 1:length(AllMutDataUsh)
        age = AllMutTimingUsh(i,1);
        data = AllMutDataUsh{i};
        for n = 1:length(data)
            position = data(n,1);
            if position>=d &&  position<d+3
                PositionData = [PositionData; age, data(n,:)];
            else
                continue
            end
        end
    end
    MutNucByPositionUsh{counter} = PositionData;
end

ZenNucByPositionUsh = cell(1,21);
counter = 0;
for d = 0:3:60
    counter = counter+1;
    PositionData = [];
    for i = 1:length(AllZenDataUsh)
        age = AllZenTimingUsh(i,1);
        data = AllZenDataUsh{i};
        for n = 1:length(data)
            position = data(n,1);
            if position>=d &&  position<d+3
                PositionData = [PositionData; age, data(n,:)];
            else
                continue
            end
        end
    end
    ZenNucByPositionUsh{counter} = PositionData;
end


%% Find the fraction of nuclei on over time in each position bin

cmap=cbrewer2('Blues',25);

FractionOn = [];
figure
hold on
for x = 1:21
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
    [logitCoef,dev] = glmfit(collectedtimes,collectedfracs,'binomial','logit');
    logitFit = glmval(logitCoef,1:60,'logit');
    logitFit(logitFit(:,1)>0.9999,1)=1;
    FractionOn = [FractionOn,logitFit];
    c = cmap(25-x,:);
    scatter(collectedtimes,collectedfracs,'MarkerEdgeColor',c,'MarkerFaceColor','none')
    plot(1:60,logitFit,'Color',c)
    ylim([0 1])
end
hold off

MutFractionOn = [];
% figure
% hold on
for x = 1:21
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
    [logitCoef,dev] = glmfit(collectedtimes,collectedfracs,'binomial','logit');
    logitFit = glmval(logitCoef,1:60,'logit');
    logitFit(logitFit(:,1)>0.9999,1)=1;
    MutFractionOn = [MutFractionOn,logitFit];
    c = cmap(25-x,:);
    % scatter(collectedtimes,collectedfracs,'MarkerEdgeColor',c,'MarkerFaceColor','none')
    % plot(1:60,logitFit,'Color',c)
end
% hold off

ZenFractionOn = [];
% figure
% hold on
for x = 1:21
    nucatpos = ZenNucByPositionUsh{x};
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
    [logitCoef,dev] = glmfit(collectedtimes,collectedfracs,'binomial','logit');
    logitFit = glmval(logitCoef,1:60,'logit');
    logitFit(logitFit(:,1)>0.9999,1)=1;
    ZenFractionOn = [ZenFractionOn,logitFit];
    c = cmap(25-x,:);
    % scatter(collectedtimes,collectedfracs,'MarkerEdgeColor',c,'MarkerFaceColor','none')
    % plot(1:60,logitFit,'Color',c)
end
% hold off

FractionOn(isnan(FractionOn))=0;
MutFractionOn(isnan(MutFractionOn))=0;
ZenFractionOn(isnan(ZenFractionOn))=0;

figure
h = heatmap(FractionOn,'Colormap',mybluemap,'ColorLimits',[0,1]);
h.YDisplayLabels = {'1','','','','','','','','',...
    '10','','','','','','','','','',...
    '20','','','','','','','','','',...
    '30','','','','','','','','','',...
    '40','','','','','','','','','',...
    '50','','','','','','','','','',...
    '60'};
h.XDisplayLabels = {'0','3','6','9','12','15','18','21','24','27','30','33','36','39','42','45','48','51','54','57','60'};
h.FontSize = 14;
h.GridVisible = 'off';
ylabel('Minutes');
xlabel('Microns');

figure
h = heatmap(MutFractionOn,'Colormap',mygreenmap,'ColorLimits',[0,1]);
h.YDisplayLabels = {'1','','','','','','','','',...
    '10','','','','','','','','','',...
    '20','','','','','','','','','',...
    '30','','','','','','','','','',...
    '40','','','','','','','','','',...
    '50','','','','','','','','','',...
    '60'};
h.XDisplayLabels = {'0','3','6','9','12','15','18','21','24','27','30','33','36','39','42','45','48','51','54','57','60'};
h.FontSize = 14;
h.GridVisible = 'off';
ylabel('Minutes');
xlabel('Microns');

figure
h = heatmap(ZenFractionOn,'Colormap',mymagentamap,'ColorLimits',[0,1]);
h.YDisplayLabels = {'1','','','','','','','','',...
    '10','','','','','','','','','',...
    '20','','','','','','','','','',...
    '30','','','','','','','','','',...
    '40','','','','','','','','','',...
    '50','','','','','','','','','',...
    '60'};
h.XDisplayLabels = {'0','3','6','9','12','15','18','21','24','27','30','33','36','39','42','45','48','51','54','57','60'};
h.FontSize = 14;
h.GridVisible = 'off';
ylabel('Minutes');
xlabel('Microns');

%% Next we identify the time on, which is when 50% of the max of nuclei in the bin are on
% Find on time at target percentile

% targetp = max(max(FractionOn))*0.5;
% 
% TimeOn = zeros(1,60);
% 
% ProbOnOverTime = [];
% for x = 1:21
%     onovertime = FractionOn(:,x);
%     ProbOnOverTime = [ProbOnOverTime,onovertime];
% 
%     for t = 1:60
%         isiton = onovertime(t,1);
%         if isiton > targetp
%             TimeOn(1,x) = t;
%             break
%         else
%             continue
%         end
%     end
% end
% 
% MutTimeOn = zeros(1,60);
% 
% MutProbOnOverTime = [];
% for x = 1:21
%     onovertime = MutFractionOn(:,x);
%     MutProbOnOverTime = [MutProbOnOverTime,onovertime];
% 
%     for t = 1:60
%         isiton = onovertime(t,1);
%         if isiton > targetp
%             MutTimeOn(1,x) = t;
%             break
%         else
%             continue
%         end
%     end
% end
% 
% ZenTimeOn = zeros(1,60);
% 
% ZenProbOnOverTime = [];
% for x = 1:21
%     onovertime = ZenFractionOn(:,x);
%     ZenProbOnOverTime = [ZenProbOnOverTime,onovertime];
% 
%     for t = 1:60
%         isiton = onovertime(t,1);
%         if isiton > targetp
%             ZenTimeOn(1,x) = t;
%             break
%         else
%             continue
%         end
%     end
% end
% 
% 
% %% Now let's find the integration window
% 
% MadIntAtOnPerBin = zeros(60,21);
% for x = 1:21
%     madintperbin = zeros(60,1);
%     nucatpos = MutNucByPosition{x};
%     timeon=TimeOn(1,x);
%     xdata = nucatpos(:,1);
%     ydata = nucatpos(:,3);
%     [f gof] = fit(xdata,ydata,'SmoothingSpline','SmoothingParam',0.005);
%     for t = 1:timeon-1
%         timeint = integrate(f,timeon,timeon-t);
%         madintperbin(t:end,1)=timeint;
%     end
%     MadIntAtOnPerBin(:,x)=madintperbin;
% end
% 
% SlopePerTau = [];
% figure
% hold on
% for inttime = 1:60
%     ydata = nonzeros(MadIntAtOnPerBin(inttime,:));
%     xdata = 1:length(ydata);
%     plot(xdata,ydata)
%     slope = polyfit(xdata,ydata,1);
%     coef2 = polyfit(xdata,ydata,1);
%     y2 = polyval(coef2,xdata);
%     SlopePerTau = [SlopePerTau;slope(:,1)./inttime];
% end
% hold off
% 
% counter = 0;
% for t = 1:60
%     slope = SlopePerTau(t,1);
%     if slope(:,1) < 0
%         counter = counter+1;
%         optinttime = 60;
%         continue
%     else
%         optinttime = counter+1;
%         break
%     end
% end
% % 
% figure
% hold on
% scatter((1:60),SlopePerTau,'MarkerFaceColor','red','MarkerEdgeColor','red')
% plot((0:60),zeros(1,61),'--k')
% %ylim([-0.01,0.005])
% hold off

%% Now I need to find the integral range in WT and Mutant data


MeanCumInt = [];
inttime = 40; %optinttime;
offtime = 0;
for x = 1:21
    avgdata = BinnedMed(x,:)';
    %avgdata = avgdata.^2;
    int = [zeros(offtime,1);cumtrapz(avgdata(offtime+1:inttime,1))];
    for t = (inttime+1):60
        intwindow = t-inttime;
        timeint = max(cumtrapz(avgdata(intwindow:t,1)));
        int = [int;timeint];
    end
    MeanCumInt = [MeanCumInt,int];
end


SogMeanCumInt = [];
%inttime = 10; %optinttime;
offtime = 0;
for x = 1:21
    avgdata = BinnedMedSog(x,:)';
    %avgdata = avgdata.^2;
    int = [zeros(offtime,1);cumtrapz(avgdata(offtime+1:inttime,1))];
    for t = (inttime+1):60
        intwindow = t-inttime;
        timeint = max(cumtrapz(avgdata(intwindow:t,1)));
        int = [int;timeint];
    end
    SogMeanCumInt = [SogMeanCumInt,int];
end

ZenMeanCumInt = [];
%inttime = 10; %optinttime;
offtime = 0;
for x = 1:21
    avgdata = BinnedMedZen(x,:)';
    %avgdata=avgdata.^2;
    int = [zeros(offtime,1);cumtrapz(avgdata(offtime+1:inttime,1))];
    for t = (inttime+1):60
        intwindow = t-inttime;
        timeint = max(cumtrapz(avgdata(intwindow:t,1)));
        int = [int;timeint];
    end
    ZenMeanCumInt = [ZenMeanCumInt,int];
end


%% Now I need to find the level of integrated Mad that relates to On Time for each position in WT

% IntegralRange = [];
% InstMad = [];
% for t = 1:21
%     timeon = TimeOn(1,t);
%     if timeon > 0
%         upperint = UpperCumInt(t,timeon-2);
%         lowerint = LowerCumInt (t,timeon-2);
%         IntegralRange = [IntegralRange;upperint,lowerint];
%         InstMad = [InstMad;BinnedUpMed(t,timeon-2),BinnedLowMed(t,timeon-2)];
%     else
%         continue
%     end
% end
% 
% % What about this range in sog mutants? Do they fall on top of each other?
% SogIntegralRange = [];
% SogInstMad = [];
% for t = 1:21
%     timeon = MutTimeOn(1,t);
%     if timeon > 0
%         upperint = SogUpperCumInt(t,timeon-2);
%         lowerint = SogLowerCumInt (t,timeon-2);
%         SogIntegralRange = [SogIntegralRange;upperint,lowerint];
%         SogInstMad = [SogInstMad;BinnedUpMedSog(t,timeon-2),BinnedLowMedSog(t,timeon-2)];
%     else
%         continue
%     end
% end
% 
% AllIntegralRange = [IntegralRange;SogIntegralRange];
% maxint = max(IntegralRange(:,1));
% minint = min(IntegralRange(:,2));
% 
% sogmaxint = mean(SogIntegralRange(:,1));
% sogminint = mean(SogIntegralRange(:,2));
% 
% figure
% hold on
% x = length(IntegralRange);
% x2 = length(SogIntegralRange);
% plot(1:x,IntegralRange);
% plot(1:x2,SogIntegralRange);
% ylim([0,12])
% hold off
% 
% figure
% hold on
% for i = 1:21
%     x = 1:60;
%     y1 = UpperCumInt(i,:);
%     y2 = LowerCumInt(i,:);
%     y = mean([y1;y2],1);
%     plot(x,y)
% end
% ylim([0,10])
% hold off
% 
% figure
% hold on
% for i = 1:21
%     x = 1:60;
%     y1 = SogUpperCumInt(i,:);
%     y2 = SogLowerCumInt(i,:);
%     y = mean([y1;y2],1);
%     plot(x,y)
% end
% ylim([0,10])
% hold off
% 
% f=figure();
% hold on
% fontsize(gcf,14,'points');
% x = 0:3:60;
% [lendata wid] = size(nonzeros(IntegralRange(:,1)));
% x = x(1,1:lendata);
% fill([x, flip(x)], [nonzeros(IntegralRange(:,1))', flip(nonzeros(IntegralRange(:,2))')], [0.85, 0.11, 0.38],'FaceAlpha',0.2, 'EdgeColor','none')
% scatter(x,mean(IntegralRange,2),'MarkerFaceColor',mymagenta,'MarkerEdgeColor',mymagenta)
% 
% x = 0:3:60;
% [lendata wid] = size(nonzeros(SogIntegralRange(:,1)));
% x = x(1,1:lendata);
% fill([x, flip(x)], [nonzeros(SogIntegralRange(:,1))', flip(nonzeros(SogIntegralRange(:,2))')], [0, 0.3, 0.25],'FaceAlpha',0.2, 'EdgeColor','none')
% scatter(x,nonzeros(mean(SogIntegralRange,2)),'MarkerFaceColor',mygreen,'MarkerEdgeColor',mygreen)
% %ylim([0,4])
% ylabel('Integrated Mad')
% xlabel('Microns')
% hold off
% 
% f=figure();
% hold on
% fontsize(gcf,14,'points');
% x = 0:3:60;
% [lendata wid] = size(nonzeros(InstMad(:,1)));
% x = x(1,1:lendata);
% fill([x, flip(x)], [nonzeros(InstMad(:,1))', flip(nonzeros(InstMad(:,2))')], [0.85, 0.11, 0.38],'FaceAlpha',0.2, 'EdgeColor','none')
% scatter(x,nonzeros(mean(InstMad,2)),'MarkerFaceColor',mymagenta,'MarkerEdgeColor',mymagenta)
% 
% x = 0:3:60;
% [lendata wid] = size(nonzeros(SogInstMad(:,1)));
% x = x(1,1:lendata);
% fill([x, flip(x)], [nonzeros(SogInstMad(:,1))', flip(nonzeros(SogInstMad(:,2))')], [0, 0.3, 0.25],'FaceAlpha',0.2, 'EdgeColor','none')
% scatter(x,nonzeros(mean(SogInstMad,2)),'MarkerFaceColor',mygreen,'MarkerEdgeColor',mygreen)
% %ylim([0,0.4])
% ylabel('Instantaneous Mad')
% xlabel('Microns')
% hold off

%% Express probability as probability of turning on if you haven't already

ProbabilityOn = [];

cmap = cbrewer2('Blues',25);
figure
hold on
for x = 1:21
    position = [0;round(FractionOn(:,x).*100)];
    totaloff = 100;
    probability = zeros(60,1);
    for t = 2:61
        prevon = position(t-1,1);
        currentlyon = position(t,1);
        newlyon = currentlyon-prevon;
        if newlyon >0 
            prob = newlyon/totaloff;
            probability(t:end,1) = prob;
            totaloff = totaloff-newlyon;
        else
            continue
        end
    end
    xdata = (1:60)';
    ydata = probability;
    xdata(isnan(probability)) = [];
    ydata(isnan(probability)) = [];
    c = cmap(25-x,:);
    scatter(xdata,ydata,'MarkerEdgeColor',c);
    [logitCoef,dev] = glmfit(xdata,ydata,'binomial','logit');
    logitFit = glmval(logitCoef,1:60,'logit');
    logitFit(logitFit(:,1)>0.9999,1)=0;
    ProbabilityOn = [ProbabilityOn,logitFit];
    plot(1:60,logitFit,'Color',c)
end
%ylim([0,0.5])
hold off


SogProbabilityOn = [];

cmap = cbrewer2('Blues',25);
figure
hold on
for x = 1:21
    position = [0;round(MutFractionOn(:,x).*100)];
    totaloff = 100;
    probability = zeros(60,1);
    for t = 2:61
        prevon = position(t-1,1);
        currentlyon = position(t,1);
        newlyon = currentlyon-prevon;
        if newlyon >0 
            prob = newlyon/totaloff;
            probability(t:end,1) = prob;
            totaloff = totaloff-newlyon;
        else
            continue
        end
    end
    xdata = (1:60)';
    ydata = probability;
    xdata(isnan(probability)) = [];
    ydata(isnan(probability)) = [];
    c = cmap(25-x,:);
    scatter(xdata,ydata,'MarkerEdgeColor',c);
    [logitCoef,dev] = glmfit(xdata,ydata,'binomial','logit');
    logitFit = glmval(logitCoef,1:60,'logit');
    logitFit(logitFit(:,1)>0.9999,1)=0;
    SogProbabilityOn = [SogProbabilityOn,logitFit];
    plot(1:60,logitFit,'Color',c)
end
%ylim([0,0.5])
hold off


ZenProbabilityOn = [];

cmap = cbrewer2('Blues',25);
figure
hold on
for x = 1:21
    position = [0;round(ZenFractionOn(:,x).*100)];
    totaloff = 100;
    probability = zeros(60,1);
    for t = 2:61
        prevon = position(t-1,1);
        currentlyon = position(t,1);
        newlyon = currentlyon-prevon;
        if newlyon >0 
            prob = newlyon/totaloff;
            probability(t:end,1) = prob;
            totaloff = totaloff-newlyon;
        else
            continue
        end
    end
    xdata = (1:60)';
    ydata = probability;
    xdata(isnan(probability)) = [];
    ydata(isnan(probability)) = [];
    c = cmap(25-x,:);
    scatter(xdata,ydata,'MarkerEdgeColor',c);
    [logitCoef,dev] = glmfit(xdata,ydata,'binomial','logit');
    logitFit = glmval(logitCoef,1:60,'logit');
    logitFit(logitFit(:,1)>0.9999,1)=0;
    ZenProbabilityOn = [ZenProbabilityOn,logitFit];
    plot(1:60,logitFit,'Color',c)
end
%ylim([0,0.5])
hold off


%% Plot probability versus integrated or instantaneous Mad for each position


ProbForLogistic = [];
SogProbForLogistic = [];
ZenProbForLogistic = [];

MadForLogistic = [];
SogMadForLogistic = [];
ZenMadForLogistic = [];

lagtime = 0;

figure 
hold on
fontsize(14,"points")
for x = 1:21     
    probability = FractionOn(:,22-x);
    ProbForLogistic =[ProbForLogistic;probability];
    instmad = [zeros(1,lagtime),BinnedMed(22-x,1:(end-lagtime))];
    MadForLogistic = [MadForLogistic;instmad'];
    cmap = cbrewer2('Blues',30);
    c=cmap(x+5,:);
    scatter(instmad,probability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

    sogprobability = MutFractionOn(:,22-x);
    SogProbForLogistic =[SogProbForLogistic;sogprobability];
    soginstmad = [zeros(1,lagtime),BinnedMedSog(22-x,1:(end-lagtime))];
    SogMadForLogistic = [SogMadForLogistic;soginstmad'];
    cmap = cbrewer2('PuRd',30);
    c = cmap(x+5,:);
    scatter(soginstmad,sogprobability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

    zenprobability = ZenFractionOn(:,22-x);
    ZenProbForLogistic =[ZenProbForLogistic;zenprobability];
    zeninstmad = [zeros(1,lagtime),BinnedMedZen(22-x,1:(end-lagtime))];
    ZenMadForLogistic = [ZenMadForLogistic;zeninstmad'];
    cmap = cbrewer2('Greens',30);
    c = cmap(x+5,:);
    scatter(zeninstmad,zenprobability,'MarkerEdgeColor',c,'MarkerFaceColor',c)
end
startPoints = [1 1];
s = fitoptions('Method','NonlinearLeastSquares',... %
    'Lower',[0    0],...
    'Upper',[inf  inf],...
    'Startpoint',startPoints);
xdata = MadForLogistic;
ydata = ProbForLogistic;
HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
[ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
%fitmax = round(max(xdata),2,'significant');
yfit=feval(ffun,0:0.001:5); %Fitted function
plot(0:0.001:5,yfit,'Color',myblue,'LineWidth',2);
txt = [num2str(gofrr.rmse)];
text(0.4,0.8,txt,'FontSize',14,'Color',myblue)

xdata = SogMadForLogistic;
ydata = SogProbForLogistic;
HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
[ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
%fitmax = round(max(xdata),2,'significant');
yfit=feval(ffun,0:0.001:5); %Fitted function
plot(0:0.001:5,yfit,'Color',mymagenta,'LineWidth',2);
txt = [num2str(gofrr.rmse)];
text(0.4,0.6,txt,'FontSize',14,'Color',mymagenta)

xdata = ZenMadForLogistic;
ydata = ZenProbForLogistic;
HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
[ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
%fitmax = round(max(xdata),2,'significant');
yfit=feval(ffun,0:0.001:5); %Fitted function
plot(0:0.001:5,yfit,'Color',mygreen,'LineWidth',2);
txt = [num2str(gofrr.rmse)];
text(0.4,0.4,txt,'FontSize',14,'Color',mygreen)
xlim([0,0.5])
ylabel('Probability')
xlabel('Instantaneous Mad')
hold off
xlim([0 0.1])
hold off

IntForLogistic = [];
SogIntForLogistic = [];
ZenIntForLogistic = [];
figure 
hold on
fontsize(14,"points")
for x = 1:21
    probability = FractionOn(:,22-x);
    instmad = [zeros(1,lagtime),MeanCumInt(1:(end-lagtime),22-x)];
    IntForLogistic = [IntForLogistic;instmad];
    cmap = cbrewer2('Blues',30);
    c=cmap(x+5,:);
    scatter(instmad,probability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

    sogprobability = MutFractionOn(:,22-x);
    soginstmad = [zeros(1,lagtime),SogMeanCumInt(1:(end-lagtime),22-x)];
    SogIntForLogistic = [SogIntForLogistic;soginstmad];
    cmap = cbrewer2('PuRd',30);
    c = cmap(x+5,:);
    scatter(soginstmad,sogprobability,'MarkerEdgeColor',c,'MarkerFaceColor',c)

    zenprobability = ZenFractionOn(:,22-x);
    zeninstmad = [zeros(1,lagtime),ZenMeanCumInt(1:(end-lagtime),22-x)];
    ZenIntForLogistic = [ZenIntForLogistic;zeninstmad];
    cmap = cbrewer2('Greens',30);
    c = cmap(x+5,:);
    scatter(zeninstmad,zenprobability,'MarkerEdgeColor',c,'MarkerFaceColor',c)
end

startPoints = [1 1];
s = fitoptions('Method','NonlinearLeastSquares',... %
    'Lower',[0    0],...
    'Upper',[inf  inf],...
    'Startpoint',startPoints);
xdata = IntForLogistic;
ydata = ProbForLogistic;
HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
[ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
%fitmax = round(max(xdata),2,'significant');
yfit=feval(ffun,0:0.001:50); %Fitted function
plot(0:0.001:50,yfit,'Color',myblue,'LineWidth',2);
txt = [num2str(ffun.a1)];
text(0.4,0.8,txt,'FontSize',14,'Color',myblue)
txt = [num2str(ffun.a2)];
text(0.6,0.8,txt,'FontSize',14,'Color',myblue)

xdata = SogIntForLogistic;
ydata = SogProbForLogistic;
HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
[ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
%fitmax = round(max(xdata),2,'significant');
yfit=feval(ffun,0:0.001:50); %Fitted function
plot(0:0.001:50,yfit,'Color',mymagenta,'LineWidth',2);
txt = [num2str(ffun.a1)];
text(0.4,0.6,txt,'FontSize',14,'Color',mymagenta)
txt = [num2str(ffun.a2)];
text(0.6,0.6,txt,'FontSize',14,'Color',mymagenta)

xdata = ZenIntForLogistic;
ydata = ZenProbForLogistic;
HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
[ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
%fitmax = round(max(xdata),2,'significant');
yfit=feval(ffun,0:0.001:50); %Fitted function
plot(0:0.001:50,yfit,'Color',mygreen,'LineWidth',2);
txt = [num2str(ffun.a1)];
text(0.4,0.4,txt,'FontSize',14,'Color',mygreen)
txt = [num2str(ffun.a2)];
text(0.6,0.4,txt,'FontSize',14,'Color',mygreen)

txt = ['Int Time: ' num2str(inttime)];
text(0,0,txt,'FontSize',14,'Color','black')
xlim([0,1])
ylabel('Probability')
xlabel('Integrated Mad')
hold off

figure
hold on
fontsize(14,"points")
startPoints = [1 1];
s = fitoptions('Method','NonlinearLeastSquares',... %
    'Lower',[0 0],...
    'Upper',[inf inf],...
    'Startpoint',startPoints);
xdata = [IntForLogistic]';
ydata = [ProbForLogistic]';
HillEqn = fittype( 'x.^a1./(x.^a1+a2.^a1)','options',s);
%Logistic = fittype('1/(1+exp(-(x-a1)/a2))');
[ffun,gofrr] = fit(xdata(:),ydata(:),HillEqn);
yfit=feval(ffun,0:0.0001:5); %Fitted function
scatter(SogIntForLogistic,SogProbForLogistic,'MarkerEdgeColor',mymagenta,'MarkerFaceColor','none');
scatter(IntForLogistic,ProbForLogistic,'MarkerEdgeColor',myblue,'MarkerFaceColor','none')
plot(0:0.0001:5,yfit,'Color','black','LineWidth',2);
txt = ['Error: ' num2str(gofrr.rmse)];
text(0.1,0.05,txt,'FontSize',14)
gofrr
xlim([0,1])
xlabel('Integrated Mad')
ylabel('Probability')
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
yfit=feval(ffun,0:0.00001:0.5); %Fitted function
scatter(SogMadForLogistic,SogProbForLogistic,'MarkerEdgeColor',mymagenta,'MarkerFaceColor','none');
scatter(MadForLogistic,ProbForLogistic,'MarkerEdgeColor',myblue,'MarkerFaceColor','none')
plot(0:0.00001:0.5,yfit,'Color','black','LineWidth',2);
txt = ['Error: ' num2str(gofrr.rmse)];
text(0.05,0.1,txt,'FontSize',14)
%gofrr
xlim([0,0.1])
xlabel('Instantaneous Mad')
ylabel('Probability')
hold off

%% Now predict when you get in this range

% TimeOnInterval=zeros(21,2);
% for x = 1:21
%     uppersogint = UpperCumInt(x,:);
%     lowersogint = LowerCumInt(x,:);
%     timeoninterval = [];
%     for t = 1:60
%         upperintatt = uppersogint(1,t);
%         lowerintatt = lowersogint(1,t);
%         if upperintatt>=minint && lowerintatt<maxint
%             timeoninterval = [timeoninterval,t];
%         else
%             continue
%         end
%     end
%     if ~isempty(timeoninterval);
%         TimeOnInterval(x,:) = [min(timeoninterval),max(timeoninterval)];
%     end
% end
% 
% SogTimeOnInterval=zeros(21,2);
% for x = 1:21
%     uppersogint = SogUpperCumInt(x,:);
%     lowersogint = SogLowerCumInt(x,:);
%     avgsogint = (uppersogint+lowersogint)./2;
%     timeoninterval = [];
%     for t = 1:60
%         upperintatt = uppersogint(1,t);
%         lowerintatt = lowersogint(1,t);
%         meanintatt = avgsogint(1,t);
%         if upperintatt>=minint && lowerintatt<maxint
%             timeoninterval = [timeoninterval,t];
%         else
%             continue
%         end
%     end
%     if ~isempty(timeoninterval);
%         SogTimeOnInterval(x,:) = [min(timeoninterval),max(timeoninterval)];
%     end
% end
% 
% AvgTimeOn = mean(TimeOnInterval,2);
% SogAvgTimeOn = mean(SogTimeOnInterval,2);

% f=figure();
% hold on
% fontsize(gcf,14,'points');
% avg_data = AvgTimeOn';
% x = 0:3:60;
% target = TimeOnInterval(:,1)~=0;
% x = x(target');
% fill([x, flip(x)], [nonzeros(TimeOnInterval(:,1))', flip(nonzeros(TimeOnInterval(:,2))')], [0.85, 0.11, 0.38],'FaceAlpha',0.2, 'EdgeColor','none')
% target2 = avg_data(1,:)~=0;
% x = 0:3:60;
% wtx = x(target2);
% avgdata = nonzeros(avg_data);
% plot(wtx, avgdata,'LineWidth',2,'Color',mymagenta)
% 
% sogavg_data = SogAvgTimeOn';
% x = 0:3:60;
% target = SogTimeOnInterval(:,1)~=0;
% x = x(target');
% fill([x, flip(x)], [nonzeros(SogTimeOnInterval(:,1))', flip(nonzeros(SogTimeOnInterval(:,2))')], [0, 0.3, 0.25],'FaceAlpha',0.2, 'EdgeColor','none')
% target2 = sogavg_data(1,:)~=0;
% x = 0:3:60;
% sogx = x(target2);
% sogavgdata = nonzeros(sogavg_data);
% plot(sogx, sogavgdata,'LineWidth',2,'Color',mygreen)
% 
% x = 0:3:60;
% scatter(x(1,MutTimeOn(1,1:21)~=0),nonzeros(MutTimeOn(1,1:21)),'MarkerEdgeColor',mygreen,'MarkerFaceColor',mygreen)
% scatter(x(1,TimeOn(1,1:21)~=0),nonzeros(TimeOn(1,1:21)),'MarkerEdgeColor',mymagenta,'MarkerFaceColor',mymagenta)
% ylim([20,60])
% xlim([0,60])
% xlabel("Position (um from midline)")
% ylabel("Time On")
% saveas(f,strcat(FigurePath,'sogpredictionmodel_ush'),'svg');
% hold off

%% 
% cmap = cbrewer2('Spectral',101);
% figure
% hold on
% for x = 1:21
%     xdata = 1:60;
%     ydata = MeanCumInt(x,:);
%     c = cmap(101-x,:);
%     plot(xdata,ydata,'Color',c)
% end
% hold off
% 
% cmap = cbrewer2('Spectral',101);
% figure
% hold on
% for x = 1:21
%     xdata = 1:60;
%     ydata = SogMeanCumInt(x,:);
%     c = cmap(101-x,:);
%     plot(xdata,ydata,'Color',c)
% end
% hold off
% 
% cmap = cbrewer2('Spectral',101);
% figure
% hold on
% for x = 1:21
%     xdata = 1:60;
%     ydata = BinnedMed(x,:);
%     c = cmap(101-x,:);
%     plot(xdata,ydata,'Color',c)
% end
% hold off
% 
% cmap = cbrewer2('Spectral',101);
% figure
% hold on
% for x = 1:21
%     xdata = 1:60;
%     ydata = BinnedMedSog(x,:);
%     c = cmap(101-x,:);
%     plot(xdata,ydata,'Color',c)
% end
% hold off
% 
% cmap = cbrewer2('Spectral',101);
% figure
% hold on
% for x = 1:21
%     xdata = 1:60';
%     ydata = FractionOn(:,x);
%     c = cmap(101-x,:);
%     plot(xdata,ydata,'Color',c)
% end
% hold off
% 
% cmap = cbrewer2('Spectral',101);
% figure
% hold on
% for x = 1:21
%     xdata = 1:60';
%     ydata = MutFractionOn(:,x);
%     c = cmap(101-x,:);
%     plot(xdata,ydata,'Color',c)
% end
% hold off


    