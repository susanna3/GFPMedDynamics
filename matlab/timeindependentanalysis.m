%% Plot Hill Curves with no time bins

genename = 'Tup';

% Load data
AllWTData = load(strcat(genename,'WTNucleiData.mat'));
AllWTTiming = load(strcat(genename,'WTNucleiTiming.mat'));
AllSogData = load(strcat(genename,'SogNucleiData.mat'));
AllSogTiming = load(strcat(genename,'SogNucleiTiming.mat'));
AllZenData = load(strcat(genename,'ZenNucleiData.mat'));
AllZenTiming = load(strcat(genename,'ZenNucleiTiming.mat'));

ZenWTData = load('ZenWTNucleiData.mat');
ZenWTTiming = load('ZenWTNucleiTiming.mat');
ZenSogData = load('ZenSogNucleiData.mat');
ZenSogTiming = load('ZenSogNucleiTiming.mat');

AllWTData = AllWTData.AllWTNucleiData;
AllWTTiming = AllWTTiming.AllWTNucleiTiming;
AllSogData = AllSogData.AllMutantNucleiData;
AllSogTiming = AllSogTiming.AllMutantNucleiTiming;
AllZenData = AllZenData.AllMutantNucleiData;
AllZenTiming = AllZenTiming.AllMutantNucleiTiming;

ZenWTData = ZenWTData.AllWTNucleiData;
ZenWTTiming = ZenWTTiming.AllWTNucleiTiming;
ZenSogData = ZenSogData.AllMutantNucleiData;
ZenSogTiming = ZenSogTiming.AllMutantNucleiTiming;

binsize = 3;
numbins = 21;
[l1,w]=size(AllWTTiming);
[l2,w]=size(AllSogTiming);
[l3,w]=size(AllZenTiming);

[l4,w]=size(ZenWTTiming);
[l5,w]=size(ZenSogTiming);

%% Check timing
% AllTiming = [AllWTTiming;AllSogTiming];
% AllTiming = AllTiming(AllTiming(:,1)>29.5,:);
% mean(AllTiming)

%% Measure distance binned values

WTBinnedMad = [];
WTBinnedRNA = [];
WTBinnedProb = [];

for embryo = 1:l1
    embryodata = sortrows(AllWTData{embryo},1);
    [nucno,colno] = size(embryodata);
    distbins = discretize(embryodata(:,1),[0:binsize:60,200],'IncludedEdge','right');
    binnedmaddata = zeros(1,numbins);
    binnedrnadata = zeros(1,numbins);
    binnedprob = zeros(1,numbins);
    for bin = 1:numbins
        distbindata = embryodata(distbins(:,1)==bin,:);
        distbindata(isnan(distbindata))=0;
        avgmadinbin = mean(distbindata(:,2));
        avgrnainbin = mean(distbindata(:,4));
        binnedmaddata(1,bin)=avgmadinbin;
        binnedrnadata(1,bin)=avgrnainbin;
        totalnuc = length(distbindata(:,1));
        if totalnuc > 9
            positivenuc = length(find(distbindata(:,3)>0));
            dotprobinbin = positivenuc./totalnuc;
            binnedprob(1,bin)=dotprobinbin;
        elseif totalnuc < 10
            dotprobinbin = 0;
            binnedprob(1,bin)=dotprobinbin;
        end
    end
    WTBinnedMad=vertcat(WTBinnedMad,binnedmaddata);
    WTBinnedRNA=vertcat(WTBinnedRNA,binnedrnadata);
    WTBinnedProb=vertcat(WTBinnedProb,binnedprob);
end

SogBinnedMad = [];
SogBinnedRNA = [];
SogBinnedProb = [];

for embryo = 1:l2
    embryodata = sortrows(AllSogData{embryo},1);
    [nucno,colno] = size(embryodata);
    distbins = discretize(embryodata(:,1),[0:binsize:60,200],'IncludedEdge','right');
    binnedmaddata = zeros(1,numbins);
    binnedrnadata = zeros(1,numbins);
    binnedprob = zeros(1,numbins);
    for bin = 1:numbins
        distbindata = embryodata(distbins(:,1)==bin,:);
        distbindata(isnan(distbindata))=0;
        avgmadinbin = mean(distbindata(:,2));
        avgrnainbin = mean(distbindata(:,4));
        binnedmaddata(1,bin)=avgmadinbin;
        binnedrnadata(1,bin)=avgrnainbin;
        totalnuc = length(distbindata(:,1));
        if totalnuc > 9
            positivenuc = length(find(distbindata(:,3)>0));
            dotprobinbin = positivenuc./totalnuc;
            binnedprob(1,bin)=dotprobinbin;
        elseif totalnuc < 10
            dotprobinbin = 0;
            binnedprob(1,bin)=dotprobinbin;
        end
    end
    SogBinnedMad=vertcat(SogBinnedMad,binnedmaddata);
    SogBinnedRNA=vertcat(SogBinnedRNA,binnedrnadata);
    SogBinnedProb=vertcat(SogBinnedProb,binnedprob);
end

ZenBinnedMad = [];
ZenBinnedRNA = [];
ZenBinnedProb = [];

for embryo = 1:l3
    embryodata = sortrows(AllZenData{embryo},1);
    [nucno,colno] = size(embryodata);
    distbins = discretize(embryodata(:,1),[0:binsize:60,200],'IncludedEdge','right');
    binnedmaddata = zeros(1,numbins);
    binnedrnadata = zeros(1,numbins);
    binnedprob = zeros(1,numbins);
    for bin = 1:numbins
        distbindata = embryodata(distbins(:,1)==bin,:);
        distbindata(isnan(distbindata))=0;
        avgmadinbin = mean(distbindata(:,2));
        avgrnainbin = mean(distbindata(:,4));
        binnedmaddata(1,bin)=avgmadinbin;
        binnedrnadata(1,bin)=avgrnainbin;
        totalnuc = length(distbindata(:,1));
        if totalnuc > 9
            positivenuc = length(find(distbindata(:,3)>0));
            dotprobinbin = positivenuc./totalnuc;
            binnedprob(1,bin)=dotprobinbin;
        elseif totalnuc < 10
            dotprobinbin = 0;
            binnedprob(1,bin)=dotprobinbin;
        end
    end
    ZenBinnedMad=vertcat(ZenBinnedMad,binnedmaddata);
    ZenBinnedRNA=vertcat(ZenBinnedRNA,binnedrnadata);
    ZenBinnedProb=vertcat(ZenBinnedProb,binnedprob);
end

%% Save MATLAB objects of distance binned data

% save(strcat(genename,'WTDistBinnedMad'),'WTBinnedMad');
% save(strcat(genename,'WTDistBinnedRNA'),'WTBinnedRNA');
% save(strcat(genename,'WTDistBinnedProb'),'WTBinnedProb');
% 
% save(strcat(genename,'SogDistBinnedMad.mat'),'SogBinnedMad');
% save(strcat(genename,'SogDistBinnedRNA'),'SogBinnedRNA');
% save(strcat(genename,'SogDistBinnedProb'),'SogBinnedProb');
% 
% save(strcat(genename,'ZenDistBinnedMad.mat'),'ZenBinnedMad');
% save(strcat(genename,'ZenDistBinnedRNA'),'ZenBinnedRNA');
% save(strcat(genename,'ZenDistBinnedProb'),'ZenBinnedProb');
% 
% save(strcat(genename,'WTTiming'),'AllWTTiming');
% save(strcat(genename,'SogTiming'),'AllSogTiming');
% save(strcat(genename,'ZenTiming'),'AllZenTiming');


%% Plot relationship between Mad and RNA

WTBinnedMad = reshape(WTBinnedMad,[],1);
WTBinnedProb = reshape(WTBinnedProb,[],1);
WTBinnedRNA = reshape(WTBinnedRNA,[],1);
SogBinnedMad = reshape(SogBinnedMad,[],1);
SogBinnedProb = reshape(SogBinnedProb,[],1);
SogBinnedRNA = reshape(SogBinnedRNA,[],1);
ZenBinnedMad = reshape(ZenBinnedMad,[],1);
ZenBinnedProb = reshape(ZenBinnedProb,[],1);
ZenBinnedRNA = reshape(ZenBinnedRNA,[],1);

%%
figure
hold on
title(genename)
xlabel('Mad Concentration')
ylabel('RNA ')
scatter(WTBinnedMad,WTBinnedRNA,5,'MarkerFaceColor','black','MarkerEdgeColor','black')
scatter(SogBinnedMad,SogBinnedRNA,5,'MarkerFaceColor','magenta','MarkerEdgeColor','magenta')
scatter(ZenBinnedMad,ZenBinnedRNA,5,'MarkerFaceColor','green','MarkerEdgeColor','green')
hold off

figure
hold on
title(genename)
xlabel('log(Mad)')
ylabel('log(P/1-P)')
%ylim([-4,6]);
%xlim([-4,1]);
scatter(log(WTBinnedMad),log(WTBinnedRNA./(1-WTBinnedRNA)),5,'MarkerFaceColor','black','MarkerEdgeColor','black')
scatter(log(SogBinnedMad),log(SogBinnedRNA./(1-SogBinnedRNA)),5,'MarkerFaceColor','magenta','MarkerEdgeColor','magenta')
scatter(log(ZenBinnedMad),log(ZenBinnedRNA./(1-ZenBinnedRNA)),5,'MarkerFaceColor','green','MarkerEdgeColor','green')
hold off

%%

figure
hold on
title(genename)
xlabel('Avg Nuclear pSMAD')
ylabel('Fraction Nuclei On')
fontsize(14,'point')
scatter(WTBinnedMad,WTBinnedProb,30,'MarkerFaceColor','#1E88E5','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(SogBinnedMad,SogBinnedProb,30,'MarkerFaceColor','#D81B60','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
scatter(ZenBinnedMad,ZenBinnedProb,30,'MarkerFaceColor','#004D40','MarkerEdgeColor','none','MarkerFaceAlpha',0.5)
s = fitoptions('Method','NonlinearLeastSquares','Robust','LAR', ...
   'Lower',[0 0],...
   'Upper',[inf inf],...
   'Startpoint',[1 1]);
HillEqn = fittype( 'x.^a1./((a2^a1) + (x.^a1))','options',s);

[ffun,gofrr] = fit(WTBinnedMad,WTBinnedProb,HillEqn);
yfit=feval(ffun,0:0.001:1); %Fitted function
plot(0:0.001:1,yfit,'Color','#1E88E5','LineWidth',3);
ffun

[ffun,gofrr] = fit(SogBinnedMad,SogBinnedProb,HillEqn);
yfit=feval(ffun,0:0.001:1); %Fitted function
plot(0:0.001:1,yfit,'Color','#D81B60','LineWidth',3);
ffun

[ffun,gofrr] = fit(ZenBinnedMad,ZenBinnedProb,HillEqn);
yfit=feval(ffun,0:0.001:1); %Fitted function
plot(0:0.001:1,yfit,'Color','#004D40','LineWidth',3);
ffun

ylim([0,1]);
xlim([0,1]);
hold off

%%

figure
hold on
title(genename)
xlabel('log(Mad)')
ylabel('log(P/1-P)')
ylim([-2,6]);
xlim([-4,1]);
WTX = log(WTBinnedMad);
WTY = log(WTBinnedProb./(1-WTBinnedProb));
WTX = WTX(isfinite(WTY));
WTY = WTY(isfinite(WTY));
WTX = WTX(WTY(:,1)>-2,:);
WTY = WTY(WTY(:,1)>-2,:);
SogX = log(SogBinnedMad);
SogY = log(SogBinnedProb./(1-SogBinnedProb));
SogX = SogX(isfinite(SogY));
SogY = SogY(isfinite(SogY));
SogX = SogX(SogY(:,1)>-2,:);
SogY = SogY(SogY(:,1)>-2,:);
ZenX = log(ZenBinnedMad);
ZenY = log(ZenBinnedProb./(1-ZenBinnedProb));
ZenX = ZenX(isfinite(ZenY));
ZenY = ZenY(isfinite(ZenY));
ZenX = ZenX(ZenY(:,1)>-2,:);
ZenY = ZenY(ZenY(:,1)>-2,:);
scatter(WTX,WTY,5,'MarkerFaceColor','#D81B60','MarkerEdgeColor','#D81B60')
scatter(SogX,SogY,5,'MarkerFaceColor','#004D40','MarkerEdgeColor','#004D40')
scatter(ZenX,ZenY,5,'MarkerFaceColor','#FFC107','MarkerEdgeColor','#FFC107')
PWT = polyfit(WTX,WTY,1);
PS = polyfit(SogX,SogY,1);
PZ = polyfit(ZenX,ZenY,1);
yfitWT = PWT(1)*WTX+PWT(2);
yfitS = PS(1)*SogX+PS(2);
yfitZ = PZ(1)*ZenX+PZ(2);
plot(WTX,yfitWT,'r-.','Color','#D81B60','LineWidth',2);
plot(SogX,yfitS,'r-.','Color','#004D40','LineWidth',2);
plot(ZenX,yfitZ,'r-.','Color','#FFC107','LineWidth',2);
hold off

%% Plot linear fit and find n and Kd

% Create and Plot Raw Data for Zen
xdata = log(ZenBinnedMad);
ydata = log(ZenBinnedProb./(1-ZenBinnedProb));
alldata = [xdata(:),ydata(:)];
alldata = sortrows(alldata,1);
x = alldata(:,1);
y = alldata(:,2);
x = x(y > -2);
y = y(y > -2);
x = x(y < 2);
y = y(y < 2);
y = y(x > -1.5);
x = x(x > -1.5);
y = y(x < 0);
x = x(x < 0);
%scatter(x,y)
% Fit line to data using polyfit
mdl = fitlm(x,y);
mdl.Coefficients(:,"Estimate")
coefficients = table2array(mdl.Coefficients(:,"Estimate"));
logkd = coefficients(1,1);
kd = exp(-logkd)
hillco = coefficients(2,1)
figure
hold on
plot(mdl)


%% Create and Plot Raw Data for WT/Sog
xdata = [log(WTBinnedMad);log(SogBinnedMad)];
ydata = [log(WTBinnedProb./(1-WTBinnedProb));log(SogBinnedProb./(1-SogBinnedProb))];
alldata = [xdata(:),ydata(:)];
alldata = sortrows(alldata,1);
x = alldata(:,1);
y = alldata(:,2);
x = x(y > -2);
y = y(y > -2);
x = x(y < 2);
y = y(y < 2);
y = y(x < -2);
x = x(x < -2);
y = y(x > -4);
x = x(x > -4);
%scatter(x,y)
% Fit line to data using polyfit
mdl = fitlm(x,y);
mdl.Coefficients(:,"Estimate")
coefficients = table2array(mdl.Coefficients(:,"Estimate"));
logkd = coefficients(1,1);
kd = exp(-logkd)
hillco = coefficients(2,1)

plot(mdl)
hold off

%% Plot relationship between Integrated Mad and RNA



WTBinnedMad = reshape(WTBinnedMad,[],1);
WTBinnedProb = reshape(WTBinnedProb,[],1);
WTBinnedRNA = reshape(WTBinnedRNA,[],1);
SogBinnedMad = reshape(SogBinnedMad,[],1);
SogBinnedProb = reshape(SogBinnedProb,[],1);
SogBinnedRNA = reshape(SogBinnedRNA,[],1);
ZenBinnedMad = reshape(ZenBinnedMad,[],1);
ZenBinnedProb = reshape(ZenBinnedProb,[],1);
ZenBinnedRNA = reshape(ZenBinnedRNA,[],1);

% figure
% hold on
% title(genename)
% xlabel('Mad Concentration')
% ylabel('RNA ')
% scatter(WTBinnedMad,WTBinnedRNA,5,'MarkerFaceColor','black','MarkerEdgeColor','black')
% scatter(SogBinnedMad,SogBinnedRNA,5,'MarkerFaceColor','magenta','MarkerEdgeColor','magenta')
% %scatter(ZenBinnedMad,ZenBinnedRNA,5,'MarkerFaceColor','green','MarkerEdgeColor','green')
% hold off
% 
% figure
% hold on
% title(genename)
% xlabel('log(Mad)')
% ylabel('log(P/1-P)')
% %ylim([-4,6]);
% %xlim([-4,1]);
% scatter(log(WTBinnedMad),log(WTBinnedRNA./(1-WTBinnedRNA)),5,'MarkerFaceColor','black','MarkerEdgeColor','black')
% scatter(log(SogBinnedMad),log(SogBinnedRNA./(1-SogBinnedRNA)),5,'MarkerFaceColor','magenta','MarkerEdgeColor','magenta')
% %scatter(log(ZenBinnedMad),log(ZenBinnedRNA./(1-ZenBinnedRNA)),5,'MarkerFaceColor','green','MarkerEdgeColor','green')
% hold off

figure
hold on
title(genename)
xlabel('Mad Concentration')
ylabel('RNA ')
scatter(WTBinnedMad,WTBinnedProb,5,'MarkerFaceColor','#D81B60','MarkerEdgeColor','#D81B60')
scatter(SogBinnedMad,SogBinnedProb,5,'MarkerFaceColor','#004D40','MarkerEdgeColor','#004D40')
scatter(ZenBinnedMad,ZenBinnedProb,5,'MarkerFaceColor','#FFC107','MarkerEdgeColor','#FFC107')
ylim([0,1]);
xlim([0,1]);
hold off

figure
hold on
title(genename)
xlabel('log(Mad)')
ylabel('log(P/1-P)')
ylim([-4,6]);
xlim([-4,1]);
WTX = log(WTBinnedMad);
WTY = log(WTBinnedProb./(1-WTBinnedProb));
WTX = WTX(isfinite(WTY));
WTY = WTY(isfinite(WTY));
WTX = WTX(WTY(:,1)>-4,:);
WTY = WTY(WTY(:,1)>-4,:);
SogX = log(SogBinnedMad);
SogY = log(SogBinnedProb./(1-SogBinnedProb));
SogX = SogX(isfinite(SogY));
SogY = SogY(isfinite(SogY));
SogX = SogX(SogY(:,1)>-4,:);
SogY = SogY(SogY(:,1)>-4,:);
ZenX = log(ZenBinnedMad);
ZenY = log(ZenBinnedProb./(1-ZenBinnedProb));
ZenX = ZenX(isfinite(ZenY));
ZenY = ZenY(isfinite(ZenY));
ZenX = ZenX(ZenY(:,1)>-4,:);
ZenY = ZenY(ZenY(:,1)>-4,:);
scatter(WTX,WTY,5,'MarkerFaceColor','#D81B60','MarkerEdgeColor','#D81B60')
scatter(SogX,SogY,5,'MarkerFaceColor','#004D40','MarkerEdgeColor','#004D40')
scatter(ZenX,ZenY,5,'MarkerFaceColor','#FFC107','MarkerEdgeColor','#FFC107')
PWT = polyfit(WTX,WTY,1);
PS = polyfit(SogX,SogY,1);
PZ = polyfit(ZenX,ZenY,1);
yfitWT = PWT(1)*WTX+PWT(2);
yfitS = PS(1)*SogX+PS(2);
yfitZ = PZ(1)*ZenX+PZ(2);
plot(WTX,yfitWT,'r-.','Color','#D81B60','LineWidth',2);
plot(SogX,yfitS,'r-.','Color','#004D40','LineWidth',2);
plot(ZenX,yfitZ,'r-.','Color','#FFC107','LineWidth',2);
hold off
