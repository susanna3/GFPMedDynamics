% Rescale embryos based on one reference embryo

function [CorrectedNucleiRescaled, scale_factor] = rescale_fixed(CorrectedRefNuclei,CorrectedNuclei,RefMinutes,Minutes)

timeedges = 0:5:60;
distedges = 0:5:200;

RefTimeBins = discretize(RefMinutes,timeedges);
TimeBins = discretize(Minutes,timeedges);

RefD = zeros(12,15);
EmbryoD = zeros(12,15);

%% Make reference matrix
RefTimeBinArray = {};

for timebin = 1:12
    RefTimeBinData = [];
    for t = 1:length(RefMinutes);
        bin1 = RefTimeBins(t);
        if timebin == bin1
            dataintimebin = CorrectedRefNuclei{t};
            RefTimeBinData = [RefTimeBinData;dataintimebin];
        else
            continue
        end
    end
    if size(RefTimeBinData) == 0;
        RefTimeBinData = 0;
        RefTimeBinArray = [RefTimeBinArray;RefTimeBinData];
    else
        RefTimeBinArray = [RefTimeBinArray;RefTimeBinData];
    end
end

for timebin = 1:12
    timebindata = RefTimeBinArray{timebin};
    if timebindata(1,1) == 0
        continue
    else
        RefDistBins = discretize(timebindata(:,1),distedges);
        for distbin = 1:15
            RefBinnedMad = [];
            for n = 1:size(RefDistBins)
                bin2 = RefDistBins(n);
                if bin2 == distbin;
                    madinbin = timebindata(n,2);
                    RefBinnedMad = [RefBinnedMad;madinbin];
                else
                    continue
                end
            end
            RefTimeDistMad = mean(RefBinnedMad);
            RefD(timebin,distbin) = RefTimeDistMad;
        end
    end
end 

%% Make data matrix

TimeBinArray = {};

for timebin = 1:12
    TimeBinData = [];
    for t = 1:length(Minutes);
        bin1 = TimeBins(t);
        if timebin == bin1
            dataintimebin = CorrectedNuclei{t};
            TimeBinData = [TimeBinData;dataintimebin];
        else
            continue
        end
    end
    if size(TimeBinData) == 0;
        TimeBinData = 0;
        TimeBinArray = [TimeBinArray;TimeBinData];
    else
        TimeBinArray = [TimeBinArray;TimeBinData];
    end
end

for timebin = 1:12
    timebindata = TimeBinArray{timebin};
    if timebindata(1,1) == 0
        continue
    else
        DistBins = discretize(timebindata(:,1),distedges);
        for distbin = 1:15
            BinnedMad = [];
            for n = 1:size(DistBins)
                bin2 = DistBins(n);
                if bin2 == distbin;
                    madinbin = timebindata(n,2);
                    BinnedMad = [BinnedMad;madinbin];
                else
                    continue
                end
            end
            TimeDistMad = mean(BinnedMad);
            EmbryoD(timebin,distbin) = TimeDistMad;
        end
    end
end 

%% Rescale

D1_binned = RefD((RefD(:,1)>0 & EmbryoD(:,1)>0),:);
D2_binned = EmbryoD((RefD(:,1)>0 & EmbryoD(:,1)>0),:);

[distancebins timebins] = size(D1_binned);

scaled_data = [];
for i=1:5001
    scale = 0.005+0.005*(i-1);
    compare = sum(sum((D1_binned(1:distancebins,1:timebins)-(scale*D2_binned)).^2));
    scaled_data = [scaled_data; scale compare];
end
scaled_data = sortrows(scaled_data,2);
scale_factor = scaled_data(1,1);

CorrectedNucleiRescaled = {};

for i = 1:length(CorrectedNuclei)
    datatorescale = CorrectedNuclei{i};
    datatorescale(:,2) = datatorescale(:,2)*scale_factor;
    CorrectedNucleiRescaled = [CorrectedNucleiRescaled datatorescale];
end

end
