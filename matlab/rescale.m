% Rescale embryos based on one reference embryo

function [scale_factor,FullyTrackedNucleiRescaled] = rescale(D1_binned,NumTimepoints,CorrectedNuclei)

D2 = CorrectedNuclei;
edges = 0:5:50;
D2distance = D2(:,1);
D2_binned_dist = discretize(D2distance,edges);
D2_data = D2(:,2:end);

D2_bybin = {};

for bin = 1:9
    in_bin = [];
    for nucleus = 1:length(D2)
        check_bin = D2_binned_dist(nucleus,1);
        if check_bin == bin
            in_bin = [in_bin; D2_data(nucleus,:)];
        else
            continue
        end
    end
    D2_bybin = [D2_bybin; in_bin];
end  

D2_binned = zeros(9,30);

for bin = 1:9
    for time = 1:30
        lookuptime = NumTimepoints-30+time;
        D2_binned(bin,time) = mean(D2_bybin{bin}(:,lookuptime));
    end
end

scaled_data = [];
for i=1:101
    scale = 0.5+0.01*(i-1);
    compare = sum(sum((D1_binned(1:9,1:30)-(scale*D2_binned)).^2));
    scaled_data = [scaled_data; scale compare];
end
scaled_data = sortrows(scaled_data,2);
scale_factor = scaled_data(1,1);

RescaleMed = CorrectedNuclei(:,3:end).*scale_factor;
FullyTrackedNucleiRescaled = [CorrectedNuclei(:,1:2) RescaleMed];

end
