%% Compute average intensity and average X/Y for each nucleus
function NucleiStacked = avgperstack( ...
    StackedNucleusArray,StackedAreaArray,StackedMedArray,StackedXArray,StackedYArray,StackedDotCountArray,StackedDotIntArray,...
    NumTimepoints)

NucleiStacked = {};

for matrix = 1:NumTimepoints;
    TimepointInfo = [];
    stacked_nuclei = StackedNucleusArray{matrix};
    stacked_area = StackedAreaArray{matrix};
    stacked_med = StackedMedArray{matrix};
    stacked_x = StackedXArray{matrix};
    stacked_y = StackedYArray{matrix};
    stacked_dotcount = StackedDotCountArray{matrix};
    stacked_dotint = StackedDotIntArray{matrix};
    count = 1;
    
    %id_mean = sum(stacked_nuclei,2) ./ sum(stacked_nuclei~=0,2);
    id_mean = [];
    area_sum = sum(stacked_area,2);
    med_sum = sum(stacked_med,2);
    med_mean = med_sum ./ area_sum;
    x_mean = sum(stacked_x,2) ./ sum(stacked_x~=0,2);
    y_mean = sum(stacked_y,2) ./ sum(stacked_y~=0,2);
    dot_sum = sum(stacked_dotcount,2);
    dot_mean = sum(stacked_dotint,2) ./ sum(stacked_dotint~=0,2);

    [num_nuc num_zs] = size(stacked_nuclei);
    for nuc = 1:num_nuc
        nuc_id = count;
        count = count+1;
        id_mean = [id_mean; nuc_id];
    end
    TimepointInfo = [id_mean,area_sum,med_sum,med_mean,x_mean,y_mean,dot_sum,dot_mean];
    NucleiStacked = [NucleiStacked TimepointInfo];
end

cmap=cbrewer2('Blues',60);
figure
hold on
for mat = 1:NumTimepoints
    %figure
    xData = NucleiStacked{mat}(:,6);
    %xData = xData./PixelToMicron;
    yData = NucleiStacked{mat}(:,4);
    c = cmap(mat+1,:);
    sz = 10;
    %scatter(xData,yData,sz,c);
    %ylim([0 0.2]);
    %xlim([200 800]);
    ft = fit(xData,yData,'smoothingspline','SmoothingParam',0.00005);
    plot(sortrows(xData),feval(ft,sortrows(xData)),"Color",c,"LineWidth",2);
    xlabel("Pixels")
    ylabel("Nuclear Med")
end
hold off

end