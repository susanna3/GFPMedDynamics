%% Compute average intensity and average X/Y for each nucleus
function NucleiStacked = avgperstack( ...
    StackedNucleusArray,StackedAreaArray,StackedMedArray,StackedXArray,StackedYArray,StackedDotCountArray,StackedDotIntArray,...
    NumTimepoints)

NucleiStacked = {};

for matrix = 1:NumTimepoints;
    TimepointInfo = [];
    stacked_nuclei = StackedNucleusArray{matrix};
    count = 1;
    [num_nuc num_zs] = size(stacked_nuclei);
    for nuc = 1:num_nuc
        numtracked = stacked_nuclei(nuc,stacked_nuclei(nuc,:)>0);
        [length width] = size(numtracked);
        if width > 4
            stacked_area = StackedAreaArray{matrix}(nuc,:);
            stacked_med = StackedMedArray{matrix}(nuc,:);
            stacked_x = StackedXArray{matrix}(nuc,:);
            stacked_y = StackedYArray{matrix}(nuc,:);
            stacked_dotcount = StackedDotCountArray{matrix}(nuc,:);
            stacked_dotint = StackedDotIntArray{matrix}(nuc,:);
            count = count+1;

            area_sum = sum(stacked_area,2);
            med_sum = sum(stacked_med,2);
            med_mean = med_sum ./ area_sum;
            x_mean = sum(stacked_x,2) ./ sum(stacked_x~=0,2);
            y_mean = sum(stacked_y,2) ./ sum(stacked_y~=0,2);
            dot_sum = sum(stacked_dotcount,2);
            dot_mean = sum(stacked_dotint,2) ./ sum(stacked_dotint~=0,2);

            NucInfo = [count,area_sum,med_sum,med_mean,x_mean,y_mean,dot_sum,dot_mean];
            TimepointInfo = [TimepointInfo;NucInfo];
        else
            continue
        end
    end
    NucleiStacked = [NucleiStacked TimepointInfo];
end

cmap=cbrewer2('Blues',NumTimepoints+1);
figure
hold on
for mat = 1:NumTimepoints
    figure
    xData = NucleiStacked{mat}(:,6);
    %xData = xData./PixelToMicron;
    yData = NucleiStacked{mat}(:,4);
    c = cmap(mat+1,:);
    sz = 10;
    scatter(xData,yData,sz,c);
    ylim([0 0.5]);
    %xlim([200 800]);
end
hold off

end