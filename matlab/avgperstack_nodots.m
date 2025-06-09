%% Compute average intensity and average X/Y for each nucleus
function NucleiStacked = avgperstack_nodots( ...
    StackedNucleusArray,StackedAreaArray,StackedMedArray,StackedXArray,StackedYArray,...
    NumTimepoints)

NucleiStacked = {};

for matrix = 1:NumTimepoints;
    TimepointInfo = [];
    stacked_nuclei = StackedNucleusArray{matrix};
    stacked_area = StackedAreaArray{matrix};
    stacked_med = StackedMedArray{matrix};
%     stacked_mad = StackedMadArray{matrix};
    stacked_x = StackedXArray{matrix};
    stacked_y = StackedYArray{matrix};
    count = 1;
    
    %id_mean = sum(stacked_nuclei,2) ./ sum(stacked_nuclei~=0,2);
    id_mean = [];
    area_sum = sum(stacked_area,2);
    med_sum = sum(stacked_med,2);
    med_mean = med_sum ./ area_sum;
%     mad_sum = sum(stacked_mad,2);
%     mad_mean = mad_sum ./ area_sum;
    x_mean = sum(stacked_x,2) ./ sum(stacked_x~=0,2);
    y_mean = sum(stacked_y,2) ./ sum(stacked_y~=0,2);

    [num_nuc num_zs] = size(stacked_nuclei);
    for nuc = 1:num_nuc
        nuc_id = count;
        count = count+1;
        id_mean = [id_mean; nuc_id];
    end
    TimepointInfo = [id_mean,area_sum,med_sum,med_mean,x_mean,y_mean];
    NucleiStacked = [NucleiStacked TimepointInfo];
end

end