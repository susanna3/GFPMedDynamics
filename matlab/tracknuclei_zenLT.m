function [NucleiOverTime,TrackedIDArray,TrackedNucZenArray, TrackedRatioArray, TrackedXArray, TrackedYArray] = tracknuclei_zenLT( ...
    NucleiStacked,TimeDistanceThreshold, X_index,Y_index,ObjectNum_index,Nuc_index,Ratio_index,scale_factor)

%% set test parameters 

% NucleiStacked = frame_array;
% X_index = 3;
% Y_index = 4;
% ObjectNum_index = 2;
% scale_factor = 7;
% Nuc_index = 9;
% Ratio_index = 8;
% TimeDistanceThreshold = 50;
%%

NucleiOverTime = {};

NumTimepoints = length(NucleiStacked);
[NumNucs, NumObvs] = size(NucleiStacked{1});

for ii = 1:NumTimepoints-1
    timepoint_mat = NucleiStacked{ii};
    [num_nuclei num_data] = size(timepoint_mat);
    compare_mat =  NucleiStacked{ii+1}; %compare to next timepoint
    %compare_mat = timepoint_array{ii-1} %compare to previous timepoint
    tracked_mat = [];
    for nucleus = 1:num_nuclei
        nucleus_info = timepoint_mat(nucleus,:);
        centerX = timepoint_mat(nucleus,X_index);
        centerY = timepoint_mat(nucleus,Y_index);
        a = compare_mat(:,X_index)-centerX;
        b = compare_mat(:,Y_index)-centerY;
        c = sqrt(a.^2 + b.^2);
        distance_mat = [compare_mat c];
        tracked_distance = min(c);
        if tracked_distance < TimeDistanceThreshold
            tracked_nucleus = distance_mat(find(distance_mat(:,num_data + 1)==tracked_distance),:);
            tracked_info = [nucleus_info, tracked_nucleus];
            tracked_mat = [tracked_mat;tracked_info];
            continue
        else
            new = [0];
            nucleus_lost = [nucleus_info, zeros(1,NumObvs+1)];
            tracked_mat = [tracked_mat;nucleus_lost];
            continue
    
        end
    end
    NucleiOverTime = [NucleiOverTime; tracked_mat];
end

%% Order nuclei in a time course

%initialize arrays
timepoint_tracked_nucleus_array = NucleiOverTime;

TrackedIDArray = [];
TrackedNucZenArray = [];
TrackedRatioArray = [];
TrackedXArray = [];
TrackedYArray = [];

%start from first frame of movie
for tt = 1:NumTimepoints - 1
    starting_nuclei = timepoint_tracked_nucleus_array{tt};
    Tracked_NucID = []; 
    Tracked_NucZen = [];
    Tracked_RatioZen = [];
    Tracked_X = [];
    Tracked_Y = [];
    for nucleus = 1:length(starting_nuclei)
        ID_tracked = zeros(1,tt-1);
        Nuc_tracked = zeros(1,tt-1);
        Ratio_tracked = zeros(1,tt-1);
        Y_tracked = zeros(1,tt-1);
        X_tracked = zeros(1,tt-1);
        starting_nuc = starting_nuclei(nucleus,ObjectNum_index);
        current_nuc = starting_nuc;
        if starting_nuc(1,1) == 0 %if 0 then it has already been tracked starting from a previous timepoint
            continue
        else
            for xx = tt:NumTimepoints-1 %go through each timepoint for each nucleus
                nuclei_info = timepoint_tracked_nucleus_array{xx}(find( ...
                    timepoint_tracked_nucleus_array{xx}(:,ObjectNum_index)==current_nuc(1:1)),:);
                [l w] = size(nuclei_info);
                if l == 0;
                    nuc_lost = zeros(1,NumTimepoints-xx);
                    ID_tracked = [ID_tracked current_nuc nuc_lost];
                    Nuc_tracked = [Nuc_tracked current_zen nuc_lost];
                    Ratio_tracked = [Ratio_tracked current_ratio nuc_lost];
                    X_tracked = [X_tracked current_X nuc_lost];
                    Y_tracked = [Y_tracked current_Y nuc_lost];
                    break
                else
                    current_nuc = nuclei_info(1,ObjectNum_index);
                    current_zen = nuclei_info(1,Nuc_index);
                    current_ratio = nuclei_info(1,Ratio_index);
                    current_X = nuclei_info(1,X_index);
                    current_Y = nuclei_info(1,Y_index);
                    next_nuc = nuclei_info(1,NumObvs + ObjectNum_index);
                    if next_nuc(1,1) == 0
                        nuc_lost = zeros(1,NumTimepoints-xx);
                        ID_tracked = [ID_tracked current_nuc nuc_lost];
                        Nuc_tracked = [Nuc_tracked current_zen nuc_lost];
                        Ratio_tracked = [Ratio_tracked current_ratio nuc_lost];
                        X_tracked = [X_tracked current_X nuc_lost];
                        Y_tracked = [Y_tracked current_Y nuc_lost];
                        timepoint_tracked_nucleus_array{xx}(find( ...
                            timepoint_tracked_nucleus_array{xx}(:,ObjectNum_index)==current_nuc(1:1)),:) = 0;
                        break
                    else
                        if xx == NumTimepoints-1
                            next_nuc = nuclei_info(1,NumObvs + ObjectNum_index);
                            next_zen = nuclei_info(1,NumObvs + Nuc_index);
                            next_ratio = nuclei_info(1,NumObvs + Ratio_index);
                            next_X = nuclei_info(1,NumObvs + X_index);
                            next_Y = nuclei_info(1,NumObvs + Y_index);

                            ID_tracked = [ID_tracked current_nuc next_nuc];
                            Nuc_tracked = [Nuc_tracked current_zen next_zen];
                            Ratio_tracked = [Ratio_tracked current_ratio next_ratio];
                            X_tracked = [X_tracked current_X next_X];
                            Y_tracked = [Y_tracked current_Y next_Y];
                            break
                        else
                            ID_tracked = [ID_tracked current_nuc];
                            Nuc_tracked = [Nuc_tracked current_zen];
                            Ratio_tracked = [Ratio_tracked current_ratio];
                            X_tracked = [X_tracked current_X];
                            Y_tracked = [Y_tracked current_Y];
                            prev_nuc = current_nuc;
                            current_nuc = next_nuc;
                            timepoint_tracked_nucleus_array{xx}(find( ...
                                timepoint_tracked_nucleus_array{xx}(:,ObjectNum_index)==prev_nuc(1:1)),:) = 0;
                            continue
                        end
                    end
                end
            end
        end
        Tracked_NucID = [Tracked_NucID; ID_tracked];
        Tracked_NucZen = [Tracked_NucZen; Nuc_tracked];
        Tracked_RatioZen = [Tracked_RatioZen; Ratio_tracked];
        Tracked_X = [Tracked_X; X_tracked];
        Tracked_Y = [Tracked_Y; Y_tracked];
    end
    TrackedIDArray = [TrackedIDArray; Tracked_NucID];
    TrackedNucZenArray = [TrackedNucZenArray; Tracked_NucZen];
    TrackedRatioArray = [TrackedRatioArray; Tracked_RatioZen];
    TrackedXArray = [TrackedXArray; Tracked_X];
    TrackedYArray = [TrackedYArray; Tracked_Y];
end

end