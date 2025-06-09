function [NucleiOverTime,TrackedNucleusArray,TrackedMedArray,TrackedXArray,TrackedYArray,TrackedDotCountArray,TrackedDotIntArray] = tracknuclei( ...
    NucleiStacked,TimeDistanceThreshold,NumTimepoints,PixelCutoffLow,PixelCutoffHigh, ...
    X_index,Y_index,ObjectNum_index,MedIntensity_index,DotCount_index,DotInt_index, ...
    NextNuc_index,NextMed_index,NextX_index,NextY_index,NextDotCount_index,NextDotInt_index)

NucleiOverTime = {};

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
            tracked_object = tracked_nucleus(ObjectNum_index);
            tracked_centerX = tracked_nucleus(X_index);
            tracked_centerY = tracked_nucleus(Y_index);
            tracked_Medlevels = tracked_nucleus(MedIntensity_index);
            tracked_DotCount = tracked_nucleus(DotCount_index);
            tracked_DotInt = tracked_nucleus(DotInt_index);
            tracked_info = [nucleus_info tracked_object tracked_Medlevels tracked_centerX tracked_centerY tracked_DotCount tracked_DotInt];
            tracked_mat = [tracked_mat;tracked_info];
            continue
        else
            new = [0];
            nucleus_lost = [nucleus_info new new new new new new];
            tracked_mat = [tracked_mat;nucleus_lost];
            continue
    
        end
    end
    NucleiOverTime = [NucleiOverTime; tracked_mat];
end

%% Order nuclei in a time course

%initialize arrays
timepoint_tracked_nucleus_array = {};
for nuc_mat = 1:NumTimepoints-1;
    dataset = NucleiOverTime{nuc_mat};
    filterbyy = dataset(dataset(:,Y_index)>PixelCutoffLow & dataset(:,Y_index)<PixelCutoffHigh,:);
    filterbyy = filterbyy(filterbyy(:,NextY_index)>PixelCutoffLow & filterbyy(:,NextY_index)<PixelCutoffHigh,:);
    timepoint_tracked_nucleus_array = [timepoint_tracked_nucleus_array; filterbyy];
end

TrackedNucleusArray = [];
TrackedMedArray = [];
TrackedXArray = [];
TrackedYArray = [];
TrackedDotCountArray = [];
TrackedDotIntArray = [];

%start from first frame of movie
for tt = 1:NumTimepoints - 1;
    starting_nuclei = timepoint_tracked_nucleus_array{tt};
    Tracked_Nuclei = [];
    Tracked_Med = [];
    Tracked_X = [];
    Tracked_Y = [];
    Tracked_DotCount = [];
    Tracked_DotInt = [];
    [leng_start wid_start] = size(starting_nuclei);
    for nucleus = 1:leng_start
        Nuc_tracked = zeros(1,tt-1);
        Med_tracked = zeros(1,tt-1);
        X_tracked = zeros(1,tt-1);
        Y_tracked = zeros(1,tt-1);
        DotCount_tracked = zeros(1,tt-1);
        DotInt_tracked = zeros(1,tt-1);

        starting_nuc = starting_nuclei(nucleus,ObjectNum_index);
        current_nuc = starting_nuc;
        if starting_nuc(1,1) == 0; %if 0 then it has already been tracked starting from a previous timepoint
            continue
        else
            for xx = tt:NumTimepoints-1; %go through each timepoint for each nucleus
                nuclei_info = timepoint_tracked_nucleus_array{xx}(find( ...
                    timepoint_tracked_nucleus_array{xx}(:,ObjectNum_index)==current_nuc(1:1)),:);
                [l w] = size(nuclei_info);
                if l == 0;
                    nuc_lost = zeros(1,NumTimepoints-xx);
                    Nuc_tracked = [Nuc_tracked current_nuc nuc_lost];
                    Med_tracked = [Med_tracked current_Med nuc_lost];
                    X_tracked = [X_tracked current_X nuc_lost];
                    Y_tracked = [Y_tracked current_Y nuc_lost];
                    DotCount_tracked = [DotCount_tracked current_DotCount nuc_lost];
                    DotInt_tracked = [DotInt_tracked current_DotInt nuc_lost];
                    break
                else
                    current_nuc = nuclei_info(1,ObjectNum_index);
                    current_Med = nuclei_info(1,MedIntensity_index);
                    current_X = nuclei_info(1,X_index);
                    current_Y = nuclei_info(1,Y_index);
                    current_DotCount = nuclei_info(1,DotCount_index);
                    current_DotInt = nuclei_info(1,DotInt_index);
                    next_nuc = nuclei_info(1,NextNuc_index);
                    if next_nuc(1,1) == 0;
                        nuc_lost = zeros(1,NumTimepoints-xx);
                        Nuc_tracked = [Nuc_tracked current_nuc nuc_lost];
                        Med_tracked = [Med_tracked current_Med nuc_lost];
                        X_tracked = [X_tracked current_X nuc_lost];
                        Y_tracked = [Y_tracked current_Y nuc_lost];
                        DotCount_tracked = [DotCount_tracked current_DotCount nuc_lost];
                        DotInt_tracked = [DotInt_tracked current_DotInt nuc_lost];
                        timepoint_tracked_nucleus_array{xx}(find( ...
                            timepoint_tracked_nucleus_array{xx}(:,ObjectNum_index)==current_nuc(1:1)),:) = 0;
                        break
                    else
                        if xx == NumTimepoints-1;
                            next_Med = nuclei_info(1,NextMed_index);
                            next_X = nuclei_info(1,NextX_index);
                            next_Y = nuclei_info(1,NextY_index);
                            next_DotCount = nuclei_info(1,NextDotCount_index);
                            next_DotInt = nuclei_info(1,NextDotInt_index);
                            Nuc_tracked = [Nuc_tracked current_nuc next_nuc];
                            Med_tracked = [Med_tracked current_Med next_Med];
                            X_tracked = [X_tracked current_X next_X];
                            Y_tracked = [Y_tracked current_Y next_Y];
                            DotCount_tracked = [DotCount_tracked current_DotCount next_DotCount];
                            DotInt_tracked = [DotInt_tracked current_DotInt next_DotInt];
                            break
                        else
                            Nuc_tracked = [Nuc_tracked current_nuc];
                            Med_tracked = [Med_tracked current_Med];
                            X_tracked = [X_tracked current_X];
                            Y_tracked = [Y_tracked current_Y];
                            DotCount_tracked = [DotCount_tracked current_DotCount];
                            DotInt_tracked = [DotInt_tracked current_DotInt];
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
        Tracked_Nuclei = [Tracked_Nuclei; Nuc_tracked];
        Tracked_Med = [Tracked_Med; Med_tracked];
        Tracked_X = [Tracked_X; X_tracked];
        Tracked_Y = [Tracked_Y; Y_tracked];
        Tracked_DotCount = [Tracked_DotCount; DotCount_tracked];
        Tracked_DotInt = [Tracked_DotInt; DotInt_tracked];
    end
    TrackedNucleusArray = [TrackedNucleusArray; Tracked_Nuclei];
    TrackedMedArray = [TrackedMedArray; Tracked_Med];
    TrackedXArray = [TrackedXArray; Tracked_X];
    TrackedYArray = [TrackedYArray; Tracked_Y];
    TrackedDotCountArray = [TrackedDotCountArray; Tracked_DotCount];
    TrackedDotIntArray = [TrackedDotIntArray; Tracked_DotInt];
end

end