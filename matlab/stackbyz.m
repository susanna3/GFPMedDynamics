%% Takes a dataset from cell profiler and stacks nuclei in Z

function [StackedNucleusArray,StackedAreaArray,StackedMedArray,StackedXArray,StackedYArray,StackedDotCountArray,StackedDotIntArray] = stackbyz( ...
    frame_array,NumTimepoints,NumZs,DistanceThreshold, ...
    X_index,Y_index,ObjectNum_index,MedIntensity_index,DotCount_index,DotInt_index,Area_index, ...
    NextNuc_index,NextArea_index,NextMed_index,NextX_index,NextY_index,NextDotCount_index,NextDotInt_index)


%% find closest nucleus in z

StackedNucleiArray = {};

for t = 1:NumTimepoints;
    timepoint_frame_array = frame_array(((t*NumZs)-(NumZs-1)):(t*NumZs));
    for z = 1:NumZs-1;
        current_z = timepoint_frame_array{z};
        next_z = timepoint_frame_array{z+1};
        [num_nuclei num_data] = size(current_z);
        stacked_mat = [];
        for nucleus = 1:num_nuclei
            nucleus_info = current_z(nucleus,:);
            centerX = current_z(nucleus,X_index);
            centerY = current_z(nucleus,Y_index);
            a = next_z(:,X_index)-centerX;
            b = next_z(:,Y_index)-centerY;
            c = sqrt(a.^2 + b.^2);
            distance_mat = [next_z c];
            stacked_distance = min(c);
            if stacked_distance < DistanceThreshold
                stacked_nucleus = distance_mat(find(distance_mat(:,num_data+1)==stacked_distance),:);
                 stacked_object = stacked_nucleus(ObjectNum_index);
                 stacked_centerX = stacked_nucleus(X_index);
                 stacked_centerY = stacked_nucleus(Y_index);
                 stacked_Medlevels = stacked_nucleus(MedIntensity_index);
                 stacked_DotCount = stacked_nucleus(DotCount_index);
                 stacked_DotInt = stacked_nucleus(DotInt_index);
                 stacked_Area = stacked_nucleus(Area_index);
                 stacked_info = [nucleus_info stacked_object stacked_Area stacked_Medlevels  stacked_centerX  stacked_centerY stacked_DotCount stacked_DotInt];
                 stacked_mat = [stacked_mat;stacked_info];
                continue
            else
                new = [0];
                nucleus_lost = [nucleus_info new new new new new new new];
                stacked_mat = [stacked_mat;nucleus_lost];
                continue
        
            end
        end
        StackedNucleiArray = [StackedNucleiArray; stacked_mat];
    end
end

%% Remove nuclei that map to two places
for n = 1:length(StackedNucleiArray)
    stackednucs = StackedNucleiArray{n};
    checkdups = stackednucs(:,NextNuc_index);
    [~, uniqueIdx] = unique( checkdups );
    duplicateLocations = ismember( checkdups, checkdups( setdiff( 1:numel(checkdups), uniqueIdx ) ) );
    dupinds = find(duplicateLocations);
    for x = 1:length(dupinds)
        dupindex = dupinds(x,1);
        stackednucs(dupindex,NextNuc_index) = 0;
    end
    StackedNucleiArray{n} = stackednucs;
end


%% Track nuclei through frames

%initialize arrays
StackedNucleusArray = {};
StackedAreaArray = {};
StackedMedArray = {};
StackedXArray = {};
StackedYArray = {};
StackedDotCountArray = {};
StackedDotIntArray = {};

%separates out data for specific time point
for tt = 1:NumTimepoints;
    timepoint_tracked_nucleus_array = StackedNucleiArray((tt*(NumZs-1)-(NumZs-2)):(tt*(NumZs-1)));
    Stacked_Nuclei = [];
    Stacked_Area = [];
    Stacked_Med = [];
    Stacked_X = [];
    Stacked_Y = [];
    Stacked_DotCount = [];
    Stacked_DotInt = [];
    %go through each z in this time point
    for zz = 1:NumZs-1;
        nuclei = timepoint_tracked_nucleus_array{zz};
        [num_nuclei num_data] = size(nuclei);
        %go through each nucleus at this z
        for nn = 1:num_nuclei;
            Nuc_tracked = zeros(1,zz-1);
            Area_tracked = zeros(1,zz-1);
            Med_tracked = zeros(1,zz-1);
            X_tracked = zeros(1,zz-1);
            Y_tracked = zeros(1,zz-1);
            DotCount_tracked = zeros(1,zz-1);
            DotInt_tracked = zeros(1,zz-1);
            starting_nuc = nuclei(nn,ObjectNum_index);
            current_nuc = starting_nuc;
            if starting_nuc(1,1) == 0; %if 0 then it has already been tracked starting from a previous z
                continue
            else
                for xx = zz:NumZs-1; %go through each z for each nucleus
                    originalinfo = StackedNucleiArray(((tt*(NumZs-1))-(NumZs-2)):(tt*(NumZs-1)));
                     nuclei_info = originalinfo{xx}(find( ...
                        originalinfo{xx}(:,ObjectNum_index)==current_nuc),:);
                    current_nuc = nuclei_info(1,ObjectNum_index);
                    current_Area = nuclei_info(1,Area_index);
                    current_Med = nuclei_info(1,MedIntensity_index);
                    current_X = nuclei_info(1,X_index);
                    current_Y = nuclei_info(1,Y_index);
                    current_DotCount = nuclei_info(1,DotCount_index);
                    current_DotInt = nuclei_info(1,DotInt_index);
                    next_nuc = nuclei_info(1,NextNuc_index);
                    if next_nuc(1,1) == 0;
                        nuc_lost = zeros(1,NumZs-xx);
                        Nuc_tracked = [Nuc_tracked current_nuc nuc_lost];
                        Area_tracked = [Area_tracked current_Area nuc_lost];
                        Med_tracked = [Med_tracked current_Med nuc_lost];
                        X_tracked = [X_tracked current_X nuc_lost];
                        Y_tracked = [Y_tracked current_Y nuc_lost];
                        DotCount_tracked = [DotCount_tracked current_DotCount nuc_lost];
                        DotInt_tracked = [DotInt_tracked current_DotInt nuc_lost];
                        timepoint_tracked_nucleus_array{xx}(find( ...
                            timepoint_tracked_nucleus_array{xx}(:,ObjectNum_index)==current_nuc(1:1)),:) = 0;
                        break 
                    else
                        if xx == NumZs-1;
                            next_Area = nuclei_info(1,NextArea_index);
                            next_Med = nuclei_info(1,NextMed_index);
                            next_X = nuclei_info(1,NextX_index);
                            next_Y = nuclei_info(1,NextY_index);
                            next_DotCount = nuclei_info(1,NextDotCount_index);
                            next_DotInt = nuclei_info(1,NextDotInt_index);
                            Nuc_tracked = [Nuc_tracked current_nuc next_nuc];
                            Area_tracked = [Area_tracked current_Area next_Area];
                            Med_tracked = [Med_tracked current_Med next_Med];
                            X_tracked = [X_tracked current_X next_X];
                            Y_tracked = [Y_tracked current_Y next_Y];
                            DotCount_tracked = [DotCount_tracked current_DotCount next_DotCount];
                            DotInt_tracked = [DotInt_tracked current_DotInt next_DotInt];
                            break
                        else
                            Nuc_tracked = [Nuc_tracked current_nuc];
                            Area_tracked = [Area_tracked current_Area];
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
            Stacked_Nuclei = [Stacked_Nuclei; Nuc_tracked];
            Stacked_Area = [Stacked_Area; Area_tracked];
            Stacked_Med = [Stacked_Med; Med_tracked];
            Stacked_X = [Stacked_X; X_tracked];
            Stacked_Y = [Stacked_Y; Y_tracked];
            Stacked_DotCount = [Stacked_DotCount; DotCount_tracked];
            Stacked_DotInt = [Stacked_DotInt; DotInt_tracked];
        end
    end
    StackedNucleusArray = [StackedNucleusArray; Stacked_Nuclei];
    StackedAreaArray = [StackedAreaArray; Stacked_Area];
    StackedMedArray = [StackedMedArray; Stacked_Med];
    StackedXArray = [StackedXArray; Stacked_X];
    StackedYArray = [StackedYArray; Stacked_Y];
    StackedDotCountArray = [StackedDotCountArray; Stacked_DotCount];
    StackedDotIntArray = [StackedDotIntArray; Stacked_DotInt];
end
%%
end