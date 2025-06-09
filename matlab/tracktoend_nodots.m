%% Find nuclei tracked through entire movie

function [FullyTrackedNuclei,FullyTrackedX, FullyTrackedY] = tracktoend( ...
    EndTime,NumTimepoints, ...
    TrackedYArray,TrackedMedArray,TrackedXArray)

StartTime = EndTime-NumTimepoints+1;
TimeToGradient = 54-StartTime;

TrackedFromStartY = TrackedYArray(TrackedYArray(:,1)>0,:);
TrackedToEndY = TrackedFromStartY(TrackedFromStartY(:,NumTimepoints)>0,:);

TrackedFromStartMed = TrackedMedArray(TrackedYArray(:,1)>0,:);
TrackedToEndMed = TrackedFromStartMed(TrackedFromStartY(:,NumTimepoints)>0,:);

TrackedFromStartX = TrackedXArray(TrackedYArray(:,1)>0,:);
TrackedToEndX = TrackedFromStartX(TrackedFromStartY(:,NumTimepoints)>0,:);

FullyTrackedNuclei = sortrows([TrackedToEndY(:,TimeToGradient) TrackedToEndMed],1);
FullyTrackedX = sortrows([TrackedToEndY(:,TimeToGradient) TrackedToEndX], 1);
FullyTrackedY = sortrows([TrackedToEndY(:,TimeToGradient) TrackedToEndY], 1);

end