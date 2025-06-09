%% Find nuclei tracked through entire movie

function [FullyTrackedNuclei,FullyTrackedX,FullyTrackedY,FullyTrackedDotCount,FullyTrackedDotInt] = tracktoend( ...
    EndTime,NumTimepoints, ...
    TrackedYArray,TrackedMedArray,TrackedXArray,TrackedDotCountArray,TrackedDotIntArray)

StartTime = EndTime-NumTimepoints+1;
TimeToGradient = NumTimepoints-10;

TrackedFromStartY = TrackedYArray(TrackedYArray(:,1)>0,:);
TrackedToEndY = TrackedFromStartY(TrackedFromStartY(:,NumTimepoints)>0,:);

TrackedFromStartMed = TrackedMedArray(TrackedYArray(:,1)>0,:);
TrackedToEndMed = TrackedFromStartMed(TrackedFromStartY(:,NumTimepoints)>0,:);

TrackedFromStartX = TrackedXArray(TrackedYArray(:,1)>0,:);
TrackedToEndX = TrackedFromStartX(TrackedFromStartY(:,NumTimepoints)>0,:);

TrackedFromStartDotCount = TrackedDotCountArray(TrackedYArray(:,1)>0,:);
TrackedToEndDotCount = TrackedFromStartDotCount(TrackedFromStartY(:,NumTimepoints)>0,:);

TrackedFromStartDotInt = TrackedDotIntArray(TrackedYArray(:,1)>0,:);
TrackedToEndDotInt = TrackedFromStartDotInt(TrackedFromStartY(:,NumTimepoints)>0,:);

FullyTrackedNuclei = sortrows([TrackedToEndY(:,TimeToGradient) TrackedToEndMed],1);
FullyTrackedX = sortrows([TrackedToEndY(:,TimeToGradient) TrackedToEndX], 1);
FullyTrackedY = sortrows([TrackedToEndY(:,TimeToGradient) TrackedToEndY], 1);
FullyTrackedDotCount = sortrows([TrackedToEndY(:,TimeToGradient) TrackedToEndDotCount], 1);
FullyTrackedDotInt = sortrows([TrackedToEndY(:,TimeToGradient) TrackedToEndDotInt], 1);

end