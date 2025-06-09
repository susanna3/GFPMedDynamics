%Find the midline and convert Y-coord to distance from midline
%Then corrects for background based on cells far from the midline

function CorrectedNuclei =  findmidline_fixed(NucleiStacked, ...
    X_index,Y_index,MedIntensity_index,DotCount_index, DotInt_index, PixelToMicron,BackgroundCutoffLow,BackgroundCutoffHigh,NumTimepoints)

MidlineOverTime = [];
CorrectedNuclei = {};
CheckRescale = {};

for ii = 1:NumTimepoints
    timepoint_mat = NucleiStacked{ii};
    distancedata = timepoint_mat(:,Y_index);
    maddata = timepoint_mat(:,MedIntensity_index);
    ft = fittype( 'smoothingspline' );
    opts = fitoptions( 'Method', 'SmoothingSpline');
    opts.SmoothingParam = 0.000005;
    [fitresult, gof] = fit( distancedata, maddata, ft, opts );
    yfitted = feval(fitresult,distancedata);
    % figure
    % plot(distancedata,yfitted)
    fittedmed = [yfitted distancedata];
    amplitude = max(yfitted);
    midlineinfo = fittedmed(fittedmed(:,1)==amplitude,:);
    midline = midlineinfo(1,2);
    MidlineOverTime = [MidlineOverTime; midlineinfo];

    rescaley = abs(distancedata-midline)./PixelToMicron;
    CheckRescale = [CheckRescale,rescaley];
    backgroundnuclei = maddata(rescaley(:,1)>BackgroundCutoffLow & rescaley(:,1)<BackgroundCutoffHigh,:);
    background = min(backgroundnuclei);
    [numnuc wid] = size(maddata);
    timepointMad = [];
    for n = 1:numnuc
        dist = rescaley(n);
        oldmad = maddata(n);
        newmad = oldmad-background;
        %newmad=oldmad;
        dotcount = timepoint_mat(n,DotCount_index);
        dotint = timepoint_mat(n,DotInt_index);
        ap = timepoint_mat(n,X_index)./PixelToMicron;
        if newmad < 0
            newmad = 0;
            timepointMad = [timepointMad;dist newmad dotcount dotint ap];
        else
            timepointMad = [timepointMad;dist newmad dotcount dotint ap];
        end
    end
    CorrectedNuclei = [CorrectedNuclei timepointMad];
end

end