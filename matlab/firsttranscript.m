function first_transcript = firsttranscript(FullyTrackedNucleiRescaled, CorrectedDots, ...
    NumTimepoints,StartTime,EndTime,delaytime)

first_transcript = [];
[num_tracked_nuclei timepoints_tracked] = size(FullyTrackedNucleiRescaled);

for xx = 1:num_tracked_nuclei %number of nuclei that have been tracked
    counter = 4;
    dist_from_mid = FullyTrackedNucleiRescaled(xx,1);
    %nuc_index = FullyTrackedNucleiRescaled(xx,2);
    for yy = 5:NumTimepoints %start at 5 because there are background dots in the first 5min
        Med_level = FullyTrackedNucleiRescaled(xx,yy-delaytime);
        if CorrectedDots(xx,yy) == 0
            counter = counter + 1;
            if counter == NumTimepoints
                checkMed = sum(FullyTrackedNucleiRescaled(xx,2:counter));
                if checkMed > 0
                    totalMed = trapz(FullyTrackedNucleiRescaled(xx,6:counter-delaytime)); %-1 adds a one minute delay
                    %totalMed = integrateMed(end);
                    %collect_mat = [dist_index counter Med_level totalMed];
                    collect_mat = [counter Med_level totalMed xx dist_from_mid];
                    first_transcript = [first_transcript; collect_mat];
                else
                    totalMed = 0;
                    collect_mat = [counter Med_level totalMed xx dist_from_mid];
                    first_transcript = [first_transcript; collect_mat];
                end
                %break
            else
                continue
            end
        else
            counter = counter + 1;
            checkMed = sum(FullyTrackedNucleiRescaled(xx,2:counter));
            if checkMed > 0
                integrateMed = cumtrapz(FullyTrackedNucleiRescaled(xx,2:counter-delaytime));
                totalMed = integrateMed(end);
                collect_mat = [counter Med_level totalMed xx dist_from_mid];
                first_transcript = [first_transcript; collect_mat];
            else
                totalMed = 0;
                collect_mat = [counter Med_level totalMed xx dist_from_mid];
                first_transcript = [first_transcript; collect_mat];
            end
            break
        end
    end
end

cmap=cbrewer2('Blues',10);

[tracked_nuc tracked_timepoints] = size(FullyTrackedNucleiRescaled);

% figure
% hold on
% title("Med v. Time (color = transcribing or not)")
% for nuc = 1:tracked_nuc;
%     xData = [StartTime:EndTime];
%     yData = FullyTrackedNucleiRescaled(nuc,2:end);
%     if first_transcript(nuc,1) < NumTimepoints
%         continue
%         binshade = 8;
%         plot(xData,yData,'-','color',cmap(binshade,:));
%         %ylim([0 0.12]);
%     else
%         continue
%         binshade = 4;
%         plot(xData,yData,'-','color',cmap(binshade,:));
%         %ylim([0 0.12]);
%     end
% end
% hold off

end