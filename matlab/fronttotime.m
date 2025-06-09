function Minutes = fronttotime(CellularizationFront, PercentCellularized, NumTimepoints)

Time = CellularizationFront(:,1);
Percent = CellularizationFront(:,2);

ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline');
opts.SmoothingParam = 0.05;
[fitresult, gof] = fit( Time, Percent, ft, opts );

yfitted = feval(fitresult,Time);
% figure
% hold on
% xlabel("Time (min)");
% ylabel("Percent Cellularized");
% title("Dynamics of Cellularization");
% scatter(Time,Percent);
% plot(Time,yfitted);
% hold off

Minutes = [];
for i = 1:NumTimepoints
    percentcell = PercentCellularized(i)*100;
    tol = 3;
    time = Time(abs(yfitted-percentcell) < tol);
    time = mean(time);
    Minutes = [Minutes time];
end

end