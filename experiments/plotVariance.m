function plotVariance(res, fname)
%PLOTVARIANCE  Prediction band and identifiability summary.
if nargin < 2
    fname = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                     'figures','fig_variance.pdf');
end
f = figure('Color','w','Position',[100 100 860 320]);
set(f,'DefaultAxesFontSize',8,'DefaultLineLineWidth',1.1);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile; hold on
fill([res.t; flipud(res.t)], [res.bandFull(1,:)'; flipud(res.bandFull(3,:)')], ...
     [0.2 0.4 0.8],'FaceAlpha',0.18,'EdgeColor','none');
fill([res.t; flipud(res.t)], [res.bandPK(1,:)'; flipud(res.bandPK(3,:)')], ...
     [0.9 0.5 0.1],'FaceAlpha',0.25,'EdgeColor','none');
plot(res.t, res.bandFull(2,:),'k','LineWidth',1);
xlabel('time (min)'); ylabel('BIS'); ylim([0 100]);
legend({'kinetics and potency','kinetics only','population prediction'}, ...
       'Box','off','FontSize',7,'Location','northeast');

nexttile;
b = bar(categorical(res.fisher.parameters), ...
        [res.fisher.rseBisOnly(:) res.fisher.rseBisAndMap(:)]);
% linear axis from zero, so that every bar is visible and none is clipped
rse = [res.fisher.rseBisOnly(:); res.fisher.rseBisAndMap(:)];
ylim([0, 1.15*max(rse)]); ylabel('relative standard error (%)');
legend({'BIS only','BIS and MAP'},'Box','off');

exportgraphics(f, fname, 'ContentType','vector');
fprintf('  figure written to %s\n', fname);
end
