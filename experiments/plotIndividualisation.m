function plotIndividualisation(ex, fname)
%PLOTINDIVIDUALISATION  Three-panel figure of one representative individual.
if nargin < 2
    fname = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                     'figures','fig_estimation.pdf');
end
f = figure('Color','w','Position',[100 100 720 600]);
set(f,'DefaultAxesFontSize',8,'DefaultLineLineWidth',1.1);
cTruth = [0 0 0];  cPop = [0.85 0.60 0.10];  cFit = [0.45 0.20 0.65];
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile; hold on
plot(ex.t, ex.bisMeas, 'Color',[.82 .82 .82]);
plot(ex.t, ex.truth.BIS,'Color',cTruth,'LineWidth',1.2);
plot(ex.t, ex.pop.BIS,'--','Color',cPop,'LineWidth',1);
plot(ex.t, ex.fit.BIS,'Color',cFit,'LineWidth',1);
xline(ex.tFit,':'); ylabel('BIS');
legend({'monitored','individual','population twin','individualised twin'}, ...
       'Box','off','FontSize',7,'Location','northwest');

nexttile; hold on
plot(ex.t, ex.mapMeas,'Color',[.82 .82 .82]);
plot(ex.t, ex.truth.MAP,'Color',cTruth,'LineWidth',1.2);
plot(ex.t, ex.pop.MAP,'--','Color',cPop,'LineWidth',1);
plot(ex.t, ex.fit.MAP,'Color',cFit,'LineWidth',1);
xline(ex.tFit,':'); ylabel('MAP (mmHg)');

nexttile; hold on
plot(ex.t, ex.truth.CePropofol,'Color',cTruth,'LineWidth',1.2);
plot(ex.t, ex.pop.CePropofol,'--','Color',cPop,'LineWidth',1);
plot(ex.t, ex.fit.CePropofol,'Color',cFit,'LineWidth',1);
xline(ex.tFit,':'); ylabel({'C_e propofol','(\mug mL^{-1})'}); xlabel('time (min)');

exportgraphics(f, fname, 'ContentType','vector');
fprintf('  figure written to %s\n', fname);
end
