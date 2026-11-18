function plotCase(sim, fname)
%PLOTCASE  Seven-panel figure of one coupled case, one row per layer.
if nargin < 2
    fname = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                     'figures','fig_case.pdf');
end
events = [2 15 45 82 85];
labels = {'intubation','incision','haemorrhage','spontaneous vent.','sugammadex'};

f = figure('Color','w','Position',[100 100 780 880]);
tl = tiledlayout(7,1,'TileSpacing','tight','Padding','compact');
ax = gobjects(7,1);
set(f,'DefaultAxesFontSize',8,'DefaultLineLineWidth',1.1);

ax(1) = nexttile;
area(sim.t, sim.stimulus,'FaceColor',[.88 .88 .88],'EdgeColor',[.45 .45 .45]);
ylabel({'stimulus','s(t)'}); ylim([0 1.05]);
yLab = [1.30 1.30 1.30 1.55 1.30];        % stagger the two late labels
for i = 1:numel(events)
    text(events(i), yLab(i), labels{i}, 'FontSize',6.5,'Color',[.35 .35 .35], ...
         'HorizontalAlignment','center','Clipping','off');
end

ax(2) = nexttile; hold on
plot(sim.t, sim.CePropofol); plot(sim.t, sim.CeRemifentanil/5);
ylabel({'effect site'});
lg2 = legend({'propofol C_e (\mug mL^{-1})','remifentanil C_e/5 (ng mL^{-1})'}, ...
       'Box','off','FontSize',6.5,'Location','north','Orientation','horizontal');
lg2.Position(1) = lg2.Position(1) - 0.055;

ax(3) = nexttile; hold on
fill([sim.t(1) sim.t(end) sim.t(end) sim.t(1)],[40 40 60 60], ...
     [0.2 0.4 0.8],'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
plot(sim.t, sim.BIS); ylabel('L1: BIS'); ylim([0 100]);

ax(4) = nexttile;
plot(sim.t, sim.Ptol,'Color',[0.85 0.33 0.10]);
ylabel('L2: P_{tol}'); ylim([0 1.08]);

ax(5) = nexttile;
yyaxis left;  plot(sim.t, sim.PaCO2); ylabel({'L3: P_aCO_2','(mmHg)'});
yyaxis right; plot(sim.t, sim.VE);    ylabel({'V_E','(L min^{-1})'});

ax(6) = nexttile; hold on
plot(sim.t, sim.MAP); plot(sim.t, sim.HR); plot(sim.t, 10*sim.CO,'--');
ylabel('L4'); ylim([20 108]);
lg4 = legend({'MAP (mmHg)','HR (min^{-1})','CO\times10 (L min^{-1})'}, ...
       'Box','off','FontSize',6.5,'Location','north','Orientation','horizontal');
lg4.Position(2) = lg4.Position(2) + 0.006;

ax(7) = nexttile; hold on
plot(sim.t, sim.TOF,'Color',[.35 .35 .35]);
if isfield(sim,'arousal'), plot(sim.t, sim.arousal,'--','Color',[0.47 0.67 0.19]); end
ylabel('L5: TOF'); ylim([0 1.08]); xlabel('time (min)');
lg5 = legend({'TOF ratio','P(arousal)'},'Box','off','FontSize',6.5, ...
       'Location','northwest','Orientation','horizontal');
lg5.Position(1) = lg5.Position(1) + 0.075;

for k = 1:7
    hold(ax(k),'on');
    for i = 1:numel(events)
        xline(ax(k), events(i), '-', 'Color',[.8 .8 .8],'HandleVisibility','off');
    end
    xlim(ax(k), [0 sim.t(end)]);
    box(ax(k),'off');
    if k < 7, set(ax(k),'XTickLabel',[]); end
end
linkaxes(ax,'x');

exportgraphics(f, fname, 'ContentType','vector');
fprintf('  figure written to %s\n', fname);
end
