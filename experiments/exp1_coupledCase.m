function res = exp1_coupledCase(doPlot)
%EXP1_COUPLEDCASE  Demonstration 1: one coupled anaesthetic.
%
%   A 95-minute case in a 45-year-old, 75 kg man: induction, intubation,
%   incision, a controlled haemorrhage between 45 and 55 min, spontaneous
%   ventilation resumed at 82 min and sugammadex at 85 min.
%
%   The purpose is the couplings, which are absent from any layer
%   simulated alone: cardiac output scaling the clearances, surgical
%   stimulus driving the haemodynamic layer, neuromuscular recovery
%   admitting electromyographic contamination of the bispectral index, and
%   reversal feeding an arousal probability gated by hypnotic depth.

if nargin < 1, doPlot = true; end
setupPaths();
rng(20260915);

p    = Patient(1, 45, 172, 75, 'M');
twin = twinInit(p, struct('emgBias',6));

opts.spontaneousFcn = @(t) t > 82;
opts.bloodLossFcn   = @(t) 0.012*(t >= 45 && t < 55);
sim = twinSimulate(twin, 95, 5/60, @(t,o) regimenCase(t,o,75,85), opts);

sim.arousal = zeros(size(sim.t));
rev = sim.t >= 85;
sim.arousal(rev) = arousalProbability(sim.BIS(rev));

k = @(tt) find(abs(sim.t - tt) == min(abs(sim.t - tt)), 1);
res.minMAP            = min(sim.MAP);
res.minCO             = min(sim.CO);
res.coDropPercent     = 100*(1 - min(sim.CO)/sim.CO(1));
res.maxPaCO2          = max(sim.PaCO2);
res.minBIS            = min(sim.BIS);
res.bisAtReversal     = sim.BIS(k(85));
res.arousalAtReversal = sim.arousal(k(85.5));
res.tofBeforeReversal = sim.TOF(k(84));
res.tofAfterReversal  = sim.TOF(k(87));
res.PtolAtIncision    = sim.Ptol(k(16));
res.sim = sim;

fprintf('\nDemonstration 1, coupled case\n');
fprintf('  BIS range during maintenance   %.0f to %.0f\n', ...
        min(sim.BIS(sim.t>10 & sim.t<78)), max(sim.BIS(sim.t>10 & sim.t<78)));
fprintf('  MAP nadir                      %.1f mmHg\n', res.minMAP);
fprintf('  cardiac output fall            %.0f %%\n', res.coDropPercent);
fprintf('  peak arterial CO2              %.0f mmHg\n', res.maxPaCO2);
fprintf('  TOF 84 min -> 87 min           %.2f -> %.2f\n', ...
        res.tofBeforeReversal, res.tofAfterReversal);
fprintf('  BIS at reversal                %.0f  (arousal probability %.2f)\n', ...
        res.bisAtReversal, res.arousalAtReversal);

if doPlot, plotCase(sim); end
end
