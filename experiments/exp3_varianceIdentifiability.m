function res = exp3_varianceIdentifiability(nPatients, doPlot)
%EXP3_VARIANCEIDENTIFIABILITY  Demonstration 3: variance propagation and
%   practical identifiability.
%
%   Part one propagates the published between-subject distributions
%   through the twin for a cohort receiving an identical regimen, and
%   separates the pharmacokinetic from the pharmacodynamic contribution by
%   re-running the cohort with only one of the two varied.
%
%   Part two computes the Fisher information of an induction-and-
%   maintenance sequence for the potency, the effect-site rate constant
%   and clearance, observed through the bispectral index alone and through
%   the bispectral index together with mean arterial pressure.

if nargin < 1 || isempty(nPatients), nPatients = 300; end
if nargin < 2, doPlot = true; end
setupPaths();

p        = Patient(1, 45, 172, 75, 'M');
dt       = 10/60;
duration = 30;
truncate = 1.5;

fprintf('\nDemonstration 3, variance propagation (n = %d per cohort)\n', nPatients);
rng(11);  full = cohort(p, nPatients, duration, dt, true,  true,  truncate);
rng(12);  pkOnly = cohort(p, nPatients, duration, dt, true,  false, truncate);
rng(13);  pdOnly = cohort(p, nPatients, duration, dt, false, true,  truncate);

t   = (0:size(full,2)-1)'*dt;
q   = @(M) prctile(M, [2.5 50 97.5], 1);
res.t = t;  res.bandFull = q(full);  res.bandPK = q(pkOnly);  res.bandPD = q(pdOnly);
sel = t > 10;
res.widthFull = mean(res.bandFull(3,sel) - res.bandFull(1,sel));
res.widthPK   = mean(res.bandPK(3,sel)   - res.bandPK(1,sel));
res.widthPD   = mean(res.bandPD(3,sel)   - res.bandPD(1,sel));
k20 = find(t >= 20, 1);
res.bisAt20 = res.bandFull(:,k20)';

fprintf('  BIS at 20 min: median %.0f, 95%% interval %.0f to %.0f\n', ...
        res.bisAt20(2), res.bisAt20(1), res.bisAt20(3));
fprintf('  mean 95%% band width: both %.0f, kinetics only %.0f, potency only %.0f\n', ...
        res.widthFull, res.widthPK, res.widthPD);

fprintf('\nPractical identifiability (Fisher information)\n');
res.fisher = fisherAnalysis(p);
fprintf('  %-10s %14s %16s\n','parameter','BIS only (%)','BIS and MAP (%)');
for k = 1:3
    fprintf('  %-10s %14.1f %16.1f\n', res.fisher.parameters{k}, ...
            res.fisher.rseBisOnly(k), res.fisher.rseBisAndMap(k));
end
fprintf('  max |correlation|  %8.3f %16.3f\n', ...
        res.fisher.maxAbsCorrBisOnly, res.fisher.maxAbsCorrBisAndMap);
fprintf('  condition number   %8.0f %16.0f\n', ...
        res.fisher.condBisOnly, res.fisher.condBisAndMap);

if doPlot, plotVariance(res); end
end

% -------------------------------------------------------------------------
function M = cohort(p, n, duration, dt, varyPK, varyPD, truncate)
M = [];
for i = 1:n
    [twin, sc] = sampleIndividual(p, struct(), truncate);
    if ~varyPK
        base = twinInit(p);
        twin.prop.V1 = base.prop.V1;  twin.prop.V2 = base.prop.V2;
        twin.prop.V3 = base.prop.V3;  twin.prop.CL = base.prop.CL;
        twin.prop.Q2 = base.prop.Q2;  twin.prop.Q3 = base.prop.Q3;
        twin.prop.ke0 = base.prop.ke0;
    end
    if ~varyPD
        twin.L1.C50p = twin.L1.C50p / sc.Ce50;
    end
    twin.x0 = twinInitialState(twin);
    twin.cNoreBase = twin.x0(twin.ix.noreA1)/twin.nore.V1;
    o = twinOutputs(0, twin.x0, twin, 0, false);  twin.CO0 = o.CO;
    sim = twinSimulate(twin, duration, dt, @regimenCaseNoReversal);
    if all(isfinite(sim.BIS)), M = [M; sim.BIS']; end %#ok<AGROW>
    if mod(i,50) == 0, fprintf('    %d of %d\n', i, n); end
end
end
