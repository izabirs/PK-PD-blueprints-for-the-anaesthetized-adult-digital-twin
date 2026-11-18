function res = exp2_individualisation(nPatients, doPlot)
%EXP2_INDIVIDUALISATION  Demonstration 2: what measurement feedback adds.
%
%   A monitored signal is reported directly by its monitor, so predicting
%   its present value adds nothing. The quantity of interest is the
%   predicted future trajectory.
%
%   Virtual individuals are drawn from the published between-subject
%   distributions, so each differs from the population in all six
%   disposition parameters, the effect-site rate constant, the potency and
%   the three haemodynamic baselines. Each receives an identical regimen.
%   The twin is individualised on the first 20 minutes of bispectral index
%   and mean arterial pressure and then predicts the next 15 minutes.
%
%   Runtime is about 10 s per individual.

if nargin < 1 || isempty(nPatients), nPatients = 25; end
if nargin < 2, doPlot = true; end
setupPaths();
rng(101);

p       = Patient(1, 52, 170, 74, 'M');
dt      = 10/60;
tFit    = 20;
horizon = 15;
duration = tFit + horizon;
noiseSd  = 3;

rows = struct('BIS',{},'MAP',{},'Ce',{});
ce50True = zeros(nPatients,1); ce50Fit = zeros(nPatients,1); ce50Pop = zeros(nPatients,1);
for i = 1:nPatients
    [truth, sc] = sampleIndividual(p, struct(), 1);
    sim   = twinSimulate(truth, duration, dt, @regimenEstimation);
    bis   = sim.BIS + noiseSd*randn(size(sim.BIS));
    mapm  = sim.MAP + noiseSd*randn(size(sim.MAP));

    pop     = twinInit(p);
    simPop  = twinSimulate(pop, duration, dt, @regimenEstimation);

    theta   = mapBayesFit(p, bis, mapm, dt, tFit);
    ce50True(i) = truth.L1.C50p;
    ce50Fit(i)  = pop.L1.C50p * exp(theta(1));
    ce50Pop(i)  = pop.L1.C50p;
    fitTwin = applyTheta(twinInit(p), theta);
    simFit  = twinSimulate(fitTwin, duration, dt, @regimenEstimation);

    k = sim.t >= tFit;
    rows(i).BIS = [metrics(simPop.BIS(k), sim.BIS(k)); ...
                   metrics(simFit.BIS(k), sim.BIS(k))];
    rows(i).MAP = [metrics(simPop.MAP(k), sim.MAP(k)); ...
                   metrics(simFit.MAP(k), sim.MAP(k))];
    rows(i).Ce  = [metrics(simPop.CePropofol(k), sim.CePropofol(k)); ...
                   metrics(simFit.CePropofol(k), sim.CePropofol(k))];
    if i == 1
        res.example = struct('t',sim.t,'truth',sim,'pop',simPop,'fit',simFit, ...
                             'bisMeas',bis,'mapMeas',mapm,'tFit',tFit);
    end
    fprintf('  individual %2d of %d\n', i, nPatients);
end

res.n = nPatients;
res.ce50.fittedAbsErrorPercent     = 100*median(abs(ce50Fit - ce50True)./ce50True);
res.ce50.populationAbsErrorPercent = 100*median(abs(ce50Pop - ce50True)./ce50True);
for sig = {'BIS','MAP','Ce'}
    for f = {'MDAPE','RMSE'}
        popv = arrayfun(@(r) r.(sig{1})(1).(f{1}), rows);
        fitv = arrayfun(@(r) r.(sig{1})(2).(f{1}), rows);
        res.(sig{1}).(f{1}).population    = quartiles(popv);
        res.(sig{1}).(f{1}).individualised = quartiles(fitv);
    end
end

fprintf('\nDemonstration 2, prediction of the next %d minutes (n = %d)\n', ...
        horizon, nPatients);
fprintf('  %-22s %20s %20s\n','quantity','population twin','individualised twin');
show('BIS, MDAPE (%)',  res.BIS.MDAPE);
show('BIS, RMSE (units)',res.BIS.RMSE);
show('MAP, MDAPE (%)',  res.MAP.MDAPE);
show('MAP, RMSE (mmHg)', res.MAP.RMSE);
show('propofol Ce, MDAPE (%)', res.Ce.MDAPE);
fprintf(['  potency recovery, median absolute error: fitted %.1f%%, ' ...
         'population value %.1f%%\n'], res.ce50.fittedAbsErrorPercent, ...
        res.ce50.populationAbsErrorPercent);

if doPlot, plotIndividualisation(res.example); end
end

% -------------------------------------------------------------------------
function m = metrics(pred, truth)
pe = 100*(pred - truth)./max(abs(truth), 1e-6);
m.MDPE = median(pe);  m.MDAPE = median(abs(pe));
m.RMSE = sqrt(mean((pred - truth).^2));
end
function q = quartiles(v)
q = struct('median',median(v),'q1',prctile(v,25),'q3',prctile(v,75));
end
function show(name, s)
fprintf('  %-22s %8.1f [%4.1f, %4.1f] %8.1f [%4.1f, %4.1f]\n', name, ...
        s.population.median, s.population.q1, s.population.q3, ...
        s.individualised.median, s.individualised.q1, s.individualised.q3);
end
