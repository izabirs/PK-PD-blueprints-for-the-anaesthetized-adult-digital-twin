function [theta, info] = mapBayesFit(patient, bisMeas, mapMeas, dtMin, tFitMin, opts)
%MAPBAYESFIT  Bounded maximum a posteriori individualisation of the twin.
%
%   Estimates theta = log-scale factors on [Ce50, clearance, ke0] from the
%   monitored bispectral index and mean arterial pressure recorded so far:
%
%     theta_hat = argmin  sum (BIS_i - BIS_hat)^2/sigmaBIS^2
%                       + sum (MAP_i - MAP_hat)^2/sigmaMAP^2
%                       + sum theta_k^2/omega_k^2
%     subject to |theta_k| <= 2*omega_k
%
%   The prior keeps the estimate inside the region in which the model was
%   identified and the bounds stop it drifting along the ridge created by
%   non-identifiability (see fisherAnalysis).
%
%   This is the same class of method as the bounded optimisation compared
%   by Tu Z, et al. J Pharmacokinet Pharmacodyn 2026;53:24, and as the
%   adaptive schemes of Martin-Mateos I, Mendez Perez JA, Reboso Morales
%   JA, Gomez-Gonzalez JF. Comput Biol Med 2016;75:173-180.
%
%   No optimisation toolbox is required: the bounds are imposed by a
%   tanh transformation and the search uses fminsearch.

if nargin < 6, opts = struct; end
g = @(f,v) ternary(isfield(opts,f), @() opts.(f), @() v);
sigmaBIS = g('sigmaBIS',3);  sigmaMAP = g('sigmaMAP',3);
regimen  = g('regimen', @regimenEstimation);
cfg      = g('cfg', struct());
dtFit    = g('dtFit', 20/60);        % integration step used while fitting
nStarts  = g('nStarts', 2);
maxFun   = g('maxFunEvals', 250);

om = pkEleveldPropofol(patient).omega2;
priorSd = [sqrt(om.Ce50), sqrt(om.CL), sqrt(om.ke0)];
bound   = 2*priorSd;

% measurement grid and the (possibly coarser) grid used while fitting
tMeas = (0:numel(bisMeas)-1)'*dtMin;
tFit  = (0:round(tFitMin/dtFit))'*dtFit;
bisOn = interp1(tMeas, bisMeas, tFit, 'linear', 'extrap');
mapOn = interp1(tMeas, mapMeas, tFit, 'linear', 'extrap');

obj = @(z) objective(toTheta(z, bound));
starts = {[0 0 0], [0.4 -0.3 0.2], [-0.4 0.3 -0.2]};
best = [];  bestVal = inf;
for iStart = 1:min(nStarts, numel(starts))
    z0 = fromTheta(starts{iStart}, bound);
    [z, val] = fminsearch(obj, z0, optimset('Display','off', ...
                          'MaxIter',maxFun,'MaxFunEvals',maxFun,'TolX',5e-3));
    if val < bestVal, bestVal = val; best = z; end
end
theta = toTheta(best, bound);
info.objective = bestVal;
info.priorSd   = priorSd;
info.bounds    = bound;

    function r = objective(th)
        twin = twinInit(patient, cfg);
        twin = applyTheta(twin, th);
        sim  = twinSimulate(twin, tFitMin, dtFit, regimen);
        k    = min(numel(bisOn), numel(sim.BIS));
        r = sum((sim.BIS(1:k) - bisOn(1:k)).^2)/sigmaBIS^2 ...
          + sum((sim.MAP(1:k) - mapOn(1:k)).^2)/sigmaMAP^2 ...
          + sum((th./priorSd).^2);
    end
end

function th = toTheta(z, bound), th = bound .* tanh(z); end
function z  = fromTheta(th, bound), z = atanh(min(max(th./bound,-0.99),0.99)); end
function out = ternary(c,a,b), if c, out = a(); else, out = b(); end, end
