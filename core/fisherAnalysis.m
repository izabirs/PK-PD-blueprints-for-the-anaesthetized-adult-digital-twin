function res = fisherAnalysis(patient, opts)
%FISHERANALYSIS  Practical identifiability of (Ce50, ke0, clearance) from
%   the monitored signals, through the Fisher information of a given
%   induction-and-maintenance sequence.
%
%   Sensitivities are computed by finite differences of the output
%   trajectories with respect to the log-parameters, and the information
%   matrix is formed as sum_signals S'S/sigma^2. Reported are the relative
%   standard errors, the largest absolute parameter correlation and the
%   condition number of the information matrix, for the bispectral index
%   alone and for the bispectral index together with mean arterial
%   pressure.
%
%   The general result that common pharmacodynamic structures are not
%   structurally identifiable in their natural parameterisation is set out
%   by Janzen DLI, et al. Front Physiol 2016;7:590.

if nargin < 2, opts = struct; end
g = @(f,v) ternary(isfield(opts,f), @() opts.(f), @() v);
dt       = g('dtMin', 10/60);
duration = g('durationMin', 20);
sigmaBIS = g('sigmaBIS', 3);
sigmaMAP = g('sigmaMAP', 3);
regimen  = g('regimen', @regimenCaseNoReversal);
step     = g('step', 0.02);

[b0, m0] = run([0 0 0]);
S_bis = zeros(numel(b0),3);  S_map = zeros(numel(m0),3);
for k = 1:3
    th = zeros(1,3); th(k) = step;
    [b1, m1] = run(th);
    S_bis(:,k) = (b1 - b0)/step;      % sensitivity wrt the log-parameter
    S_map(:,k) = (m1 - m0)/step;
end

res.parameters = {'Ce50','ke0','CL'};
[res.rseBisOnly, cBis]   = stats({S_bis}, sigmaBIS);
[res.rseBisAndMap, cAll] = stats({S_bis, S_map}, [sigmaBIS sigmaMAP]);
res.maxAbsCorrBisOnly   = cBis.maxCorr;
res.maxAbsCorrBisAndMap = cAll.maxCorr;
res.condBisOnly         = cBis.cond;
res.condBisAndMap       = cAll.cond;

    function [bis, mapv] = run(theta)
        twin = twinInit(patient);
        % order here is [Ce50, ke0, CL]; applyTheta expects [Ce50, CL, ke0]
        twin = applyTheta(twin, [theta(1) theta(3) theta(2)]);
        sim  = twinSimulate(twin, duration, dt, regimen);
        bis  = sim.BIS;  mapv = sim.MAP;
    end
end

function [rse, c] = stats(Slist, sigmas)
fim = zeros(3);
for k = 1:numel(Slist)
    S = Slist{k};
    fim = fim + (S'*S)/sigmas(k)^2;
end
cov = inv(fim);
rse = 100*sqrt(diag(cov))';
d   = sqrt(diag(cov));
corr = cov ./ (d*d');
c.maxCorr = max(abs(corr(~eye(3))));
c.cond    = cond(fim);
end

function out = ternary(c,a,b), if c, out = a(); else, out = b(); end, end
