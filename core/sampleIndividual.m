function [twin, scales] = sampleIndividual(patient, cfg, truncateSd)
%SAMPLEINDIVIDUAL  Draw one virtual individual from the published
%   between-subject distributions.
%
%   Propofol disposition, effect-site rate constant and potency are drawn
%   from the omega values of Eleveld 2018 (Tables 2 and 3); the
%   haemodynamic baselines from the between-subject variability of
%   Su 2023 (Table 2).
%
%   Those omegas describe a population spanning neonates to morbidly obese
%   adults, so draws are truncated at TRUNCATESD standard deviations
%   (default 1) to keep the virtual cohort inside a plausible adult range.

if nargin < 2, cfg = struct; end
if nargin < 3 || isempty(truncateSd), truncateSd = 1; end

twin = twinInit(patient, cfg);
om   = twin.prop.omega2;
names = {'V1','V2','V3','CL','Q2','Q3','ke0','Ce50'};
scales = struct();
for k = 1:numel(names)
    sd = sqrt(om.(names{k}));
    e  = min(max(randn*sd, -truncateSd*sd), truncateSd*sd);
    scales.(names{k}) = exp(e);
end
twin.prop.V1  = twin.prop.V1  * scales.V1;
twin.prop.V2  = twin.prop.V2  * scales.V2;
twin.prop.V3  = twin.prop.V3  * scales.V3;
twin.prop.CL  = twin.prop.CL  * scales.CL;
twin.prop.Q2  = twin.prop.Q2  * scales.Q2;
twin.prop.Q3  = twin.prop.Q3  * scales.Q3;
twin.prop.ke0 = twin.prop.ke0 * scales.ke0;
twin.L1.C50p  = twin.L1.C50p  * scales.Ce50;

% haemodynamic baselines, CV per cent on a log scale
bsv = twin.L4.bsv;
for f = {'BaseSV','BaseTPR','BaseHR'}
    sd = sqrt(log(1 + (bsv.(f{1})/100)^2));
    e  = min(max(randn*sd, -truncateSd*sd), truncateSd*sd);
    twin.L4.(f{1}) = twin.L4.(f{1}) * exp(e);
    scales.(f{1})  = exp(e);
end

twin.x0  = twinInitialState(twin);
twin.cNoreBase = twin.x0(twin.ix.noreA1)/twin.nore.V1;
o        = twinOutputs(0, twin.x0, twin, 0, false);
twin.CO0 = o.CO;
end
