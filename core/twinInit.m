function twin = twinInit(patient, cfg)
%TWININIT  Assemble the coupled twin for one patient.
%
%   twin = twinInit(patient) builds every layer from the published models
%   in ../models and returns a struct that twinDerivatives and twinOutputs
%   consume. PATIENT is a Patient object (see initialization/Patient.m).
%
%   cfg fields, all optional:
%     .coCoupling   true (default): clearances scale with cardiac output
%     .coBounds     [0.6 1.4] (default) clip on the cardiac output ratio
%     .emgBias      BIS units of electromyographic contamination at full
%                   stimulus and full neuromuscular recovery (default 6,
%                   an ASSUMPTION; see README)
%     .respSynergy  true (default): Nieuwenhuijs interaction; false uses a
%                   multiplicative combination of the single-agent terms

if nargin < 2, cfg = struct; end
d = @(f,v) ternary(isfield(cfg,f), @() cfg.(f), @() v);

twin.patient = patient;
twin.prop  = pkEleveldPropofol(patient);
twin.remi  = pkEleveldRemifentanil(patient);
twin.nmb   = pkKleijnRocSugammadex(patient);
twin.nore  = pkLiNorepinephrine(patient);
twin.L1    = pdBisSurface();
twin.L2    = pdNociception();
twin.L3    = pdRespiratory();
twin.L4    = pdHemodynamics(patient);

twin.cfg.coCoupling  = d('coCoupling', true);
twin.cfg.coBounds    = d('coBounds', [0.6 1.4]);
twin.cfg.emgBias     = d('emgBias', 6);
twin.cfg.respSynergy = d('respSynergy', true);
twin.cfg.clScale     = d('clScale', 1);
twin.confounders     = d('confounders', struct());

twin.ix = twinIndex();
twin.x0 = twinInitialState(twin);
twin.cNoreBase = twin.x0(twin.ix.noreA1)/twin.nore.V1;   % nmol/L, endogenous
o       = twinOutputs(0, twin.x0, twin, 0, false);
twin.CO0 = o.CO;
end

function out = ternary(c, a, b)
if c, out = a(); else, out = b(); end
end
