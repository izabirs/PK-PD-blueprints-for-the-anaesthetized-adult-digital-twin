function o = twinOutputs(tMin, x, twin, stimulus, spontaneous, noiseSd)
%TWINOUTPUTS  Output map g(x,u,theta) of the coupled twin.
%
%   Returns the monitored signals: BIS, probability of tolerance, minute
%   ventilation and arterial CO2, mean arterial pressure, heart rate,
%   cardiac output and the train-of-four ratio.

if nargin < 4 || isempty(stimulus),   stimulus   = surgicalStimulus(tMin); end
if nargin < 5 || isempty(spontaneous), spontaneous = false; end
if nargin < 6, noiseSd = 0; end

ix = twin.ix;
cePro = x(ix.propCe);   ceRem = x(ix.remiCe);
cNore = x(ix.noreA1)/twin.nore.V1;                 % nmol/L

% ---- Layer 5, train of four
tof = 1 - hill(x(ix.rocCe), twin.nmb.roc.EC50, twin.nmb.roc.gamma);
tof = min(max(tof,0),1);

% ---- Layer 1, BIS as an observation of the latent state
conf = twin.confounders;
conf.bEMG = twin.cfg.emgBias;
bis = obsBis(x(ix.hBis), twin.L1, stimulus, tof, conf, noiseSd);

% ---- Layer 2, probability of tolerating the prevailing stimulus
o.Ptol = nociceptionTolerance(cePro, ceRem, stimulus, twin.L2);

% ---- Layer 4, haemodynamics
hemo = hemoOutputs(x, twin.L4, twin.nore, cNore, twin.cNoreBase);
mapv = hemo.MAP;
if noiseSd > 0, mapv = mapv + noiseSd*randn; end

% ---- Layer 3, respiratory
if spontaneous
    ve = ventilation(x(ix.PaCO2), x(ix.ceCO2prop), x(ix.ceCO2remi), twin);
else
    ve = NaN;
end

o.t = tMin;  o.BIS = bis;  o.MAP = mapv;  o.HR = hemo.HR;
o.SV = hemo.SV;  o.TPR = hemo.TPR;  o.CO = hemo.CO;
o.TOF = tof;  o.VE = ve;  o.PaCO2 = x(ix.PaCO2);
o.CePropofol = cePro;  o.CeRemifentanil = ceRem;  o.CeRocuronium = x(ix.rocCe);
o.CNorepinephrine = cNore;  o.stimulus = stimulus;
o.hypnoticDepth = x(ix.hBis);
end
