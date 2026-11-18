function dx = twinDerivatives(tMin, x, u, twin)
%TWINDERIVATIVES  Right-hand side f(x,u,theta) of the coupled system.
%
%   u is a struct of input rates:
%     .propofol        mg/min
%     .remifentanil    ug/min
%     .rocuronium      mg/min
%     .sugammadex      mg/min
%     .norepinephrine  ug/min (converted internally to nmol/min)
%     .stimulus        surgical stimulus in [0,1]
%     .spontaneous     logical, spontaneous ventilation
%     .bloodLoss       fractional stroke-volume loss per minute (optional)
%
%   Three cross-layer couplings are encoded here:
%     1. cardiac output scales every clearance (Bienert 2020)
%     2. surgical stimulus multiplies the haemodynamic production rates
%        (Garraza-Obaldia 2024)
%     3. neuromuscular reversal feeds hypnosis, handled at the output
%        (see arousalProbability) and through the BIS observation model

ix = twin.ix;
dx = zeros(ix.n,1);
g  = @(f,v) ternary(isfield(u,f), @() u.(f), @() v);

rProp = g('propofol',0);  rRemi = g('remifentanil',0);
rRoc  = g('rocuronium',0); rSug  = g('sugammadex',0);
rNore = g('norepinephrine',0);
s     = g('stimulus', surgicalStimulus(tMin));
spont = g('spontaneous', false);
bleed = g('bloodLoss', 0);

cpProp = x(ix.propA1)/twin.prop.V1;        % ug/mL
cpRemi = x(ix.remiA1)/twin.remi.V1;        % ng/mL
cNore  = x(ix.noreA1)/twin.nore.V1;        % nmol/L

% ---- coupling 1: cardiac output scales the clearances
hemo = hemoOutputs(x, twin.L4, twin.nore, cNore, twin.cNoreBase);
clScale = twin.cfg.clScale;
if twin.cfg.coCoupling && isfield(twin,'CO0') && twin.CO0 > 0
    r = hemo.CO / twin.CO0;
    clScale = clScale * min(max(r, twin.cfg.coBounds(1)), twin.cfg.coBounds(2));
end

% ---- propofol and remifentanil, three-compartment disposition
dx(ix.propA1:ix.propCe) = threeComp(x(ix.propA1:ix.propCe), rProp, twin.prop, clScale);
dx(ix.remiA1:ix.remiCe) = threeComp(x(ix.remiA1:ix.remiCe), rRemi, twin.remi, clScale);

% ---- rocuronium and sugammadex, in umol, with 1:1 complexation
nmb = twin.nmb;
roc1 = x(ix.rocR1);  roc2 = x(ix.rocR2);
sug1 = x(ix.sugS1);  sug2 = x(ix.sugS2);
rs   = x(ix.complexRS);
cR   = roc1/nmb.roc.V1;   cS = sug1/nmb.sug.V1;
bind   = nmb.bind.k1 * cR * cS * nmb.roc.V1;     % umol/min
unbind = nmb.bind.k2 * rs;
dx(ix.rocR1) = rRoc*1000/nmb.MW_roc - nmb.roc.CL*cR - nmb.roc.Q*cR ...
               + nmb.roc.Q*roc2/nmb.roc.V2 - bind + unbind;
dx(ix.rocR2) = nmb.roc.Q*cR - nmb.roc.Q*roc2/nmb.roc.V2;
dx(ix.sugS1) = rSug*1000/nmb.MW_sug - nmb.sug.CL*cS - nmb.sug.Q*cS ...
               + nmb.sug.Q*sug2/nmb.sug.V2 - bind + unbind;
dx(ix.sugS2) = nmb.sug.Q*cS - nmb.sug.Q*sug2/nmb.sug.V2;
dx(ix.complexRS) = bind - unbind - nmb.sug.CL*rs/nmb.sug.V1;
dx(ix.rocCe) = nmb.roc.ke0*(cR - x(ix.rocCe)) - nmb.bind.ks*cS*x(ix.rocCe);

% ---- norepinephrine, with the propofol effect on clearance
clNore = twin.nore.CL * exp(twin.nore.thetaProp/100 * cpProp);
k10 = clNore/twin.nore.V1;  k12 = twin.nore.Q/twin.nore.V1;
k21 = twin.nore.Q/twin.nore.V2;
dx(ix.noreA1) = rNore*1000/twin.nore.MW + twin.nore.uEndo - (k10 + k12)*x(ix.noreA1) + k21*x(ix.noreA2);
dx(ix.noreA2) = k12*x(ix.noreA1) - k21*x(ix.noreA2);

% ---- Layer 4, haemodynamic turnover with baroreflex feedback
dx(ix.TPR:ix.tdeHR) = hemoDerivativesPublic(x, twin.L4, cpProp, cpRemi, s, bleed);

% ---- Layer 3, CO2 effect sites and the CO2 balance
dx(ix.ceCO2prop) = twin.L3.ke0CO2prop*(cpProp - x(ix.ceCO2prop));
dx(ix.ceCO2remi) = twin.L3.ke0CO2remi*(cpRemi - x(ix.ceCO2remi));
if spont
    ve  = ventilation(x(ix.PaCO2), x(ix.ceCO2prop), x(ix.ceCO2remi), twin);
    vst = twin.L3.VEbase/twin.L3.kelCO2;                 % effective CO2 store, L
    vco2 = twin.L3.VEbase*twin.L3.PaCO2base/863;          % L/min
    dx(ix.PaCO2) = (863*vco2 - ve*x(ix.PaCO2))/vst;
else
    dx(ix.PaCO2) = (twin.L3.PaCO2base - x(ix.PaCO2))/2;   % ventilator-controlled
end

% ---- Layer 1, latent hypnotic depth
U = max(x(ix.propCe),0)/twin.L1.C50p + max(x(ix.remiCe),0)/twin.L1.C50r;
if U > 0, hInf = U^twin.L1.gamma/(1 + U^twin.L1.gamma); else, hInf = 0; end
dx(ix.hBis) = (hInf - x(ix.hBis))/(twin.L1.tau/60);
end

% -------------------------------------------------------------------------
function d = threeComp(xs, rate, pk, clScale)
%THREECOMP  [A1 A2 A3 Ce] derivatives for a mammillary model with an
%   effect compartment. CLSCALE implements the cardiac-output coupling.
a1 = xs(1); a2 = xs(2); a3 = xs(3); ce = xs(4);
k10 = pk.CL*clScale/pk.V1;  k12 = pk.Q2*clScale/pk.V1;
k21 = pk.Q2*clScale/pk.V2;  k13 = pk.Q3*clScale/pk.V1;
k31 = pk.Q3*clScale/pk.V3;
d = zeros(4,1);
d(1) = rate - (k10 + k12 + k13)*a1 + k21*a2 + k31*a3;
d(2) = k12*a1 - k21*a2;
d(3) = k13*a1 - k31*a3;
d(4) = pk.ke0*(a1/pk.V1 - ce);
end

% -------------------------------------------------------------------------
function out = ternary(c,a,b)
if c, out = a(); else, out = b(); end
end
