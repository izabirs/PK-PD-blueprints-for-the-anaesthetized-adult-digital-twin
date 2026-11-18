function m = pkKleijnRocSugammadex(patient)
%PKKLEIJNROCSUGAMMADEX  Rocuronium and sugammadex PK/PD with 1:1 binding.
%
%   Source: Kleijn HJ, Zollinger DP, van den Heuvel MW, Kerbusch T.
%   Population pharmacokinetic-pharmacodynamic analysis for
%   sugammadex-mediated reversal of rocuronium-induced neuromuscular
%   blockade. Br J Clin Pharmacol 2011;72:415-433.
%
%   Amounts are handled in micromoles so that the encapsulation reaction
%   can be written directly. Volumes in L, clearances in L/min.
%   Molar masses: rocuronium 529.8 g/mol, sugammadex 2178 g/mol.
%
%   Not represented: the sevoflurane effects on ke0, EC50 and ks, the
%   Asian-race effects on rocuronium Q and sugammadex V1, and distribution
%   of the complex to the peripheral compartment.

bw  = patient.weight;
age = patient.age;
cr  = 119;                 % creatinine clearance, mL/min: model reference value

% Rocuronium
m.roc.CL = 0.269 * (bw/70)^0.75 * (1 - 0.00678*(age - 43));
m.roc.V1 = 4.73  * (bw/70) * exp(-0.00143*(cr - 119));
m.roc.Q  = 0.279 * (bw/70)^0.75;
m.roc.V2 = 6.76  * (bw/70) * exp(0.00613*(age - 43));
m.roc.ke0    = 0.134;      % 1/min
m.roc.EC50   = 1.62;       % umol/L
m.roc.gamma  = 7.52;
m.roc.E0     = 1.0;        % T4/T1 normalised to baseline (Kleijn E0 = 1.04)

% Sugammadex
m.sug.CL = 0.093 * (1 + 0.00378*(bw - 74.5)) * ((2*cr/(cr + 119))^1.29);
m.sug.V1 = 4.42  * (bw/70) * (1 - 0.00354*(bw - 74.5));
m.sug.Q  = 0.206 * (bw/70)^0.75;
m.sug.V2 = 6.35  * (bw/70) * exp(-0.00305*(cr - 119));

% 1:1 complexation in the central compartment
m.bind.k1 = 0.61;                      % 1/min per umol/L, association
m.bind.k2 = 0.034;                     % 1/min, dissociation
m.bind.Kd = m.bind.k2 / m.bind.k1;     % 0.0559 umol/L
% additional sugammadex-driven removal from the effect compartment
m.bind.ks = 0.033 * (bw/70)^-0.25;     % 1/min per umol/L

m.MW_roc = 529.8;  m.MW_sug = 2178;
m.source = 'Kleijn 2011, Br J Clin Pharmacol 72:415-433';
end
