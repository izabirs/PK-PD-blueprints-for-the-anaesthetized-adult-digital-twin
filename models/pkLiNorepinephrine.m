function m = pkLiNorepinephrine(patient)
%PKLINOREPINEPHRINE  Norepinephrine PK with a propofol effect on clearance.
%
%   Source: Li Y, Koomen JV, Eleveld DJ, van den Berg JP, Vos JJ,
%   de Keijzer IN, Struys MMRF, Colin P. Population pharmacokinetic
%   modelling of norepinephrine in healthy volunteers prior to and during
%   general anesthesia. Clin Pharmacokinet 2024;63:1597-1608.
%
%   Parameters from Table 3: CL 2.1 L/min, Vc 2.4 L, Q 0.6 L/min, Vp 3.6 L,
%   endogenous production 497.7 ng/min, AGE~CL -0.344, CPROP~CL -3.57.
%   Covariate functions are exponentials of theta/100, so the propofol
%   coefficient of -3.57 is a percentage per ug/mL. The reference state is
%   the awake patient (Cprop = 0); at the study median of 3.53 ug/mL this
%   gives the 11.8% reduction in clearance reported in the paper.
%
%   Production is reported in ng/min and is converted here to nmol/min
%   with the molar mass of norepinephrine (169.18 g/mol). The resulting
%   endogenous steady state, prod/CL, is 1.4 nmol/L for a 70 kg adult
%   aged 25, the baseline reported by Li et al. (Table 2, 18-34 years).
%
%   ASSUMPTION: the published model has no weight covariate. Clearances,
%   volumes and production are scaled allometrically to body weight
%   (exponents 0.75 and 1), which leaves the 70 kg reference unchanged.
%
%   Concentrations are handled in nmol/L throughout, the unit of the
%   intraoperative PD model of Joachim et al. (Br J Clin Pharmacol
%   2024;90:2861-2869). Volumes in L, clearances in L/min, production in
%   nmol/min.

wgt = patient.weight;  age = patient.age;

m.CL       = 2.1 * exp(-0.344/100*(age - 35)) * (wgt/70)^0.75;
m.V1       = 2.4 * (wgt/70);
m.Q        = 0.6 * (wgt/70)^0.75;
m.V2       = 3.6 * (wgt/70);
m.MW       = 169.18;                          % g/mol
m.uEndo    = 497.7/m.MW * (wgt/70)^0.75;      % 497.7 ng/min -> nmol/min

% Pharmacodynamics: fractional rise in MAP above the endogenous baseline.
% Joachim J, Cartailler J, Vallee F, et al. Br J Clin Pharmacol
% 2024;90:2861-2869: E0 57.09 mmHg, Emax 113.52 mmHg, EC50 15.7 nmol/L,
% gamma 1.46. Referencing the effect to the endogenous concentration
% avoids counting the resting vasomotor tone twice, since the haemodynamic
% baselines of the turnover model already include it.
m.pd = struct('E0',57.09,'Emax',113.52,'EC50',15.7,'gamma',1.46);
m.pd.fracMax = (m.pd.Emax - m.pd.E0)/m.pd.E0;
m.thetaProp = -3.57;                          % percent per ug/mL propofol
m.source   = 'Li 2024, Clin Pharmacokinet 63:1597-1608';
end
