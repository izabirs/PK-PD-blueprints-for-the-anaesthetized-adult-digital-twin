function m = pkEleveldPropofol(patient, opts)
%PKELEVELDPROPOFOL  Propofol PK/PD, general-purpose covariate model.
%
%   m = pkEleveldPropofol(patient) returns a struct with the three-
%   compartment disposition parameters, the effect-site rate constant and
%   the BIS pharmacodynamic parameters for the Patient object PATIENT.
%
%   Source: Eleveld DJ, Colin P, Absalom AR, Struys MMRF.
%   Pharmacokinetic-pharmacodynamic model for propofol for broad
%   application in anaesthesia and sedation.
%   Br J Anaesth 2018;120:942-959.  doi:10.1016/j.bja.2018.01.018
%   Structural parameters from Table 2, PD parameters from Table 3.
%
%   Volumes in L, clearances in L/min, rate constants in 1/min.
%
%   opts.sampling  'arterial' (default) or 'venous'
%   opts.opiates   true (default) if a concomitant opioid is given

if nargin < 2, opts = struct; end
if ~isfield(opts,'sampling'), opts.sampling = 'arterial'; end
if ~isfield(opts,'opiates'),  opts.opiates  = true;       end

age = patient.age;  wgt = patient.weight;  bmi = patient.bmi;
sex = double(upper(char(patient.sex)) == 'M');     % 1 male, 0 female

AGE_ref = 35;  WGT_ref = 70;  HGT_ref = 1.7;  BMI_ref = WGT_ref/HGT_ref^2;

% Table 2, theta_1 .. theta_18
th = [6.2830780766822, 25.5013145036879, 272.8166615043603, ...
      1.7895836588902,  1.7500983738779,   1.1085424008536, ...
      0.191307,        42.2760190602615,   9.0548452392807, ...
     -0.015633,        -0.00285709,       33.5531248778544, ...
     -0.0138166,       68.2767978846832,   2.1002218877899, ...
      1.3042680471360,  1.4189043652084,   0.6805003109141];

faging  = @(x) exp(x*(age - AGE_ref));
fsig    = @(x,c50,g) x.^g ./ (c50.^g + x.^g);
fcentr  = @(x) fsig(x, th(12), 1);

PMA     = age + 40/52;                 % post-menstrual age, years
PMA_ref = (40 + AGE_ref*52)/52;
fCLmat     = fsig(PMA*52,     th(8), th(9));
fCLmat_ref = fsig(PMA_ref*52, th(8), th(9));
fQ3mat     = fsig(PMA*52,     th(14), 1);
fQ3mat_ref = fsig(PMA_ref*52, th(14), 1);
fsal     = alSallamiFFM(sex, wgt, age, bmi);
fsal_ref = alSallamiFFM(1,  WGT_ref, AGE_ref, BMI_ref);
if opts.opiates, fopiate = @(x) exp(x*age); else, fopiate = @(x) 1; end

m.V1 = th(1) * fcentr(wgt)/fcentr(WGT_ref);
m.V2 = th(2) * (wgt/WGT_ref) * faging(th(10));
m.V3 = th(3) * (fsal/fsal_ref) * fopiate(th(13));
m.CL = (sex*th(4) + (1-sex)*th(15)) * (wgt/WGT_ref)^0.75 * ...
       (fCLmat/fCLmat_ref) * fopiate(th(11));
m.Q2 = th(5) * (m.V2/th(2))^0.75 * (1 + th(16)*(1 - fQ3mat));
m.Q3 = th(6) * (m.V3/th(3))^0.75 * fQ3mat/fQ3mat_ref;

if strcmpi(opts.sampling,'venous')
    m.V1  = m.V1 * (1 + th(17)*(1 - fcentr(wgt)));
    m.Q2  = m.Q2 * th(18);
    m.ke0 = 1.24 * (wgt/70)^-0.25;
else
    m.ke0 = 0.146 * (wgt/70)^-0.25;
end

% Table 3, pharmacodynamics (BIS)
m.Ce50        = 3.08 * exp(-0.00635*(age - AGE_ref));   % ug/mL
m.E0          = 93.0;
m.gammaLow    = 1.89;      % Ce < Ce50
m.gammaHigh   = 1.47;      % Ce > Ce50
m.delay       = 15 + exp(0.0517*age);                   % s
m.residualBIS = 8.03;

% Inter-individual variability, reported as variance omega^2 (Tables 2, 3)
m.omega2 = struct('V1',0.610,'V2',0.565,'V3',0.597, ...
                  'CL',0.265,'Q2',0.346,'Q3',0.209, ...
                  'ke0',0.702,'Ce50',0.242);
m.source = 'Eleveld 2018, Br J Anaesth 120:942-959';
end
