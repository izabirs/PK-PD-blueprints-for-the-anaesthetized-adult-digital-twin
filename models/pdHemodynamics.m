function m = pdHemodynamics(patient)
%PDHEMODYNAMICS  Mechanism-based turnover model of total peripheral
%   resistance, stroke volume and heart rate with baroreflex feedback.
%
%   dTPR/dt = kinTPR * RMAP^(-FB) * (1 + effPropTPR) - kout*TPR*(1 - effRemiTPR)
%   dSV*/dt = kinSV  * RMAP^(-FB) * (1 + effPropSV)  - kout*SV* *(1 - effRemiSV)
%   dHR*/dt = kinHR  * RMAP^(-FB)                    - kout*HR* *(1 - effRemiHR)
%   SV  = (SV* + TDE_SV) * (1 - HR_SV*log(HR/BaseHR));   HR = HR* + TDE_HR
%   MAP = TPR*SV*HR;   CO = SV*HR
%
%   Remifentanil effects use the general pharmacodynamic interaction form
%   of Box 1, equation 5:
%       EFF_REMI = ( SL_REMI + INT_RP*Cp/(EC50_PROP + Cp) ) * Cr
%
%   Sources: Su H, Eleveld DJ, Struys MMRF, Colin PJ. Br J Anaesth
%   2022;128:806-816 (structure); Su H, Koomen JV, Eleveld DJ,
%   Struys MMRF, Colin PJ. Br J Anaesth 2023;131:222-233 (interaction
%   model, Table 2, values below).
%
%   The norepinephrine effect on MAP is applied at the output, see
%   hemoOutputs and pkLiNorepinephrine (Joachim et al. 2024).
%
%   Surgical stimulus enters as a multiplicative perturbation (1 + theta*s)
%   of the production rates of TPR and HR. The magnitude theta = 0.284 is
%   taken from Garraza-Obaldia M, et al. Pharmaceutics 2024;16:1615, where
%   it acts on MAP production as a pulse 1 + theta*exp(-kStim*(t - tStim))
%   with kStim = 0.021 1/min, identified for hysteroscopy only and with no
%   effect on HR. ADAPTATION: here the stimulus is a sustained profile s(t)
%   (see surgicalStimulus), so the decay is not applied, and the
%   perturbation is split between TPR and HR as in Eq. (l4) of the
%   article. The HR term is a specification choice, not a published
%   estimate.

age = patient.age;

m.kout     = 0.072;      % 1/min
m.BaseSV   = 82.2;       % mL
m.BaseHR   = 56.0;       % 1/min
m.BaseTPR  = 0.016;      % mmHg*mL^-1*min
m.FB       = 0.66;       % power of the NEGATIVE feedback
m.HR_SV    = 0.312;      % fixed
m.kTDE     = 0.067;      % 1/min
m.LTDE_SV  = 0.090;
m.LTDE_HR  = 0.12;

m.EC50propTPR = 3.21;    % ug/mL
m.EmaxpropTPR = -0.78;
m.gammaProp   = 1.83;
m.EC50propSV  = 0.44;    % ug/mL
m.EmaxpropSV  = -0.15 * exp(0.033*(age - 35));
m.EC50remiTPR = 4.59;    % ng/mL
m.EmaxremiTPR = -1.0;    % fixed
m.SLremiHR    = 0.033;   % per ng/mL
m.SLremiSV    = 0.058;   % per ng/mL
m.INThr       = -0.12;
m.EC50intHR   = 0.20;    % ug/mL
m.INTtpr      = 1.00;    % 100 per cent
m.INTsv       = -0.21;

% surgical stimulus (magnitude from Garraza-Obaldia 2024, see header)
m.thetaStim = 0.284;

% between-subject variability, CV per cent, with correlations (Table 2)
m.bsv = struct('BaseSV',18.3,'BaseTPR',23.3,'BaseHR',15.67, ...
               'C50propTPR',74.34,'EmaxremiTPR',75.30, ...
               'SLremiHR',6.18,'SLremiSV',9.34);
m.corr = struct('SV_TPR',-0.58,'TPR_HR',-0.65,'SLhr_SLsv',0.57);

m.source = 'Su 2022/2023, Br J Anaesth 128:806-816 and 131:222-233';
end
