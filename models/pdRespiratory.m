function m = pdRespiratory()
%PDRESPIRATORY  Non-steady-state CO2-ventilation model with a published
%   propofol-remifentanil interaction.
%
%   Two sources are combined, and the split matters.
%
%   1. Single-agent non-steady-state dynamics, each with its own CO2
%      effect compartment (ke0CO2), distinct from the BIS effect site:
%      Bouillon TW, et al. A model of the ventilatory depressant potency
%      of remifentanil in the non-steady state.
%      Anesthesiology 2003;99:779-787.
%      Bouillon TW, et al. Mixed-effects modeling of the intrinsic
%      ventilatory depressant potency of propofol in the non-steady state.
%      Anesthesiology 2004;100:240-250.
%
%   2. The joint depression of the two agents, which is synergistic:
%      Nieuwenhuijs DJF, Olofsen E, Romberg RR, Sarton E, Ward D,
%      Engbers F, Vuyk J, Mooren R, Teppema LJ, Dahan A. Response surface
%      modeling of remifentanil-propofol interaction on cardiorespiratory
%      control and bispectral index.
%      Anesthesiology 2003;98:312-322.
%
%      depression = [ (Cr/C50r + Cp/C50p) * I(Q) ]^gamma * 0.5
%      I(Q)       = 1 + (Imax - 1) * 4 * Q * (1 - Q),   Qmax = 0.5
%      Q          = (Cr/C50r) / (Cr/C50r + Cp/C50p)
%      Apnoea when the depression reaches 1.
%
%      C50 values are recovered from the reported lambda, the percentage
%      depression at 1 ng/mL remifentanil or 1 ug/mL propofol, because the
%      authors state that the C50 values themselves are extrapolated:
%          C50 = 1 / (2*lambda)^(1/gamma)
%      Verification: at 1 ug/mL propofol plus 1 ng/mL remifentanil this
%      reproduces the 82% depression of ventilation reported in the paper.

% --- CO2 kinetics, Bouillon 2003/2004 (means of the two analyses where
%     the two drugs give different values for the same quantity)
m.VEbase    = 6.75;      % L/min
m.PaCO2base = 41.6;      % mmHg
m.F         = 4.34;      % gain of the CO2 response
m.kelCO2    = 0.095;     % 1/min
m.ke0CO2prop = 0.95;     % 1/min  (Bouillon 2004)
m.ke0CO2remi = 0.92;     % 1/min  (Bouillon 2003)

% --- Nieuwenhuijs 2003, Table 2, resting minute ventilation
m.gamma      = 0.5;
m.Imax       = 1.9;
m.lambdaRemi = 0.277;    % fraction depression at 1 ng/mL
m.lambdaProp = 0.126;    % fraction depression at 1 ug/mL
m.C50r = 1 / (2*m.lambdaRemi)^(1/m.gamma);   % ng/mL
m.C50p = 1 / (2*m.lambdaProp)^(1/m.gamma);   % ug/mL

% --- Airway obstruction in obstructive sleep apnoea, a separate output
%     Koo BN, et al. Yonsei Med J 2015;56:1408-1414.
%     Ce50 = 1.98 - 0.012*(AHI - 46) ug/mL for obstruction score >= 4
m.obstr = struct('theta1',1.98,'theta2',0.012,'AHIref',46,'gamma',1);

m.source = ['Bouillon 2003/2004 (dynamics) + Nieuwenhuijs 2003 ' ...
            '(interaction) + Koo 2015 (obstruction)'];
end
