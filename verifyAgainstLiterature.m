function results = verifyAgainstLiterature(verbose)
%VERIFYAGAINSTLITERATURE  Check that each implemented layer reproduces a
%   number reported in its source publication.
%
%   Run this after any change to the model files. Every check states the
%   published value, the value the implementation produces, and the
%   tolerance applied. Two kinds of check are included and labelled:
%     [source]    a value reported by the publication that the code has to
%                 reproduce from the covariate equations (Eleveld reference
%                 individual, Li clearance reduction);
%     [structure] a property the implementation must have by construction
%                 (fixed point of the turnover model, reproduction of a
%                 parameter used to calibrate a term, evaluation of a
%                 published regression). These guard against coding errors
%                 but are not independent evidence.

if nargin < 1, verbose = true; end
setupPaths();
results = struct('name',{},'published',{},'obtained',{},'tol',{},'pass',{});

add = @(r,n,p,o,t) [r, struct('name',n,'published',p,'obtained',o,'tol',t, ...
                              'pass',abs(o-p) <= t)];

% --- 1. Eleveld 2018, reference individual (35 yr, 70 kg, 170 cm, male)
ref  = Patient(0, 35, 170, 70, 'M');
prop = pkEleveldPropofol(ref, struct('opiates',false));
results = add(results,'[source] Eleveld V1ref (L)',            6.28,  prop.V1,   0.02);
results = add(results,'[source] Eleveld V2ref (L)',           25.5,   prop.V2,   0.1);
results = add(results,'[source] Eleveld V3ref (L)',          273,     prop.V3,   1.0);
results = add(results,'[source] Eleveld CLref (L/min)',        1.79,  prop.CL,   0.01);
results = add(results,'[source] Eleveld ke0 (1/min)',          0.146, prop.ke0,  0.001);
results = add(results,'[source] Eleveld Ce50 (ug/mL)',         3.08,  prop.Ce50, 0.01);

% --- 2. Nieuwenhuijs 2003, single-drug depression of resting ventilation
p    = Patient(1, 45, 172, 75, 'M');
twin = twinInit(p);
veBase = twin.L3.VEbase;
depRemi = 1 - ventilation(twin.L3.PaCO2base, 0, 1, twin)/veBase;   % 1 ng/mL
depProp = 1 - ventilation(twin.L3.PaCO2base, 1, 0, twin)/veBase;   % 1 ug/mL
results = add(results,'[structure] Nieuwenhuijs remi (%)', 27.7, 100*depRemi, 0.5);
results = add(results,'[structure] Nieuwenhuijs prop (%)',     12.6, 100*depProp, 0.5);

% --- 3. Su 2023, structural checks on the haemodynamic layer.
%     The baseline must be an exact fixed point of the turnover equations,
%     and the baseline MAP must be the product of the three baselines.
d0 = hemoDerivativesPublic(twin.x0, twin.L4, 0, 0, 0, 0);
results = add(results,'[structure] Su 2023 fixed point', 0, max(abs(d0)), 1e-12);
o0 = hemoOutputs(twin.x0, twin.L4);
mapExpected = twin.L4.BaseTPR*twin.L4.BaseSV*(1+twin.L4.LTDE_SV)* ...
              twin.L4.BaseHR*(1+twin.L4.LTDE_HR);
results = add(results,'[structure] Su 2023 baseline MAP (mmHg)', mapExpected, o0.MAP, 1e-9);

%     Informational: the Discussion of Su 2023 quotes a MAP fall of 11.7%
%     and an HR fall of 2.7% at 1 ug/mL propofol with 1 ng/mL
%     remifentanil, obtained from a simulation of the full model rather
%     than from an infinite-time steady state. The steady state of the
%     equations as published is reported here for comparison only.
base = steadyStateHemo(twin, 0, 0);
comb = steadyStateHemo(twin, 1, 1);
infoMAP = 100*(comb.MAP - base.MAP)/base.MAP;
infoHR  = 100*(comb.HR  - base.HR )/base.HR;

% --- 4. Kleijn 2011, TOF ratio at the reported EC50
tofAtEC50 = 1 - hill(twin.nmb.roc.EC50, twin.nmb.roc.EC50, twin.nmb.roc.gamma);
results = add(results,'[structure] Kleijn TOF at EC50', 0.5, tofAtEC50, 1e-6);

% --- 5. Li 2024, clearance reduction at the study median propofol
red = 100*(1 - exp(twin.nore.thetaProp/100 * 3.53));
results = add(results,'[source] Li 2024 NE CL reduction (%)', 11.8, red, 0.3);

% --- 6. Kim 2023, arousal probability at the midpoint of the light arm
results = add(results,'[structure] Kim 2023 arousal at BIS 50', 0.364, arousalProbability(50), 0.01);

% --- 7. Perrin 2022, fitted regression for the personalised remifentanil Ce
results = add(results,'[structure] Perrin 2022 Ce, dNOL 20 (ng/mL)', ...
              0.069*20 + 3.1, perrinRemifentanilCe(20), 1e-9);

if verbose
    fprintf('\n%-44s %10s %10s %6s\n','check','published','obtained','pass');
    for k = 1:numel(results)
        r = results(k);
        fprintf('%-44s %10.3f %10.3f %6s\n', r.name, r.published, r.obtained, ...
                string(r.pass));
    end
    fprintf('\n%d of %d checks passed.\n', sum([results.pass]), numel(results));
    fprintf(['\nInformational (Su 2023 Discussion quotes -11.7%% MAP and -2.7%% HR\n' ...
             'from a full-model simulation): steady state of the published\n' ...
             'equations gives %.1f%% MAP and %+.1f%% HR.\n'], infoMAP, infoHR);
end
end

% -------------------------------------------------------------------------
function o = steadyStateHemo(twin, cpProp, cpRemi)
%STEADYSTATEHEMO  Integrate the haemodynamic layer alone to steady state
%   at fixed plasma concentrations.
ix = twin.ix;
x  = twin.x0;
x(ix.tdeSV) = 0;  x(ix.tdeHR) = 0;        % transient terms fully decayed
dt = 0.05;
for k = 1:round(120/dt)
    d = zeros(ix.n,1);
    d(ix.TPR:ix.tdeHR) = hemoDerivativesPublic(x, twin.L4, cpProp, cpRemi, 0, 0);
    x = x + dt*d;
end
o = hemoOutputs(x, twin.L4);
end
