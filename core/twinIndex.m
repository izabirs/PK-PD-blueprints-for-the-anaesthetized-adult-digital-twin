function ix = twinIndex()
%TWININDEX  Index map of the 25-state coupled system.
%
%   Drug quantities are carried as amounts, not concentrations, so that a
%   change of a volume parameter cannot silently rescale a stored state.
%   Rocuronium and sugammadex are carried in micromoles because their
%   interaction is a 1:1 chemical reaction.
%
%   state                         units
%   propA1..A3   propofol         mg
%   propCe       propofol         ug/mL
%   remiA1..A3   remifentanil     ug
%   remiCe       remifentanil     ng/mL
%   rocR1,R2     rocuronium       umol
%   sugS1,S2     sugammadex       umol
%   complexRS    1:1 complex      umol
%   rocCe        rocuronium       umol/L
%   noreA1,A2    norepinephrine   nmol
%   TPR                           mmHg*mL^-1*min
%   SV                            mL
%   HR                            1/min
%   tdeSV, tdeHR time-dependent   mL, 1/min
%   PaCO2                         mmHg
%   hBis         latent hypnosis  dimensionless, [0,1]
%   ceCO2prop    propofol at the CO2 effect site   ug/mL
%   ceCO2remi    remifentanil at the CO2 effect site ng/mL

names = {'propA1','propA2','propA3','propCe', ...
         'remiA1','remiA2','remiA3','remiCe', ...
         'rocR1','rocR2','sugS1','sugS2','complexRS','rocCe', ...
         'noreA1','noreA2', ...
         'TPR','SV','HR','tdeSV','tdeHR','PaCO2','hBis', ...
         'ceCO2prop','ceCO2remi'};
ix = struct();
for k = 1:numel(names), ix.(names{k}) = k; end
ix.n       = numel(names);
ix.nonNeg  = 1:16;        % amounts and concentrations that cannot go below zero
ix.names   = names;
end
