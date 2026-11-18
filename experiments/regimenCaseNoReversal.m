function u = regimenCaseNoReversal(t, o)
%REGIMENCASENOREVERSAL  The case regimen without neuromuscular reversal,
%   used where the reversal event would confound the analysis.
u = regimenCase(t, o, 75, []);
end
