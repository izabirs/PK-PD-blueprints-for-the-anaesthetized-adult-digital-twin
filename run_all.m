%RUN_ALL  Reproduce every number and figure of the supplementary material.
%
%   Verifies each layer against its source publication, then runs the
%   three demonstrations. Total runtime is roughly 10 minutes, dominated
%   by the individualisation experiment.

setupPaths();
fprintf('=== Verification against the source publications ===\n');
verifyAgainstLiterature(true);

fprintf('\n=== Demonstration 1 ===\n');
res1 = exp1_coupledCase(true);

fprintf('\n=== Demonstration 2 ===\n');
res2 = exp2_individualisation(25, true);

fprintf('\n=== Demonstration 3 ===\n');
res3 = exp3_varianceIdentifiability(300, true);

save(fullfile(fileparts(mfilename('fullpath')),'figures','results.mat'), ...
     'res1','res2','res3');
fprintf('\nResults saved to figures/results.mat\n');
