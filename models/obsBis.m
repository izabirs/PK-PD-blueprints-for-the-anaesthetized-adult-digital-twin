function bis = obsBis(h, pd, stimulus, tofRatio, confounders, noiseSd)
%OBSBIS  Observation model linking the latent hypnotic state to BIS.
%
%   BIS = E0*(1 - h) + bEMG*s(t)*TOF(t) + sum_j b_j*1_j(t) + v(t)
%
%   The electromyographic term requires both a stimulus and recovered
%   neuromuscular transmission, so it is gated by the product of the
%   stimulus and the train-of-four ratio and vanishes under deep block.
%   The remaining bias terms are the documented confounder catalogue of
%   Dahaba AA. Anesth Analg 2005;101:765-773: ketamine, nitrous oxide,
%   dexmedetomidine, hypothermia, cerebral hypoperfusion and electrical
%   interference. Their weights are not established quantitatively; the
%   fields default to zero and exist so that the structure is explicit.
%
%   h           latent hypnotic depth in [0,1]
%   pd          struct from pdBisSurface
%   stimulus    surgical stimulus in [0,1]
%   tofRatio    train-of-four ratio in [0,1]
%   confounders struct, optional, fields: bEMG, ketamine, n2o,
%               dexmedetomidine, hypothermia, hypoperfusion, interference
%   noiseSd     measurement noise standard deviation, BIS units

if nargin < 5 || isempty(confounders), confounders = struct; end
if nargin < 6, noiseSd = 0; end
f = @(name,default) ternary(isfield(confounders,name), ...
                            @() confounders.(name), @() default);

bEMG = f('bEMG', 0);
bias = bEMG * stimulus * tofRatio;
for name = {'ketamine','n2o','dexmedetomidine','hypothermia', ...
            'hypoperfusion','interference'}
    bias = bias + f(name{1}, 0);
end

bis = pd.E0*(1 - h) + bias;
if noiseSd > 0, bis = bis + noiseSd*randn; end
bis = min(max(bis, 0), 100);
end

function out = ternary(cond, fTrue, fFalse)
if cond, out = fTrue(); else, out = fFalse(); end
end
