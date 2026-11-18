function out = twinSimulate(twin, durationMin, dtMin, controlFcn, opts)
%TWINSIMULATE  Integrate the coupled twin with fixed-step RK4.
%
%   out = twinSimulate(twin, duration, dt, controlFcn) runs the twin for
%   DURATION minutes with step DT. CONTROLFCN is a function handle
%   u = controlFcn(t, o) returning the input struct given the current time
%   and the current outputs.
%
%   opts.spontaneousFcn  @(t) logical, spontaneous ventilation
%   opts.bloodLossFcn    @(t) fractional stroke-volume loss per minute
%   opts.noiseSd         measurement noise on BIS and MAP, default 0
%   opts.x0              initial state, default twin.x0

if nargin < 5, opts = struct; end
g = @(f,v) ternary(isfield(opts,f), @() opts.(f), @() v);
spontFcn = g('spontaneousFcn', @(t) false);
bleedFcn = g('bloodLossFcn',   @(t) 0);
noiseSd  = g('noiseSd', 0);
x        = g('x0', twin.x0);

ix = twin.ix;
n  = round(durationMin/dtMin);
fields = {'t','BIS','MAP','HR','SV','TPR','CO','TOF','VE','PaCO2','Ptol', ...
          'CePropofol','CeRemifentanil','CeRocuronium','CNorepinephrine', ...
          'stimulus','hypnoticDepth'};
for f = fields, out.(f{1}) = zeros(n+1,1); end
out.uPropofol = zeros(n+1,1);  out.uRemifentanil = zeros(n+1,1);
out.uNorepinephrine = zeros(n+1,1);
out.X = zeros(n+1, ix.n);

for k = 0:n
    t  = k*dtMin;
    sp = spontFcn(t);
    o  = twinOutputs(t, x, twin, [], sp, noiseSd);
    u  = controlFcn(t, o);
    if ~isfield(u,'stimulus'),    u.stimulus    = o.stimulus; end
    if ~isfield(u,'spontaneous'), u.spontaneous = sp;         end
    u.bloodLoss = bleedFcn(t);

    for f = fields, out.(f{1})(k+1) = o.(f{1}); end
    out.uPropofol(k+1)       = getfielddef(u,'propofol',0);
    out.uRemifentanil(k+1)   = getfielddef(u,'remifentanil',0);
    out.uNorepinephrine(k+1) = getfielddef(u,'norepinephrine',0);
    out.X(k+1,:) = x(:)';
    if k == n, break; end

    f1 = twinDerivatives(t,            x,               u, twin);
    f2 = twinDerivatives(t + dtMin/2,  x + dtMin/2*f1,  u, twin);
    f3 = twinDerivatives(t + dtMin/2,  x + dtMin/2*f2,  u, twin);
    f4 = twinDerivatives(t + dtMin,    x + dtMin*f3,    u, twin);
    x  = x + dtMin/6*(f1 + 2*f2 + 2*f3 + f4);
    x(ix.nonNeg) = max(x(ix.nonNeg), 0);
end
end

function v = getfielddef(s,f,d)
if isfield(s,f), v = s.(f); else, v = d; end
end
function out = ternary(c,a,b)
if c, out = a(); else, out = b(); end
end
