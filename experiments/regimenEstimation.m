function u = regimenEstimation(t, o, weight)
%REGIMENESTIMATION  Induction, two maintenance levels and emergence.
%
%   The step change at 25 min keeps the hypnotic signal away from
%   saturation, so that the fitting window contains information about the
%   individual parameters.

if nargin < 3 || isempty(weight), weight = 70; end
u = struct('propofol',0,'remifentanil',0,'rocuronium',0, ...
           'sugammadex',0,'norepinephrine',0);
if t < 0.5
    u.propofol = 1.8*weight/0.5;
elseif t < 25
    u.propofol = 4.0*weight/60;
elseif t < 50
    u.propofol = 6.0*weight/60;
end
if t < 52, u.remifentanil = 0.15*weight; end
if t >= 2 && t < 2.5, u.rocuronium = 0.6*weight/0.5; end
end
