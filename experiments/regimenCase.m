function u = regimenCase(t, o, weight, sugammadexTime)
%REGIMENCASE  Open-loop clinical regimen for the illustrative case.
%
%   Propofol 2 mg/kg induction then 5 mg/kg/h; remifentanil
%   0.12 ug/kg/min, raised to 0.22 during deep manipulation; rocuronium
%   0.6 mg/kg with a top-up; norepinephrine titrated to the MAP deficit;
%   sugammadex 4 mg/kg at reversal.

if nargin < 3 || isempty(weight), weight = 75; end
if nargin < 4, sugammadexTime = 85; end
stopPropofol = 78;  stopRemi = 80;

u = struct('propofol',0,'remifentanil',0,'rocuronium',0, ...
           'sugammadex',0,'norepinephrine',0);

if t < 0.5
    u.propofol = 2.0*weight/0.5;                 % 2 mg/kg over 30 s
elseif t < stopPropofol
    u.propofol = 5.0*weight/60;                  % 5 mg/kg/h
end
if t < stopRemi
    if t >= 44 && t < 56, rate = 0.22; else, rate = 0.12; end
    u.remifentanil = rate*weight;                % ug/min
end
if t >= 2 && t < 2.5,  u.rocuronium = 0.6*weight/0.5;  end
if t >= 40 && t < 40.5, u.rocuronium = 0.15*weight/0.5; end
if ~isempty(sugammadexTime) && t >= sugammadexTime && t < sugammadexTime + 0.5
    u.sugammadex = 4.0*weight/0.5;               % 4 mg/kg over 30 s
end
if t > 3 && o.MAP < 65
    u.norepinephrine = min(max(0.02*(65 - o.MAP),0),0.15)*weight;
end
end
