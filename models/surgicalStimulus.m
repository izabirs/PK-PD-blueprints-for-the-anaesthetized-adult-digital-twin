function s = surgicalStimulus(tMin, profile)
%SURGICALSTIMULUS  Normalised surgical stimulus s(t) in [0,1].
%
%   s(t) is the seventh element of the input vector. It drives the
%   nociceptive layer through the stimulus intensity of the response
%   surface, and the haemodynamic layer through a multiplicative
%   perturbation of the production rates.
%
%   PROFILE is an n-by-3 matrix of [tStart tEnd amplitude] in minutes.
%   The default is a generic case: intubation, incision, a surgical
%   plateau, deep manipulation and closure.

if nargin < 2 || isempty(profile)
    profile = [ 5   7  0.85     % laryngoscopy and intubation
               15  17  0.90     % skin incision
               17  45  0.55     % surgical plateau
               45  55  0.95     % deep manipulation
               55  75  0.50     % surgical plateau
               75  85  0.25];   % closure
end

s = 0;
for i = 1:size(profile,1)
    if tMin >= profile(i,1) && tMin < profile(i,2)
        ramp = min((tMin - profile(i,1))/0.5, 1);
        s = max(s, profile(i,3)*ramp);
    end
end
end
