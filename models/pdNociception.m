function m = pdNociception()
%PDNOCICEPTION  Hierarchical surface for the probability of tolerating a
%   noxious stimulus of a given intensity.
%
%       post = iota * (1 - Ce_remi^gr / ((C50r*iota)^gr + Ce_remi^gr))
%       Ptol = Ce_prop^gp / ((C50p*post)^gp + Ce_prop^gp)
%
%   Source: Bouillon TW, et al. Anesthesiology 2004;100:1353-1372
%   (tolerance of laryngoscopy).
%
%   iota is the stimulus intensity. iota0 is the pre-stimulus intensity
%   estimated for laryngoscopy; kappa scales surgical stimuli relative to
%   laryngoscopy and is NOT identified by any published study. It is the
%   least constrained quantity in the specification.

m.C50p   = 8.04;    % ug/mL
m.C50r   = 1.07;    % ng/mL
m.gammaP = 5.1;
m.gammaR = 0.97;
m.iota0  = 1.05;    % laryngoscopy
m.kappa  = 1.3;     % ASSUMPTION, see README
m.source = 'Bouillon 2004, Anesthesiology 100:1353-1372';
end
