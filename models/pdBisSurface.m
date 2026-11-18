function m = pdBisSurface()
%PDBISSURFACE  Propofol-remifentanil response surface for the latent
%   hypnotic state, read out as the bispectral index.
%
%   Minto-type surface on normalised effect-site concentrations:
%       U     = Ce_prop/C50p + Ce_remi/C50r
%       h     = U^gamma / (1 + U^gamma)        (latent hypnotic depth)
%       BIS   = E0 * (1 - h)                   (see obsBis for the
%                                               measurement model)
%
%   Source: Bouillon TW, Bruhn J, Radulescu L, Andresen C, Shafer TJ,
%   Cohane C, Shafer SL. Pharmacodynamic interaction between propofol and
%   remifentanil regarding hypnosis, tolerance of laryngoscopy, bispectral
%   index, and electroencephalographic approximate entropy.
%   Anesthesiology 2004;100:1353-1372.
%
%   Interaction structure is a configuration choice; a reduced Greco or a
%   three-drug MAZA surface (Liou 2018) can be substituted at this
%   interface without changing the rest of the twin.

m.C50p  = 4.47;    % ug/mL
m.C50r  = 19.3;    % ng/mL
m.gamma = 1.43;
m.E0    = 97.4;
m.tau   = 15;      % s, cortical and processing delay as a first-order lag
m.cv    = struct('C50p',0.182,'C50r',0.888,'gamma',0.304);
m.source = 'Bouillon 2004, Anesthesiology 100:1353-1372';
end
