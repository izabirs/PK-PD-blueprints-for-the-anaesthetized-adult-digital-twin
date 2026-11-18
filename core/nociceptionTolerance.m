function p = nociceptionTolerance(cePro, ceRem, stimulus, L2)
%NOCICEPTIONTOLERANCE  Probability of tolerating a stimulus of intensity
%   iota = iota0*kappa*s(t), hierarchical surface of Bouillon 2004.
iota = L2.iota0 * L2.kappa * max(stimulus, 1e-6);
post = iota * (1 - hill(ceRem, L2.C50r*iota, L2.gammaR));
p    = hill(cePro, L2.C50p*post, L2.gammaP);
end
