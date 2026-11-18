function p = arousalProbability(bis)
%AROUSALPROBABILITY  Probability of clinical arousal at neuromuscular
%   reversal, gated by the prevailing hypnotic depth.
%
%   Two-point logistic calibrated to the two arms of
%   Kim J, Kim JA, Jang JN, Yang M, Ahn HJ, Choi J, Jo S. The arousal
%   effect of sugammadex reversal of neuromuscular blockade differs with
%   anesthetic depth in propofol-remifentanil anesthesia.
%   Sci Rep 2023;13:20776:  12 of 33 (36.4%) at BIS 45-55, 0 of 33 at
%   BIS 25-35, after 4 mg/kg sugammadex from deep block.
%
%   Two points do not identify a logistic curve. This function carries the
%   correct qualitative behaviour rather than an estimated relationship;
%   prospective calibration is an item of the research agenda.

p = 0.728 ./ (1 + exp(-(bis - 50)/3));
end
