function ce = perrinRemifentanilCe(dNOL, rounded)
%PERRINREMIFENTANILCE  Personalised remifentanil effect-site concentration
%   for skin incision, from the nociception level index response to a
%   standardised tetanic stimulus.
%
%   Source: Perrin L, Bisdorff M, Saxena S, et al. Predicting personalised
%   remifentanil effect site concentration for surgical incision using the
%   nociception level index. Eur J Anaesthesiol 2022;39:918-927.
%
%   Fitted:  Ce = 0.069*dNOL + 3.1   (95% CI 0.069 [0.03-0.11]; 3.1 [2.3-4.4])
%   The authors round this to Ce = 0.1*dNOL + 3.1 for clinical use.
%
%   dNOL is the maximum index value minus the value before the stimulus,
%   measured during a 30 s, 100 Hz, 50 mA tetanic stimulation of the distal
%   ulnar nerve at a remifentanil Ce of 4 ng/mL with propofol titrated to a
%   bispectral index of 40-60.
%
%   This is a titration heuristic, not a forward pharmacodynamic model: no
%   model of the form NOL = f(Ce_prop, Ce_remi, stimulus) has been
%   published, which is an item of the research agenda.

if nargin < 2, rounded = false; end
if rounded, ce = 0.1*dNOL + 3.1; else, ce = 0.069*dNOL + 3.1; end
end
