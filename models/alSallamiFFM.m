function ffm = alSallamiFFM(isMale, weight, age, bmi)
%ALSALLAMIFFM  Al-Sallami fat-free mass, the size descriptor used by the
%   Eleveld propofol and remifentanil models.
%
%   Source: Al-Sallami HS, et al. Prediction of fat-free mass in children.
%   Clin Pharmacokinet 2015;54:1169-1178, as applied in
%   Eleveld 2018 (Br J Anaesth 120:942-959) and Eleveld 2017
%   (Anesthesiology 126:1005-1018).

if isMale
    ffm = (0.88 + (1-0.88)/(1 + (age/13.4)^-12.7)) * (9270*weight)/(6680 + 216*bmi);
else
    ffm = (1.11 + (1-1.11)/(1 + (age/7.1)^-1.1))   * (9270*weight)/(8780 + 244*bmi);
end
end
