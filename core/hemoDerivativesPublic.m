function d = hemoDerivativesPublic(x, L4, cpProp, cpRemi, s, bleed)
%HEMODERIVATIVESPUBLIC  Haemodynamic turnover derivatives, exposed so that
%   the layer can be exercised on its own (see verifyAgainstLiterature).
ix = twinIndex();
effPropTPR = (L4.EmaxpropTPR + L4.INTtpr*hill(cpRemi, L4.EC50remiTPR, 1)) * ...
             hill(cpProp, L4.EC50propTPR, L4.gammaProp);
effPropSV  = L4.EmaxpropSV * hill(cpProp, L4.EC50propSV, 1);
effRemiTPR = L4.EmaxremiTPR * hill(cpRemi, L4.EC50remiTPR, 1);
effRemiSV  = (L4.SLremiSV + L4.INTsv*hill(cpProp, L4.EC50propSV, 1)) * cpRemi;
effRemiHR  = (L4.SLremiHR + L4.INThr*hill(cpProp, L4.EC50intHR, 1)) * cpRemi;
clip = @(v) min(max(v,-0.999),0.999);
effPropTPR = clip(effPropTPR); effPropSV = clip(effPropSV);
effRemiTPR = clip(effRemiTPR); effRemiSV = clip(effRemiSV);
effRemiHR  = clip(effRemiHR);

hemo  = hemoOutputs(x, L4);
baseMAP = L4.BaseTPR*L4.BaseSV*(1+L4.LTDE_SV)*L4.BaseHR*(1+L4.LTDE_HR);
rmap  = max(hemo.MAP/baseMAP, 1e-6);
fb    = rmap^(-L4.FB);
fstim = 1 + L4.thetaStim*s;
kinTPR = L4.kout*L4.BaseTPR; kinSV = L4.kout*L4.BaseSV; kinHR = L4.kout*L4.BaseHR;

d = zeros(5,1);
d(1) = kinTPR*fstim*fb*(1 + effPropTPR) - L4.kout*x(ix.TPR)*(1 - effRemiTPR);
d(2) = kinSV *fb*(1 + effPropSV)        - L4.kout*x(ix.SV) *(1 - effRemiSV) - bleed*L4.BaseSV;
d(3) = kinHR *fstim*fb                  - L4.kout*x(ix.HR) *(1 - effRemiHR);
if cpProp > 0
    d(4) = -L4.kTDE*x(ix.tdeSV);
    d(5) = -L4.kTDE*x(ix.tdeHR);
end
end
