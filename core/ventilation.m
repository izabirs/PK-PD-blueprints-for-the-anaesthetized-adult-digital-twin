function ve = ventilation(paco2, ceCO2prop, ceCO2remi, twin)
%VENTILATION  Minute ventilation from arterial CO2 and the CO2-effect-site
%   concentrations, with the Nieuwenhuijs 2003 interaction.
L3 = twin.L3;
ur = max(ceCO2remi,0)/L3.C50r;
up = max(ceCO2prop,0)/L3.C50p;
u  = ur + up;
if twin.cfg.respSynergy
    if u > 0, q = ur/u; else, q = 0.5; end
    iq  = 1 + (L3.Imax - 1)*4*q*(1 - q);       % Qmax = 0.5, symmetric
    dep = min(u^L3.gamma * 0.5 * iq, 1);
else
    er  = min((ur)^L3.gamma * 0.5, 1);
    ep  = min((up)^L3.gamma * 0.5, 1);
    dep = 1 - (1 - er)*(1 - ep);
end
drive = 1 + L3.F*(paco2 - L3.PaCO2base)/L3.PaCO2base;
ve    = L3.VEbase * max(drive,0) * (1 - dep);
ve    = max(ve, 0);
end
