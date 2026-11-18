function twin = applyTheta(twin, theta)
%APPLYTHETA  Apply log-scale individual parameters to a twin.
%   theta = [log Ce50 scale, log clearance scale, log ke0 scale]
twin.L1.C50p     = twin.L1.C50p * exp(theta(1));
twin.cfg.clScale = exp(theta(2));
twin.prop.ke0    = twin.prop.ke0 * exp(theta(3));
twin.x0          = twinInitialState(twin);
twin.cNoreBase   = twin.x0(twin.ix.noreA1)/twin.nore.V1;
o        = twinOutputs(0, twin.x0, twin, 0, false);
twin.CO0 = o.CO;
end
