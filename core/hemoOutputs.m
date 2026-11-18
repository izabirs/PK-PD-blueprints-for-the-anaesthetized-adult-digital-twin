function o = hemoOutputs(x, L4, nore, cNore, cNoreBase)
%HEMOOUTPUTS  MAP, cardiac output and the apparent SV and HR from the
%   turnover states (Su 2022/2023).
ix  = twinIndex();
tpr = x(ix.TPR);
sv  = x(ix.SV) + x(ix.tdeSV);
hr  = x(ix.HR) + x(ix.tdeHR);
baseHRapp = L4.BaseHR*(1 + L4.LTDE_HR);
svApp = sv * (1 - L4.HR_SV*log(max(hr,1e-3)/baseHRapp));
mapv  = tpr * svApp * hr;
if nargin > 3 && ~isempty(nore) && cNore > cNoreBase
    rise = hill(cNore,     nore.pd.EC50, nore.pd.gamma) - ...
           hill(cNoreBase, nore.pd.EC50, nore.pd.gamma);
    mapv = mapv * (1 + nore.pd.fracMax * max(rise, 0));
end
o.TPR = tpr;  o.SV = svApp;  o.HR = hr;  o.MAP = mapv;  o.CO = svApp*hr/1000;
end
