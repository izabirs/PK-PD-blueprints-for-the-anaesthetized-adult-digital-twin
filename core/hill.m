function y = hill(c, c50, gamma)
%HILL  Sigmoid Emax fraction c^gamma/(c50^gamma + c^gamma), zero at c<=0.
c = max(c, 0);
if c == 0, y = 0; else, y = c^gamma / (c50^gamma + c^gamma); end
end
