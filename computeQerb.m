function [Qerb,ERB] = computeQerb(f,TC,CF)
% function [Qerb,ERB] = computeQerb(f,TC,?CF?)
% Computes Qerb=CF/ERB for the tuning curve TC(f).
% The TC is assumed to be in dB. The CF is found
% from the minimum of the TC if not specified.
% C.A.Shera

[minTC,idx] = min(TC);
if (nargin<3)
  CF = f(idx);                          
end

b = f./CF;                              % normalize the frequency
TC = minTC - TC;                        % invert and normalize the tuning curve
t = 10.^(TC/20);                        % t is 1 at its peak
Qerb = 1./trapz(b,t.^2);
if (nargout>1)
  ERB = CF/Qerb;
end