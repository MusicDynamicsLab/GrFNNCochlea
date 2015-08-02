function [rStarBM, rStarOC, psiStarBM, psiStarOC, stability, stabType] ...
  = rStarCochlea(alphaBM, alphaOC, beta1OC, delta1OC, F, c21, c12, ...
  Omega, All)
%%
% [rStarBM, rStarOC, psiStarBM, psiStarOC] = ...
%   rStarCochlea(alphaBM, alphaOC, beta1OC, delta1OC, F, c21, c12, ...
%   Omega, All)
%
% Finds rBM*, rOC*, psiBM*, psiOC*, stability (1 or 0), and stability
% type (0-4) numerically for the cochlear model with bidirectional
% coupling (see the equations in the code below). When frequency scaled,
% use 2*pi*(f-f0)/f for Omega where f is the natural frequency of the
% oscillator pair and f0 is forcing frequency. When not scaled, use
% 2*pi*(f-f0) instead. Set the optional input argument 'All' to 1 (or
% any nonzero value) to get both stable and unstable fixed points.
% (Default for All is 0, that is, rStarCochlea outputs only stable
% fixed points when All is specified.)
%
% stability: 1 = stable, 0 = unstable
% stabType: 4 = stable node, 3 = stable spiral, 2 = unstable node,
% 1 = unstable spiral, 0 = saddle point

%% Equation
% $$\frac{1}{f}\dot{z}_{bm} = z_{bm}\left(\alpha_{bm} +
% 2\pi\textrm{i}\right) + Fe^{\textrm{i}2\pi f_0t} + c_{12}z_{oc}$$
% 
% $$\frac{1}{f}\dot{z}_{oc} = z_{oc}\left(\alpha_{oc} +
% 2\pi\textrm{i} + \left(\beta_{1oc} +
% \textrm{i}\delta_{1oc}\right)|z_{oc}|^2\right) + c_{21}z_{bm}$$
%
% where $z_x = r_xe^{\textrm{i}\phi_x}$, $\psi_{bm} = \phi_{bm} - 
% 2\pi f_0t$, $\psi_{oc} = \phi_{oc} - \phi_{bm}$.

%% Check input arguments
if nargin < 9
  All = 0; % default: only stable points
end
if F <= 0
  error('F must be positive.')
end

abm = alphaBM; aoc = alphaOC; b1oc = beta1OC; d1oc = delta1OC;
W = Omega;

%% Get rOC*
roc = sqrt(roots([ ...
  power(power(b1oc,2) + power(d1oc,2),2)*(power(abm,2) + ...
  power(W,2)),...
  2*(power(b1oc,2) + power(d1oc,2))*(-(abm*b1oc*c12*c21) + ...
  2*power(abm,2)*(aoc*b1oc + d1oc*W) + W*(c12*c21*d1oc + ...
  2*W*(aoc*b1oc + d1oc*W))),...
  -2*abm*c12*c21*(aoc*(3*power(b1oc,2) + power(d1oc,2)) + ...
  2*b1oc*d1oc*W) + 4*aoc*b1oc*d1oc*W*(c12*c21 + 2*power(W,2)) + ...
  2*power(abm,2)*(power(aoc,2)*(3*power(b1oc,2) + power(d1oc,2)) + ...
  4*aoc*b1oc*d1oc*W + (power(b1oc,2) + 3*power(d1oc,2))*power(W,2)) + ...
  power(b1oc,2)*(power(c12,2)*power(c21,2) + ...
  6*power(aoc,2)*power(W,2) + 2*c12*c21*power(W,2) + 2*power(W,4)) + ...
  power(d1oc,2)*(power(c12,2)*power(c21,2) + ...
  2*power(aoc,2)*power(W,2) + 6*c12*c21*power(W,2) + 6*power(W,4)),...
  -(power(b1oc,2)*power(c21,2)*power(F,2)) - ...
  power(c21,2)*power(d1oc,2)*power(F,2) + ...
  2*power(c12,2)*power(c21,2)*d1oc*W + 4*power(aoc,3)*b1oc*power(W,2) + ...
  6*c12*c21*d1oc*power(W,3) + 4*d1oc*power(W,5) + ...
  4*power(abm,2)*(aoc*b1oc + d1oc*W)*(power(aoc,2) + power(W,2)) + ...
  2*power(aoc,2)*d1oc*W*(c12*c21 + 2*power(W,2)) - ...
  2*abm*c12*c21*(3*power(aoc,2)*b1oc + 2*aoc*d1oc*W + ...
  b1oc*power(W,2)) + 2*aoc*b1oc*(power(c12,2)*power(c21,2) + ...
  2*c12*c21*power(W,2) + 2*power(W,4)),...
  -2*aoc*b1oc*power(c21,2)*power(F,2) + power(aoc,4)*power(W,2) - ...
  2*abm*aoc*c12*c21*(power(aoc,2) + power(W,2)) + ...
  power(abm,2)*power(power(aoc,2) + power(W,2),2) + ...
  power(aoc,2)*(power(c12,2)*power(c21,2) + 2*c12*c21*power(W,2) + ...
  2*power(W,4)) + W*(2*c12*c21*power(W,3) + power(W,5) + ...
  power(c21,2)*(-2*d1oc*power(F,2) + power(c12,2)*W)),...
  -(power(c21,2)*power(F,2)*(power(aoc,2) + power(W,2)))]));

roc = roc(abs(imag(roc)) < eps('single')); % take only real roots
roc = roc(abs(roc) > 0); % take only positive roots
roc = real(roc);
roc = sort(unique(roc),'descend'); % remove multiple roots

%% Get rBM*
rbm = sqrt((aoc*roc+b1oc*roc.^3).^2 + (W*roc + d1oc*roc.^3).^2)/c21;
ind = find(abs(imag(rbm)) < eps('single'));
rbm = real(rbm(ind));
roc = roc(ind);

%% Get psiOC*
signPsiOC = ((W + d1oc*roc.^2)*c21 >= 0)*2 - 1; % sign of psi*oc
psioc = signPsiOC.*acos(-(aoc*roc+b1oc*roc.^3)./(c21*rbm));
ind = find(abs(imag(psioc)) < 10^(-5));
psioc = real(psioc(ind));
rbm = rbm(ind);
roc = roc(ind);

%% Get psiBM*
signPsiBM = (W + c12*roc.*sin(psioc)./rbm >= 0)*2 - 1; % sign of psi*oc
psibm = signPsiBM.*acos(-(abm*rbm+c12*roc.*cos(psioc))/F);
ind = find(abs(imag(psibm)) < 10^(-5));
psibm = real(psibm(ind));
rbm = rbm(ind);
roc = roc(ind);
psioc = psioc(ind);

%% Calculate stability type using Jacobian matrix
stabType = zeros(size(roc));
for n = 1:length(roc)
  J = JacobianCochlea(abm,aoc,b1oc,d1oc,F,c21,c12,rbm(n),roc(n),...
    psibm(n),psioc(n));
  ev = eig(J);
  if isreal(ev) && all(ev < 0)
    stabType(n) = 4; % stable node
  elseif all(real(ev) < 0)
    stabType(n) = 3; % stable spiral
  elseif isreal(ev) && all(ev > 0)
    stabType(n) = 2; % unstable node
  elseif all(real(ev) > 0)
    stabType(n) = 1; % unstable spiral
  end % saddle pt otherwise
end
stability = (stabType >= 3); % 1 = stable, 0 = unstable

%% Prepare output
if All % both stable and unstable solutions
  rStarBM = rbm;
  rStarOC = roc;
  psiStarBM = psibm;
  psiStarOC = psioc;
else % only stable solutions
  indStab = find(stability);
  rStarBM = rbm(indStab);
  rStarOC = roc(indStab);
  psiStarBM = psibm(indStab);
  psiStarOC = psioc(indStab);
  stability = stability(indStab);
  stabType = stabType(indStab);
end

%% Calculate Jacobian matrix
function J = JacobianCochlea(abm,aoc,b1oc,d1oc,F,c21,c12,rbm,roc,...
  psibm,psioc)
J = NaN(4);
J(1,1) = abm;
J(1,2) = -F*sin(psibm);
J(1,3) = c12*cos(psioc);
J(1,4) = -c12*roc*sin(psioc);
J(2,1) = (F*sin(psibm)-c12*roc*sin(psioc))/rbm^2;
J(2,2) = -F*cos(psibm)/rbm;
J(2,3) = c12*sin(psioc)/rbm;
J(2,4) = c12*roc*cos(psioc)/rbm;
J(3,1) = c21*cos(psioc);
J(3,2) = 0;
J(3,3) = aoc + 3*b1oc*roc^2;
J(3,4) = -c21*rbm*sin(psioc);
% J(4,1) = -c21*sin(psioc)/roc;
% J(4,2) = 0;
% J(4,3) = 2*d1oc*roc + c21*rbm*sin(psioc)/roc^2;
% J(4,4) = -c21*rbm*cos(psioc)/roc;
J(4,1) = -c21*sin(psioc)/roc - (F*sin(psibm)-c12*roc*sin(psioc))/rbm^2;
J(4,2) = F*cos(psibm)/rbm;
J(4,3) = 2*d1oc*roc + (c21*rbm/roc^2 - c12/rbm)*sin(psioc);
J(4,4) = -(c21*rbm/roc + c12*roc/rbm)*cos(psioc);
