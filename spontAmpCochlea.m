function [rStarBM, rStarOC, psiStar, stability, stabType] = ...
  spontAmpCochlea(alphaBM, alphaOC, beta1OC, delta1OC, c21, c12, All)
%%
% [rStarBM, rStarOC, psiStar, stability, stabType] = ...
%   spontAmpCochlea(alphaBM, alphaOC, beta1OC, delta1OC, c21, c12, All)
%
% Finds rBM*, rOC*, psi*, stability (1 or 0), and stability
% type (0-4) numerically for the cochlear model with bidirectional
% coupling but no external forcing (see the equations in the code below).
% Set the optional input argument 'All' to 1 (or any nonzero value)
% to get both stable and unstable solutions. (Default for All is 0,
% that is, rStarCochlea outputs only stable fixed points when All is not
% specified.)
%
% stability: 1 = stable, 0 = unstable
% stabType: 4 = stable node, 3 = stable spiral, 2 = unstable node,
% 1 = unstable spiral, 0 = saddle point

%% Equation
% $$\frac{1}{f}\dot{z}_{bm} = z_{bm}\left(\alpha_{bm} +
% 2\pi\textrm{i}\right) + c_{12}z_{oc}$$
% 
% $$\frac{1}{f}\dot{z}_{oc} = z_{oc}\left(\alpha_{oc} +
% 2\pi\textrm{i} + \left(\beta_{1oc} +
% \textrm{i}\delta_{1oc}\right)|z_{oc}|^2\right) + c_{21}z_{bm}$$
%
% where $z_x = r_xe^{\textrm{i}\phi_x}$ and $\psi = \phi_{oc} - \phi_{bm}$.

%% Check input arguments
if c12*c21 == 0
  error('Use spontAmp and rStarDriven11 when c12 or c21 is 0.')
end
if nargin < 7
  All = 0; % default: only stable points
end

abm = alphaBM; aoc = alphaOC; b1oc = beta1OC; d1oc = delta1OC;

%% Get rOC*
roc = sqrt(roots([...
  abm*b1oc*(power(b1oc,2) + power(d1oc,2)),...
  2*power(abm,2)*power(b1oc,2) - power(b1oc,2)*c12*c21 + ...
  abm*aoc*(3*power(b1oc,2) + power(d1oc,2)),...
  (abm + aoc)*b1oc*(power(abm,2) + 3*abm*aoc - 2*c12*c21),...
  power(abm + aoc,2)*(abm*aoc - c12*c21)]));
roc = roc(abs(imag(roc)) < eps('single')); % take only real roots
roc = roc(abs(roc) > 0); % take only positive roots
roc = real(roc);
roc = sort(unique(roc),'descend'); % remove multiple roots

%% Get rBM*
rbm = sqrt(c12/(c21*abm)*(aoc*roc.^2 + b1oc*roc.^4));
ind = find(abs(imag(rbm)) < eps('single'));
rbm = real(rbm(ind));
roc = roc(ind);

%% Get psi*
signPsi = (d1oc*roc.^2./(c12*roc./rbm + c21*rbm./roc) >= 0)*2 - 1;
                                                % sign of psi*oc
psi = signPsi.*acos(-abm*rbm./(c12*roc));
ind = find(abs(imag(psi)) < 10^(-5));
psi = real(psi(ind));
rbm = rbm(ind);
roc = roc(ind);

%% Calculate stability type using Jacobian matrix
stabType = zeros(size(roc));
for n = 1:length(roc)
  J = JacobianSpontAmpCochlea(abm,aoc,b1oc,d1oc,c21,c12,...
    rbm(n),roc(n),psi(n));
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
  psiStar = psi;
else % only stable solutions
  indStab = find(stability);
  rStarBM = rbm(indStab);
  rStarOC = roc(indStab);
  psiStar = psi(indStab);
  stability = stability(indStab);
  stabType = stabType(indStab);
end

%% Calculate Jacobian matrix
function J = JacobianSpontAmpCochlea(abm,aoc,b1oc,d1oc,c21,c12,...
  rbm,roc,psi)
J = NaN(3);
J(1,1) = abm;
J(1,2) = c12*cos(psi);
J(1,3) = -c12*roc*sin(psi);
J(2,1) = c21*cos(psi);
J(2,2) = aoc + 3*b1oc*roc^2;
J(2,3) = -c21*rbm*sin(psi);
J(3,1) = (c12*roc/rbm^2 - c21/roc)*sin(psi);
J(3,2) = 2*d1oc*roc + (c21*rbm/roc^2 - c12/rbm)*sin(psi);
J(3,3) = -(c12*roc/rbm + c21*rbm/roc)*cos(psi);
