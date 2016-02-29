%% Define cochlear parameters =================================================
abm  = -412; aoc  =      0;
b1bm = 0;    b1oc =  -40816;
b2bm = 0;    b2oc =      0;
d1bm = 0;    d1oc =      0;
d2bm = 0;    d2oc =      0;
ebm  = 0;    eoc  =      0;

c21     = 197960;
rThresh = 0.1;

%% Make a cochlea network =====================================================
n1 = networkMake(1, 'hopf', abm, b1bm, b2bm, d1bm, d2bm, ebm, ...
                    'log', 30, 10000, 831, ...
                    'display', 50, 'save', 1, 'noScale');
                
n2 = networkMake(2, 'hopf', aoc, b1oc, b2oc, d1oc, d2oc, eoc, ...
                    'log', 30, 10000, 831, ...
                    'display', 50, 'save', 1, 'noScale');
                


%% Make a stimulus ============================================================
F = dB2Pa(9);
freqs = 1850;
index = freqToIndex(n1, freqs);
freqs = n1.f(index);

s = stimulusMake(1, 'fcn', [0 .1], 100000, {'exp'}, freqs, F, 0, ...
                 'ramp', 0.010, 1, 'display', 200);

s.x = midearfilt(s.x, s.fs);

n1 = connectAdd(s, n1, 1, 'noScale');

%% Add connections from bm to oc
bm2oc = diag(c21 * n2.f);

n2    = connectAdd(n1, n2, bm2oc, 'type', '1freq', 'noScale');

%% Add connections from oc to bm
oc2bm = (real(n1.a) ./ (c21 * n2.f) .* (real(n2.a) + real(n2.b1) * (0.5 * rThresh)^2));
oc2bm = diag(oc2bm);

n1    = connectAdd(n2, n1, oc2bm, 'type', '1freq', 'noScale');

ind = round(n1.N/2);
[rStarBM, rStarOC] = spontAmpCochlea(real(n1.a(ind)), real(n2.a(ind)), real(n2.b1),...
    imag(n2.b1), bm2oc(ind, ind), oc2bm(ind, ind));

n1.z0 = rStarBM * ones(n1.N,1);
n1.z  = n1.z0;
n2.z0 = rStarOC * ones(n2.N,1);
n2.z  = n2.z0;

%% Run the network
M = modelMake(@zdot, @cdot, s, n1, n2);

tic;
M = odeRK4fs(M);
clc; toc;

%%
figure(101);
plot(s.t, real(M.n{2}.Z(index,:)));zoom xon
