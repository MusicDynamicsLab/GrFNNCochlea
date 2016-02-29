%% Define cochlear parameters =================================================
abm  = -412; aoc  =      0;
b1bm = 0;    b1oc =  -40816;
b2bm = 0;    b2oc =      0;
d1bm = 0;    d1oc =      0;
d2bm = 0;    d2oc =      0;
ebm  = 0;    eoc  =      0;

c21 = 197960;

%% Make a cochlea network =====================================================
n1 = networkMake(1, 'hopf', abm, b1bm, b2bm, d1bm, d2bm, ebm, ...
                    'log', 30, 10000, 831, ...
                    'display', 50, 'save', 1, 'znaught', 0, 'noScale');
                
n2 = networkMake(2, 'hopf', aoc, b1oc, b2oc, d1oc, d2oc, eoc, ...
                    'log', 30, 10000, 831, ...
                    'display', 50, 'save', 1, 'znaught', 0, 'noScale');

%% Make a stimulus ============================================================
F = dB2Pa(10);
freqs = 1850;
index = freqToIndex(n1, freqs);
freqs = n1.f(index);

s = stimulusMake(1, 'fcn', [0 .10], 100000, {'exp'}, freqs, F, 0, ...
                 'ramp', 0.010, 1, 'display', 200);

s.x = midearfilt(s.x, s.fs);

n1 = connectAdd(s, n1, 1, 'noScale');

%% Add connections from bm to oc
bm2oc = diag(c21 * n2.f);

n2    = connectAdd(n1, n2, bm2oc, 'type', '1freq', 'noScale');

%% Run the network
M = modelMake(@zdot, @cdot, s, n1, n2);

tic;
M = odeRK4fs(M);
clc; toc;

%%
figure(101);
plot(s.t, real(M.n{2}.Z(index,:)));zoom xon
