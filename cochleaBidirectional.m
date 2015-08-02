%% Make a cochlea network =====================================================
n1 = networkMake(1, 'hopf', -500, -1000000, 0, 0, 0, 0.0, ...
                    'log', 30, 10000, 831, 'channel', 1, ...
                    'display', 200, 'save', 1);
n2 = networkMake(2, 'hopf',   -0,   -40000, 0, 0, 0, 0.0, ...
                    'log', 30, 10000, 831, ...
                    'display', 200, 'save', 1);

n1.z0 = 0*n1.z;
n2.z0 = .025 + 0*n2.z;
n1.z  = 0*n1.z;
n2.z  = .025 + 0*n2.z;

%% Make a stimulus ============================================================
F = dB2Pa(20);
freqs = 1000;
index = freqToIndex(n1, freqs);
freqs = n1.f(index);

s = stimulusMake('fcn', [0 .10], 100000, {'exp'}, freqs, F, 0, ...
                 'ramp', 0.010, 1, 'display', 200);

stim_rms = rms(s.x);
ratio = (F/sqrt(2))/stim_rms;
s.x = ratio*s.x;

s.x=midearfilt(s.x,s.fs);

%%
abm  = -412; aoc  =      0;
b1bm = 0;    b1oc =  -40e3;
b2bm = 0;    b2oc = -1.6e6;
d1bm = 0;    d1oc =      0;
d2bm = 0;    d2oc =      0;
ebm  = 0;    eoc  =    .04;

c21 = 200e3;
c12 =     0;
r_thresh = 0.1;

%% Add connections from bm to oc

bm2oc = connectMake(n1, n2, 'one', 1, 1, 0, 1);
oc2oc = diag(ones(length(bm2oc)-1,1), 1) + diag(ones(length(bm2oc)-1,1), -1);

n2    = connectAdd(n1, n2, bm2oc, 'weight', c21, 'type', '1freq');
n1    = connectAdd(n2, n1, bm2oc, 'weight', c12, 'type', '1freq');

%% Fitted cochlear params 
n1.w  = 1;
n1.a  = abm + i*2*pi.*n1.f;
n2.a  = aoc + i*2*pi.*n1.f;
n2.b1 = b1oc;
n2.b2 = b2oc;
n1.con{1}.w = (real(n1.a)./n2.con{1}.w(end)) .* (real(n2.a) +  real(n2.b1) * (0.50*r_thresh).^2);

%% Run the network
M = modelMake(@zdot, @cdot, s, n1, n2);

tic;
M = odeRK4fs(M, s);
clc; toc;

%%
% figure(101);
% plot(s.t, real(M.n{1}.Z(index,:)));
