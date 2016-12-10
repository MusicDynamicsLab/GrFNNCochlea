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
                    'display', 0, 'save', 1, 'znaught', 0, 'noScale');
                
n2 = networkMake(2, 'hopf', aoc, b1oc, b2oc, d1oc, d2oc, eoc, ...
                    'log', 30, 10000, 831, ...
                    'display', 0, 'save', 1, 'znaught', 0, 'noScale');

%% Make a stimulus ============================================================
F = dB2Pa(10);
freqs = 1850;
index = freqToIndex(n1, freqs);
freqs = n1.f(index);

s = stimulusMake(1, 'fcn', [0 .10], 100000, {'exp'}, freqs, F, 0, ...
                 'ramp', 0.010, 1, 'display', 0);

s.x = midearfilt(s.x, s.fs);

n1 = connectAdd(s, n1, 1, 'noScale');

%% Add connections from bm to oc
bm2oc = diag(c21 * n2.f);

n2    = connectAdd(n1, n2, bm2oc, 'type', '1freq', 'noScale');

%% Run the network
M = modelMake(@zdot, @cdot, s, n1, n2);

tic


%% Uncomment this block to use ode45
ind = 1;
zColumnLength = 0;
for i=1:length(M.n)
    zColumnLength = zColumnLength + M.n{i}.N;
end
z = NaN(zColumnLength,1);
for i=1:length(M.n)
    z(ind:ind + M.n{i}.N - 1) = M.n{i}.z;
    ind = ind + M.n{i}.N;
end
options=odeset('RelTol',1e-2);
M.odefun = @ode45;
[T,Z] = M.odefun(@zdotAdaptive,s.t,z,options,M);

ind = 1;
for i = 1:length(M.n)
    M.n{i}.t=T;
    M.n{i}.Z=Z(:,ind:ind+M.n{i}.N - 1).';
    ind = ind + M.n{i}.N;
end


%% Uncomment this line to use odeRK4fs
% M = odeRK4fs(M);

%%
toc

%%
figure
plot(s.t, real(M.n{2}.Z(index,:)));zoom xon
