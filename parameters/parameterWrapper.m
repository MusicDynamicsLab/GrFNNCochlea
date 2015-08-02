%% Fitted cochlear params 
% These curves look weird when scaled by f
% What can we make of them? No frequency scaling?

clear all

n1 = networkMake(1, 'hopf', -0.10,       0, 0, 0, 0, 0.0000, 'log', 30, 20000, 831, 'channel', 1, 'display', 200, 'save', 10);
n1.df = 1;
n2 = networkMake(2, 'hopf',  0.00, -200000, 0, 0, 0, 0.0025, 'log', 30, 20000, 831, 'display', 200, 'save', 10);
A = 1;
bm2oc = connectMake(n1, n2, 'one', 1, 1, 0, 1);
n2    = connectAdd(n1, n2, bm2oc, 'weight', A);

parameters;
figure(11); semilogx(f_stim, alpha_bm.*f_stim, 'o', n1.f, real(n1.a));

b = regress(log(-alpha_bm.*f_stim), [log(f_stim), ones(size(f_stim))]);
logalpha_bm = b(2) + b(1)*log(n1.f); 
newalpha_bm = -exp(logalpha_bm);

hold on; semilogx(n1.f, newalpha_bm, 'k'); hold off;
set(gca, 'FontSize', 20)
grid
title('\alpha_b_m', 'FontSize', 20)

figure(12); semilogx(f_stim, beta_oc.*f_stim, 'o', n2.f, real(n2.b1));

b = regress(log(-beta_oc.*f_stim), [log(f_stim), ones(size(f_stim))]);
logbeta_oc = b(2) + b(1)*log(n1.f); 
newbeta_oc = -exp(logbeta_oc);

hold on; semilogx(n2.f, newbeta_oc, 'k'); hold off;
set(gca, 'FontSize', 20)
grid
title('\beta_o_c', 'FontSize', 20)


figure(13); semilogx(f_stim, A_conn.*f_stim, 'o', n1.f, n2.con{1}.w, 'b');
n2.con{1}.w = 2000; % interp1(f_stim, A_conn, n1.f).*n1.f; % -- looks flat
hold on; semilogx(n1.f, n2.con{1}.w, 'k'); hold off;

b = regress(log(A_conn.*f_stim), [log(f_stim), ones(size(f_stim))]);
logA_conn = b(2) + b(1)*log(n1.f); 
newA_conn = exp(logA_conn);
hold on; semilogx(n1.f, newA_conn, 'm'); hold off;

set(gca, 'FontSize', 20)
grid
title('A', 'FontSize', 20)

