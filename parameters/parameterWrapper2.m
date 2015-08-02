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
figure(11); semilogx(f_stim, alpha_bm, '-o');
% hold on; loglog(n1.f, newalpha_bm, 'k'); hold off;
set(gca, 'FontSize', 20)
grid
title('\alpha_b_m', 'FontSize', 20)

figure(12); semilogx(f_stim, beta_oc, '-o');
% n2.beta1 = interp1(f_stim, beta_oc, n2.f);
% n2.b1 = n2.beta1.*n2.f;
set(gca, 'FontSize', 20)
grid
title('\beta_o_c', 'FontSize', 20)


figure(13); semilogx(f_stim, A_conn, '-o');
n2.con{1}.w = 2000; % interp1(f_stim, A_conn, n1.f).*n1.f; % -- looks flat

% hold on; loglog(n1.f, n2.con{1}.w, 'k'); hold off;
set(gca, 'FontSize', 20)
grid
title('A', 'FontSize', 20)

