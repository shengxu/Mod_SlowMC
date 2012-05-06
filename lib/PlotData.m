clc;
clear all;

xs = load('U235_300K_18');

xs2 = load('../../SlowMC/lib/U235_300K_18');
% xs2 = [xs(:,1) 11.5860*ones(length(xs), 1)];

figure;
loglog(xs(:,1), xs(:,2), '-k');
hold on;
loglog(xs2(:,1), xs2(:,2), '-r');
legend('NJOY', 'Broadened');
hold off;