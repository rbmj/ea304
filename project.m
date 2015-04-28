close all;

% Note: naca_coords() generates functions for a chord length of 1, so
% all of our points are already in xbar, zbar space (no need to normalize)
[~, mcl] = naca_coords(2412);

% FIXME for now hardcode n, just generate a lift curve for 0-10 deg.
n = 20; % number of lumped vortex panels
alpha = linspace(0, 10, 11);
[cl, Dcp, cm] = arrayfun(@(a) vpm(a*pi/180, n, mcl), alpha);
figure();
plot(alpha, cl);
axis([0 10 0 1.5]);
figure();
plot(alpha, cm);
axis([0 10 -0.5 0.25]);
figure();
cp_distribution = Dcp{1};
plot(cp_distribution(:,1), cp_distribution(:,2));