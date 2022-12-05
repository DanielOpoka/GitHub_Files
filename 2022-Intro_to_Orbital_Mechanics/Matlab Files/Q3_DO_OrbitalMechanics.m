
%% ------------------------------------------------------------------------
%                          ORBITAL MECHANICS - PART III
%  ------------------------------------------------------------------------

% Analysis of the time-frequency domain of signals through FFT and STFT.

clear; clc; 

close all;
n_fontsize = 14;
nb_fig = 0; % For an ever increasing number of plots
new_dir = "images/Q3/";
line_thickness = 2.0;

format long;

%% 3a) % b)

R_E = 6378;                                             % Radius of the Earth [km]
miu = 3.986e5;                                          % Gravitational parameter of the Earth [km^3 / s^2]

alt_ISS = 404;                                          % ISS altitude [km]
nu_0_cha = 180;                                         % True Anomaly of Progress (chaser) [º]
delta_Omega = 100;                                      % Angular Separation between ISS and Progress (chaser) [º]

% ISS Orbital Parameters
a_ISS = (R_E + alt_ISS) * 1e3;                          % Semi-major axis [m]
tau_ISS = 2 * pi * sqrt((a_ISS / 1e3)^3 / miu);         % Orbital period [s] 

[r_ISS, v_ISS] = keplerian2ijk(a_ISS, 0, 0, 0, 0, 0, 'truelon', delta_Omega);

nb_days = 5;
time = 0:10:(nb_days * (24 * 60^2));
X_ISS_0 = [r_ISS; v_ISS] ./ 1e3;

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); 
[time_ISS_out, X_ISS_out] = ode45(@twobody, time, X_ISS_0, options);

X_ISS_out_pos = X_ISS_out(:, 1:3);
X_ISS_out_vel = X_ISS_out(:, 4:6);

% For Earth sphere
[X_sph, Y_sph, Z_sph] = sphere;
X_sph_E = X_sph * R_E;
Y_sph_E = Y_sph * R_E;
Z_sph_E = Z_sph * R_E;

N_rev_cha = 12;                                         % Number of Revolutions for Progress (chaser)[#] 
% Semi-major axis for Progress (chaser) [m]
a_cha = (1 - (deg2rad(delta_Omega) / (2 * pi * N_rev_cha)))^(2/3) * a_ISS;
ecc_cha = (a_ISS / a_cha) - 1;                          % Eccentricity for Progress (chaser) []
tau_cha = 2 * pi * sqrt((a_cha / 1e3)^3 / miu);         % Orbital period [s] 

[r_cha, v_cha] = keplerian2ijk(a_cha, ecc_cha, 0, 0, 0, 180, 'lonper', 180);
X_cha_0 = [r_cha; v_cha] ./ 1e3;

[time_cha_out, X_cha_out] = ode45(@twobody, time, X_cha_0, options);

X_cha_out_pos = X_cha_out(:, 1:3);
X_cha_out_vel = X_cha_out(:, 4:6);

dist_ISS_cha = vecnorm(X_cha_out_pos' - X_ISS_out_pos');

%% 3 c)

r_ISS = a_ISS * 1e-3;                                    % Semi-major axis of ISS [km]
v_ISS = sqrt(miu / r_ISS);                              % Velocity of the ISS [km / s]

N_rev_cha_range = 2:30;                                 % Number of Revolutions for Progress (chaser)[#] 
% Semi-major axis for Progress (chaser) [m]
a_cha_range = (1 - (deg2rad(delta_Omega) ./ (2 .* pi .* N_rev_cha_range))).^(2/3) .* a_ISS;
ecc_cha_range = (a_ISS ./ a_cha_range) - 1;             % Eccentricity for Progress (chaser) []
% Orbital period [s] 
tau_cha_range = 2 .* pi .* sqrt((a_cha_range ./ 1e3).^3 ./ miu);  

delta_V_range = zeros(1, length(N_rev_cha_range));

for ind = 1:length(N_rev_cha_range)
    % Radius of apogee of Progress (chaser) [km]
    r_a_cha = a_cha_range(ind) * 1e-3 * (1 + ecc_cha_range(ind));
%     v_cha_ind = sqrt((2 * miu) * ((1 / r_a_cha) - (1 / (r_ISS + r_a_cha))))
    v_cha_ind = sqrt((miu) * ((2 / r_a_cha) - (1 / (a_cha_range(ind) * 1e-3))))
    delta_V_range(ind) = v_ISS - v_cha_ind;
end

%% d)

r_p_cha_range = 60:210;
r_p_cha_range = r_p_cha_range + R_E;

a_cha1_range = (r_ISS + r_p_cha_range) / 2;
tau_cha1_range = 2 .* pi .* sqrt((a_cha1_range).^3 ./ miu);  

delta_V1_range = zeros(1, length(r_p_cha_range));

for ind = 1:length(r_p_cha_range)
%     v_cha_ind = sqrt((2 * miu) * ((1 / r_ISS) - (1 / (r_ISS + r_p_cha_range(ind)))));
    v_cha_ind = sqrt((miu) * ((2 / r_ISS) - (1 / (a_cha1_range(ind)))));
    delta_V1_range(ind) = v_cha_ind - v_ISS;
end


%% Plots

close all

nb_fig = nb_fig + 1;
figure(nb_fig)
set(gcf,'position',[50 , 50, 800, 800])

% axis equal;
hold on;
grid on;
% sphere_surf = surf(X_sph_E, Y_sph_E, Z_sph_E);
% set(sphere_surf,'FaceColor', [1 0 0]);
% colormap jet;

X_out = X_ISS_out;
plot3(X_out(:, 1), X_out(:, 2), X_out(:, 3), 'b', 'LineWidth', line_thickness);
plot3(X_out(1, 1), X_out(1, 2), X_out(1, 3), 'b*', 'LineWidth', line_thickness);

X_out = X_cha_out;
plot3(X_out(:, 1), X_out(:, 2), X_out(:, 3), 'r', 'LineWidth', line_thickness);
plot3(X_out(1, 1), X_out(1, 2), X_out(1, 3), 'r*', 'LineWidth', line_thickness);

X_dist = [X_ISS_out(1, 1:3); zeros(1, 3); X_cha_out(1, 1:3)];
plot3(X_dist(:, 1), X_dist(:, 2), X_dist(:, 3), 'k--', 'LineWidth', line_thickness);

hold off;

view(0, 90)

% legend('Earth', 'Orbit')

set(gca,'fontsize',n_fontsize);
% Labelling
title("\textbf{Spacecraft's Orbit}", 'interpreter', 'latex')
xlabel('x-axis [km]', 'interpreter', 'latex');
ylabel('y-axis [km]', 'interpreter', 'latex');
zlabel('z-axis [km]', 'interpreter', 'latex');

legend('$r_{ISS}$', '$r_{ISS,0}$', '$r_{cha}$', '$r_{cha,0}$', 'interpreter', 'latex');

figure_title = new_dir + "orb_ISS.png";
saveas(gcf, figure_title)

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot(time ./ (tau_cha), dist_ISS_cha, 'b', 'LineWidth', line_thickness)
xline(N_rev_cha, 'k', 'LineWidth', line_thickness);
hold off;

set(gca,'fontsize',n_fontsize);

% Labelling
title("\textbf{Distance between}", '\textbf{ISS and Progress (chaser)}', 'interpreter', 'latex')
xlabel('Time / $\tau$ [$\#_{orb}$]', 'interpreter', 'latex');
ylabel('Distance [km]', 'interpreter', 'latex');

legend('', '$nb_{orb}$ = 12', 'interpreter', 'latex', 'location', 'best')

figure_title = new_dir + "orb_dist_ISS_cha.png";
saveas(gcf, figure_title)

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot(N_rev_cha_range, delta_V_range, 'LineWidth', line_thickness)
xline(N_rev_cha, 'k', 'LineWidth', line_thickness);
hold off;

set(gca,'fontsize',n_fontsize);

% Labelling
title("\textbf{Velocity Change}", 'interpreter', 'latex')
xlabel("Number of Revolutions [$\#_{rev}$]", 'interpreter', 'latex');
ylabel('$\Delta$V [km/s]', 'interpreter', 'latex');

legend('', '$nb_{orb}$ = 12', 'interpreter', 'latex', 'location', 'best')

figure_title = new_dir + "orb_dV_cha.png";
saveas(gcf, figure_title)

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot(N_rev_cha_range, tau_cha_range / (60^2), 'LineWidth', line_thickness)
xline(N_rev_cha, 'k', 'LineWidth', line_thickness);
hold off;

set(gca,'fontsize',n_fontsize);

% Labelling
title("\textbf{Orbit Period}", 'interpreter', 'latex')
xlabel("Number of Revolutions [$\#_{rev}$]", 'interpreter', 'latex');
ylabel('$\tau_{cha}$ [hr]', 'interpreter', 'latex');

legend('', '$nb_{orb}$ = 12', 'interpreter', 'latex', 'location', 'best')

figure_title = new_dir + "orb_tau_cha.png";
saveas(gcf, figure_title)

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot(r_p_cha_range - R_E, delta_V1_range, 'LineWidth', line_thickness)
hold off;

set(gca,'fontsize',n_fontsize);

% Labelling
title("\textbf{Velocity Change for De-Orbiting}", 'interpreter', 'latex')
xlabel("Altitude [km]", 'interpreter', 'latex');
ylabel('$\Delta$V [km/s]', 'interpreter', 'latex');

figure_title = new_dir + "orb_dV1_cha.png";
saveas(gcf, figure_title)

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot(r_p_cha_range - R_E, tau_cha1_range ./ 60, 'LineWidth', line_thickness)
hold off;

set(gca,'fontsize',n_fontsize);

% Labelling
title("\textbf{Orbit Period for De-Orbiting}", 'interpreter', 'latex')
xlabel("Altitude [km]", 'interpreter', 'latex');
ylabel('$\tau$ [min]', 'interpreter', 'latex');

figure_title = new_dir + "orb_tau1_cha.png";
saveas(gcf, figure_title)

%% Functions

function orb_acc = twobody(t, X)
    miu = 3.986e5;                                      % Gravitational parameter of the Earth [km^3 / s^2]

    orb_acc = zeros(6, 1);
    
    orb_acc(1) = X(4);
    orb_acc(2) = X(5);
    orb_acc(3) = X(6);
 
    mag = sqrt(X(1)^2 + X(2)^2 + X(3)^2);               % Magnitude of the radial elements [km]
    orb_acc(4) = -(miu / mag^3) * X(1);
    orb_acc(5) = -(miu / mag^3) * X(2);
    orb_acc(6) = -(miu / mag^3) * X(3);
end
