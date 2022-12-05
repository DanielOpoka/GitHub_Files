
%% ------------------------------------------------------------------------
%                          ORBITAL MECHANICS - PART II
%  ------------------------------------------------------------------------

% Analysis of the time-frequency domain of signals through FFT and STFT.

clear; clc; 

close all;
n_fontsize = 14;
nb_fig = 0; % For an ever increasing number of plots
new_dir = "images/Q2/";
line_thickness = 2.0;

format long;

%% 2 a), b) & c)

% Orbital Parameters
R_E = 6378;                                             % Radius of the Earth [km]
miu = 3.986e5;                                          % Gravitational parameter of the Earth [km^3 / s^2]

time = 0:10:(24 * 60^2);                                % Time [s]
% Initial State Vector for ODE [km, km, km, 
%                               km/s, km/s, km/s]
X_0 = [7115.48; 3391.696; 3492.221; 
    -3.762; 4.063; 4.184];

X_0_pos = [X_0(1), X_0(2), X_0(3)] .* 1e3;
X_0_vel = [X_0(4), X_0(5), X_0(6)] .* 1e3;
mag_X_0_pos = sqrt(X_0(1).^2 + X_0(2).^2 + X_0(3).^2);
mag_X_0_vel = sqrt(X_0(4).^2 + X_0(5).^2 + X_0(6).^2);

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); 
[time_out, X_out] = ode45(@twobody, time, X_0, options);

X_out_pos = X_out(:, 1:3);
X_out_vel = X_out(:, 4:6);

mag_X_out_pos = sqrt(X_out(:, 1).^2 + X_out(:, 2).^2 + X_out(:, 3).^2);
mag_X_out_vel = sqrt(X_out(:, 4).^2 + X_out(:, 5).^2 + X_out(:, 6).^2);

% For Earth sphere
[X_sph, Y_sph, Z_sph] = sphere;
X_sph_E = X_sph * R_E;
Y_sph_E = Y_sph * R_E;
Z_sph_E = Z_sph * R_E;

%% 2 d)

% Kinetic Energy [km^2/s^2]
eps_kin = mag_X_out_vel .^2 / 2;
eps_pot = -1 * (miu ./ mag_X_out_pos);
eps_tot = eps_kin + eps_pot;

omega_dot = (mag_X_out_vel ./ mag_X_out_pos);
h_out = mag_X_out_pos.^2 .* omega_dot;

%% 2 e)

time_unk = 0:10:(24 * 60^2);                            % Time [s]
% Initial State Vector for ODE [3x1 position vector [km],
%                               3x1 velocity vector [km/s]]
X_unk_0 = [0; 0; 8550; ...
    0; -7.0; 0];

X_0_unk_pos = [X_unk_0(1), X_unk_0(2), X_unk_0(3)] .* 1e3;
X_0_unk_vel = [X_unk_0(4), X_unk_0(5), X_unk_0(6)] .* 1e3;
mag_X_0_unk_pos = sqrt(X_unk_0(1).^2 + X_unk_0(2).^2 + X_unk_0(3).^2);
mag_X_0_unk_vel = sqrt(X_unk_0(4).^2 + X_unk_0(5).^2 + X_unk_0(6).^2);

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); 
[time_out_unk, X_unk_out] = ode45(@twobody, time_unk, X_unk_0, options);

X_unk_pos = X_unk_out(:, 1:3);
X_unk_vel = X_unk_out(:, 4:6);

mag_X_unk_pos = sqrt(X_unk_out(:, 1).^2 + X_unk_out(:, 2).^2 + X_unk_out(:, 3).^2);
mag_X_unk_vel = sqrt(X_unk_out(:, 4).^2 + X_unk_out(:, 5).^2 + X_unk_out(:, 6).^2);

% From Cartesian vector to Keplerian elements
[a_unk, ecc_unk, incl_unk, RAAN_unk, argp_unk, nu_unk, ...
    truelon_unk, arglat_unk, lonper_unk] = ijk2keplerian(X_0_unk_pos, X_0_unk_vel);

tau_unk = 2 * pi * sqrt((a_unk / 1e3)^3 / miu);          % Orbital period [s] 

%% Plots
close all;
 
nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot(time_out ./ 60^2, mag_X_out_pos, 'b', 'LineWidth', line_thickness)
hold off;

set(gca,'fontsize',n_fontsize);
% Labelling
title("\textbf{Spacecraft's Magnitude}", '\textbf{of the Position Vector}', 'interpreter', 'latex')
xlabel('Time [hr]', 'interpreter', 'latex');
ylabel('Radius [km]', 'interpreter', 'latex');

figure_title = new_dir + "orb_mag_pos.png";
saveas(gcf, figure_title)

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot(time_out ./ 60^2, mag_X_out_vel, 'r', 'LineWidth', line_thickness)
hold off;

set(gca,'fontsize',n_fontsize);
% Labelling
title("\textbf{Spacecraft's Magnitude}", '\textbf{of the Velocity Vector}', 'interpreter', 'latex')
xlabel('Time [hr]', 'interpreter', 'latex');
ylabel('Velocity [km/s]', 'interpreter', 'latex');

figure_title = new_dir + "orb_mag_vel.png";
saveas(gcf, figure_title)

nb_fig = nb_fig + 1;
figure(nb_fig)
set(gcf,'position',[50 , 50, 800, 800])

axis equal;
hold on;
grid on;
sphere_surf = surf(X_sph_E, Y_sph_E, Z_sph_E);
set(sphere_surf,'FaceColor', [1 0 0])
plot3(X_out(:, 1), X_out(:, 2), X_out(:, 3), 'LineWidth', line_thickness);
colormap jet;
hold off;

view(45, 45)

legend('Earth', 'Orbit')

set(gca,'fontsize',n_fontsize);
% Labelling
title("\textbf{Spacecraft's Orbit}", 'interpreter', 'latex')
xlabel('x-axis [km]', 'interpreter', 'latex');
ylabel('y-axis [km]', 'interpreter', 'latex');
zlabel('z-axis [km]', 'interpreter', 'latex');

figure_title = new_dir + "orb_3D.png";
saveas(gcf, figure_title)

% % Unknown Orbit
nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot(time_out ./ 60^2, eps_kin, 'r', 'LineWidth', line_thickness)
plot(time_out ./ 60^2, eps_pot, 'g', 'LineWidth', line_thickness)
plot(time_out ./ 60^2, eps_tot, 'b', 'LineWidth', line_thickness)
hold off;

set(gca,'fontsize',n_fontsize);
% Labelling
title("\textbf{Spacecraft's Orbit}", 'Specific Energy', 'interpreter', 'latex')
xlabel('Time [hr]', 'interpreter', 'latex');
ylabel('Energy [$km^2$/$s^2$]', 'interpreter', 'latex');

legend('kinetic', 'potential', 'total', 'location' ,'east')

figure_title = new_dir + "orb_eps.png";
saveas(gcf, figure_title)

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot(time_out ./ 60^2, h_out, 'g', 'LineWidth', line_thickness)
hold off;

set(gca,'fontsize',n_fontsize);
% Labelling
title("\textbf{Spacecraft's Orbit}", '\textbf{Specific Angular Momentum}', 'interpreter', 'latex')
xlabel('Time [hr]', 'interpreter', 'latex');
ylabel('Angular Momentum [$km^2$/s]', 'interpreter', 'latex');

figure_title = new_dir + "orb_ang_mom.png";
saveas(gcf, figure_title)

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot(time_out ./ 60^2, mag_X_unk_pos, 'b', 'LineWidth', line_thickness)
hold off;

set(gca,'fontsize',n_fontsize);
% Labelling
title("\textbf{Unknown Spacecraft's Magnitude}", '\textbf{of the Position Vector}', 'interpreter', 'latex')
xlabel('Time [hr]', 'interpreter', 'latex');
ylabel('Radius [km]', 'interpreter', 'latex');

figure_title = new_dir + "orb_unk_mag_pos.png";
saveas(gcf, figure_title)

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot(time_out ./ 60^2, mag_X_unk_vel, 'r', 'LineWidth', line_thickness)
hold off;

set(gca,'fontsize',n_fontsize);
% Labelling
title("\textbf{Unknown Spacecraft's Magnitude}", '\textbf{of the Position Vector}', 'interpreter', 'latex')
xlabel('Time [hr]', 'interpreter', 'latex');
ylabel('Velocity [km/s]', 'interpreter', 'latex');

figure_title = new_dir + "orb_unk_mag_vel.png";
saveas(gcf, figure_title)


nb_fig = nb_fig + 1;
figure(nb_fig)
set(gcf,'position',[50 , 50, 800, 800])

axis equal;
hold on;
grid on;
sphere_surf = surf(X_sph_E, Y_sph_E, Z_sph_E);
set(sphere_surf,'FaceColor', [1 0 0])
plot3(X_unk_out(:, 1), X_unk_out(:, 2), X_unk_out(:, 3), 'LineWidth', line_thickness);
colormap jet;
hold off;

view(90, 0)

legend('Earth', 'Orbit')

set(gca,'fontsize',n_fontsize);
% Labelling
title("\textbf{Spacecraft's Orbit}", 'interpreter', 'latex')
xlabel('x-axis [km]', 'interpreter', 'latex');
ylabel('y-axis [km]', 'interpreter', 'latex');
zlabel('z-axis [km]', 'interpreter', 'latex');

figure_title = new_dir + "orb_unk_3D.png";
saveas(gcf, figure_title)

%% Functions

function orb_acc = twobody(t, X)
% Gravitational parameter of the Earth [km^3 / s^2]
miu = 3.986e5;                   

orb_acc = zeros(6, 1);

orb_acc(1) = X(4);
orb_acc(2) = X(5);
orb_acc(3) = X(6);

% Magnitude of the radial elements [km]
mag = sqrt(X(1)^2 + X(2)^2 + X(3)^2);               
orb_acc(4) = -(miu / mag^3) * X(1);
orb_acc(5) = -(miu / mag^3) * X(2);
orb_acc(6) = -(miu / mag^3) * X(3);
end
