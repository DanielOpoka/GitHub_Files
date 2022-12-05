
%% ------------------------------------------------------------------------
%                          ORBITAL MECHANICS - PART I
%  ------------------------------------------------------------------------

% Analysis of the time-frequency domain of signals through FFT and STFT.

clear; clc; 

close all;
n_fontsize = 14;
nb_fig = 0; % For an ever increasing number of plots
new_dir = "images/Q1/";
line_thickness = 2.0;

format long;

%% 1b) Kepler
M = deg2rad(21);                                        % Mean anomaly [rad]
e = 0.25;                                               % Eccentricity []

tol = 1e-6;                                             % Tolerance []
ite_max = 100;                                          % Max iterations [#]
[E] = Kepler(M, e, tol, ite_max)                       % Eccentric anomaly [rad]
%%
tol = 1e-12;                                            % Tolerance []
ite_max = 100;                                          % Max iterations [#]
[E] = Kepler(M, e, tol, ite_max);
tht = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));   % True anomaly [rad]
a = 24000;                                              % Semi-major axis [km]
p = a * (1 - e^2);                                      % Semi-lactus rectum [km]
r = p / (1 + e * cos(tht));                             % Radial distance [km]
%%
M = deg2rad(180);                                       % Mean anomaly [rad]
[E] = Kepler(M, e, tol, ite_max);
tht = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));   % True anomaly [rad]
a = 24000;                                              % Semi-major axis [km]
p = a * (1 - e^2);                                      % Semi-lactus rectum [km]
r = p / (1 + e * cos(tht));                             % Radial distance [km]

% r becomes the radius of the apogee as tht is 180°, it points towards the
% other side of the orbit.

%% 1c) & 1d) - MEO

% Orbital Parameters
R_E = 6378;                                             % Radius of the Earth [km]
a = 24e3;                                               % Semi-major axis [km]
e = 0.72;                                               % Eccentricity []
p = a * (1 - e^2);                                      % Semi-lactus rectum [km]
miu = 3.986e5;                                          % Gravitational parameter of the Earth [km^3 / s^2]
tau = 2 * pi * sqrt(a^3 / miu);                         % Orbital period [s] 
n = (2 * pi) / tau;                                     % Mean Motion [rad / s]
r_p = a * (1 - e);                                      % Radius of perigee [km]
r_a = a * (1 + e);                                      % Radius of apogee [km]

% Calculation Parameters
tol = 1e-6;                                             % Tolerance []
ite_max = 100;                                          % Max iterations [#]

% Eclipse Parameters
alpha = R_E^2 * e^2 + p^2;
beta = 2 * R_E^2 * e;
gamma = R_E^2 - p^2;

f_ecs = @(theta)(alpha .* cos(theta).^2 + beta .* cos(theta) + gamma);
tht_ref = 0;

% For eclipse in perigee
tht_range = 0:0.001:((1 / 2) * pi);
tht_ecs = f_ecs(tht_range);
diff = abs(tht_ecs - tht_ref);
[min_val, min_ind] = min(diff);
tht_11 = tht_range(min_ind);

tht_range = (-(1 / 2) * pi):0.001:0;
tht_ecs = f_ecs(tht_range);
diff = abs(tht_ecs - tht_ref);
[min_val, min_ind] = min(diff);
tht_12 = tht_range(min_ind);

% For eclipse in apogee
tht_range = ((1 / 2) * pi):0.001:(1 * pi);
tht_ecs = f_ecs(tht_range);
diff = abs(tht_ecs - tht_ref);
[min_val, min_ind] = min(diff);
tht_21 = tht_range(min_ind);

tht_range = (-pi):0.001:(-(1 / 2) * pi);
tht_ecs = f_ecs(tht_range);
diff = abs(tht_ecs - tht_ref);
[min_val, min_ind] = min(diff);
tht_22 = tht_range(min_ind);

% % Altitude for 1 orbit
% Simulation Parameters
t0 = 0;                                                 % Time at perigee [s]
dt = 15;                                                % Time step [s]
tend = tau;                                             % End time (1 orbital period) [s]
time = t0:dt:tend;                                      % Time [s]

alt = zeros(1, length(time));                           % Altitude [km]

orb_x = zeros(1, length(time));                         % x-axis of orbit [km]
orb_y = zeros(1, length(time));                         % y-axis of orbit [km]

tht_1 = zeros(1, length(time));                         % Theta of eclipse on perigee [rad]
tht_2 = zeros(1, length(time));                         % Theta of eclipse on apogee [rad]

ite_k = 0;

for t_k = t0:dt:tend
    ite_k = ite_k + 1;
    
    M_k = n * (t_k - t0);
    [E_k] = Kepler(M_k, e, tol, ite_max);
    tht_k = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E_k / 2));
    r_k = p / (1 + e * cos(tht_k));
    
    orb_x(ite_k) = r_k * cos(tht_k);
    orb_y(ite_k) = r_k * sin(tht_k);
    
    alt(ite_k) = r_k - R_E;
end

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot((time ./ tau) .* 100, alt, 'LineWidth', line_thickness)
hold off;

set(gca,'fontsize',n_fontsize);
% Labelling
title('\textbf{Altitude of the Satellite}', 'interpreter', 'latex')
xlabel('Time / $\tau$ [\%]', 'interpreter', 'latex');
ylabel('Altitude [km]', 'interpreter', 'latex');

figure_title = new_dir + "alt_1orb.png";
saveas(gcf, figure_title)

% % Altitude for 1 day
% Simulation Parameters
t0 = 0;                                                 % Time at perigee [s]
dt = 15;                                                % Time step [s]
tend = 24 * 60 * 60;                                    % End time (1 day) [s]
time = t0:dt:tend;                                      % Time [s]

alt = zeros(1, length(time));                           % Altitude [km]

orb_x = zeros(1, length(time));                         % x-axis of orbit [km]
orb_y = zeros(1, length(time));                         % y-axis of orbit [km]

orb_x_ecs_p = zeros(1, length(time));                   % x-axis of orbit during perigee eclipse [km]
orb_y_ecs_p = zeros(1, length(time));                   % y-axis of orbit during perigee eclipse [km]
time_ecs_p = 0;                                         % Time spent in perigee eclipse [s]

orb_x_ecs_a = zeros(1, length(time));                   % x-axis of orbit during apogee eclipse [km]
orb_y_ecs_a = zeros(1, length(time));                   % y-axis of orbit during apogee eclipse [km]
time_ecs_a = 0;                                         % Time spent in apogee eclipse [s]

tht_k = zeros(1, length(time));                         % Theta of eclipse on perigee [rad]

ite_k = 0;

for t_k = t0:dt:tend
    ite_k = ite_k + 1;
    
    M_k = n * (t_k - t0);
    [E_k] = Kepler(M_k, e, tol, ite_max);
    tht_k(ite_k) = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E_k / 2));
    r_k = p / (1 + e * cos(tht_k(ite_k)));
    
    orb_x(ite_k) = r_k * cos(tht_k(ite_k));
    orb_y(ite_k) = r_k * sin(tht_k(ite_k));
    
    if (tht_k(ite_k) <= tht_11) && (tht_k(ite_k) >= tht_12)
        orb_x_ecs_p(ite_k) = orb_x(ite_k);
        orb_y_ecs_p(ite_k) = orb_y(ite_k);
        time_ecs_p = time_ecs_p + dt;
    elseif (tht_k(ite_k) >= tht_21) 
        orb_x_ecs_a(ite_k) = orb_x(ite_k);
        orb_y_ecs_a(ite_k) = orb_y(ite_k);
        time_ecs_a = time_ecs_a + dt;
    elseif(tht_k(ite_k) <= tht_22)
        orb_x_ecs_a(ite_k) = orb_x(ite_k);
        orb_y_ecs_a(ite_k) = orb_y(ite_k);
        time_ecs_a = time_ecs_a + dt;
    end
    
    alt(ite_k) = r_k - R_E;
end

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
plot(time ./ (60^2), alt, 'LineWidth', line_thickness)
hold off;

set(gca,'fontsize',n_fontsize);
% Labelling
title('\textbf{Altitude of the Satellite}', 'interpreter', 'latex')
xlabel('Time [hr]', 'interpreter', 'latex');
ylabel('Altitude [km]', 'interpreter', 'latex');

figure_title = new_dir + "alt_1day.png";
saveas(gcf, figure_title)

% % 1 orbit (x-axis and y-axis)
nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
axis equal
orbit = plot(orb_x, orb_y, 'LineWidth', 3.0);
orbit_ecs_p = plot(orb_x_ecs_p, orb_y_ecs_p, 'm--', 'LineWidth', 1.0);
orbit_ecs_a = plot(orb_x_ecs_a, orb_y_ecs_a, 'g--', 'LineWidth', 1.0);
earth = viscircles([0, 0], R_E, 'LineWidth', line_thickness);
hold off;

legend([orbit, orbit_ecs_p, orbit_ecs_a, earth], ...
    {'orbit', 'eclipse at perigee', 'eclipse at apogee', 'Earth'}, 'location' ,'best')

set(gca,'fontsize',n_fontsize);
% Labelling
title('\textbf{Orbit of the Satellite}', 'interpreter', 'latex')
xlabel('x-axis [km]', 'interpreter', 'latex');
ylabel('y-axis [km]', 'interpreter', 'latex');

figure_title = new_dir + "orb_1orb.png";
saveas(gcf, figure_title)

% The solution for the eclipse
tht_range = (-pi):0.001:(pi);
tht_ecs = f_ecs(tht_range);

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
yline(0);
plot(tht_range, tht_ecs, 'r', 'LineWidth', line_thickness)
plot(tht_11, 0, 'g*', tht_12, 0, 'g*', tht_21, 0, 'm*', tht_22, 0, 'm*', 'LineWidth', 3.0)
hold off;

xlim([-pi, pi])

legend('', '', '\theta_{11}', '\theta_{12}', '\theta_{21}', '\theta_{22}', 'location', 'best')

set(gca,'fontsize',n_fontsize);
% Labelling
title('\textbf{Eclipse Equation}', 'interpreter', 'latex')
xlabel('$\theta [^{\circ}]$', 'interpreter', 'latex');
ylabel('$\theta_{ec} [^{\circ}]$', 'interpreter', 'latex');

figure_title = new_dir + "ecs_eq.png";
saveas(gcf, figure_title)

%% 2


%% Functions

function [Ek] = Kepler(M, e, tol, ite_max)
Ek = M;
f_E = @(E)(E - e * sin(E) - M);
df_E = @(E)(1 - e * cos(E));

tol_E = tol + 1;   
ite = 0;

while (tol_E > tol && ite < ite_max)
    ite = ite + 1;
    Ek1 = Ek;
    delta_E =  (f_E(Ek)) / (df_E(Ek));
    ite
    Ek = Ek - delta_E
    tol_E = abs(Ek - Ek1);
end
end

