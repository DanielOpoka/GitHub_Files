
%% ------------------------------------------------------------------------
%                             CORRELATION RANGE
%  ------------------------------------------------------------------------

% Seeing how correlation changes from one emitter position to another. 

clear all; close all; clc; format long;
addpath("functions");
addpath("data/H_fs");

new_dir_0 = "images/EstEmtPos"; 
mkdir(new_dir_0);

new_dir_1 = new_dir_0 + "/CorrRange"; 
mkdir(new_dir_1);

nb_fig = 0;
n_fontsize = 12;
       

%% ------------------------------------------------------------------------
%                          SIMULATION (VALENTIN'S CODE)
%  ------------------------------------------------------------------------

% Definition of simulation of a Reverberant Chamber.

%% PARAMETERS

% Emitter parameters
signal_type = "wavelet";

nb_signals = 1;

d_impulse = 8e-9;                                                      % Duration of the impulse [s]
f0 = 2.4e9;                                                            % Carrier frequency of the impulse [Hz]

d_before = d_impulse * 2;                                                  % Delay before impulse [s]
d_after = 500e-9;                                                          % Delay after impulse [s]
A_impulse = 1;                                                             % Amplitude of the impulse [V/m]

% Cavity parameters
Lx=1.9; Ly=1.9;                                                            % Cavity dimensions [m, m]
sigma=(6.5)*1e2;                                                           % Conductivity of the walls [S/m]
coef_losses_ok=1;                                                          % If=1: losses, if=0: PEC condition 
coef_regular_chao=2;                                                       % If=1: regular cavity, if=2: chaotic cavity 

% WARNING: lossless chaotic cavity not coded (possible: lossless or with losses reglar cavity, or with losses chaotic cavity)
[time, e_t, param_simu] = function_parameters_simu(Lx,Ly,d_impulse,A_impulse,f0,d_before,d_after,coef_regular_chao,coef_losses_ok,sigma, signal_type);

% Wave parameters
c0 = 3e8;                                                                  % Speed of light [m/s]
f_impulse = 1 / d_impulse;                                                 % Frequency of impulse [Hz]
lambda_min = c0 / (f0 + 2 * f_impulse);                                    % Minimum wavelength of signal [m]
fact_delta=10;                                                             % Discretization []
delta=lambda_min/fact_delta;                                               % Discretized space step (dx, dy) [m]
dt=delta/c0/sqrt(2);                                                       % Discretized time step [s]
fs = 1 / dt;                                                               % Sampling Frequency [Hz]

% Receiver and emitter locations:
pos_rec = [1.0; 1.4];                                                      % Location of the receiver [m]
pos_emt = [1.0; 0.6];                                                % Location of the emitter [m]

% Display of the impulse respone
Display_rep_imp = 0;                                                       % If=1: display, if=0: no display
[r_t, Ez] = function_impulse_response(param_simu, pos_emt, pos_rec, Display_rep_imp, time, e_t);

% Frequency Range
df = 1 / time(end);
freq = [0:length(time)-1] * df;

new_dir_test = new_dir_1 + "/test-" + "emt_x=" + num2str(round(pos_emt(1), 3)) + ...
    "m_y=" + num2str(round(pos_emt(2), 3)) + "m";
mkdir(new_dir_test);


%% PLOTTING OF EMITTER(S) AND RECEIVER(S)

openfig("data/CavitySpace/cavity_space_empty.fig");
figure(1);

pos_emt_x = pos_emt(1, :);
pos_emt_y = pos_emt(2, :);

pos_rec_x = pos_rec(1, :);
pos_rec_y = pos_rec(2, :);

marker_size = 20;
line_width = 2;

hold on;
plot(pos_emt_x, pos_emt_y, 'g*', 'MarkerSize', marker_size, 'LineWidth', line_width)
plot(pos_rec_x, pos_rec_y, 'rx', 'MarkerSize', marker_size, 'LineWidth', line_width)
hold off;

set(gca,'fontsize',n_fontsize);

title({"Cavity with " + num2str(size(pos_emt, 2)) + ...
    " Emitter(s) and " + num2str(size(pos_rec, 2)) + " Receiver(s)"});
xlabel("x-axis [m]");
ylabel("y-axis [m]");

figure_title = new_dir_test + "/cavity-" + num2str(size(pos_emt, 2)) + "_emt(s)-" + ...
    num2str(size(pos_rec, 2)) + "_rec(s).png";
saveas(gcf, figure_title)


%% ------------------------------------------------------------------------
%                           CREATING CLOSE EMITTERS
%  ------------------------------------------------------------------------

% Creating the second emitter where the correlation will be done.

%% PARAMETERS

% % Fourier Transform: e(t) -> E(f) (UNKNOWN)
% For Normalization:     
N_e_t = nnz(e_t);                                                          % Number of elements != 0 in e_t
E_f = fft(e_t)';                                                           % FFT of e(t)
ph_E_f = unwrap(angle(E_f));                                               % Phase angle of e(t)
E_f_norm = abs(E_f) ./ N_e_t;                                              % Normalized E(f)
% PSD of e(t)
[P_e_t, f_e_t] = periodogram(e_t, hann(length(e_t), 'periodic'), [], fs);
% 99% Occupied Bandwidth of e(t)
[bw_e_t, flo_e_t, fhi_e_t, pow_e_t] = obw(P_e_t, f_e_t);

% % Fourier Transform: r(t) -> R(f) (KNOWN)
% For Normalization:
abs_r_t_cond = abs(r_t) > 0;                                               % Conditioner if abs(r_t) is greater than 0
N_r_t = sum(abs_r_t_cond == 1);                                            % Number of elements != 0 in r_t
R_f = fft(r_t);                                                            % FFT of r(t)
ph_R_f = unwrap(angle(R_f));                                               % Phase angle of r(t)
R_f_norm = abs(R_f) ./ N_r_t;                                              % Normalized R(f)
% PSD of r_{i}(t)
[P_r_t, f_r_t] = periodogram(r_t, hann(length(r_t), 'periodic'), [], fs);
% 99% Occupied Bandwidth of r_{i}(t) 
[bw_r_t, flo_r_t, fhi_r_t, pow_r_t] = obw(P_r_t, f_r_t);

% % Bandwidth for H_{i,j}(f)
bw_flo = mean(flo_r_t);                                                    % Lower frequency of bandwidth [Hz]
bw_fhi = mean(fhi_r_t);                                                    % Higher frequency of bandwidth [Hz]
bw_sig = bw_fhi - bw_flo;                                                  % Bandwidth of signal [Hz]
freq_car = (bw_fhi + bw_flo) / 2;                                          % *Carrier frequency [Hz]

%% ORIGINAL H(f)

[xxx, ind_bw_flo] = min(abs(freq - bw_flo));
[xxx, ind_bw_fhi] = min(abs(freq - bw_fhi));

% Indexes for H(f):
ind0_H_f = ind_bw_flo; 
ind1_H_f = ind_bw_fhi;
ind2_H_f = size(r_t, 1) - ind1_H_f;    
ind3_H_f = size(r_t, 1) - ind0_H_f; 

% Allocating memory for mapping functions
H_0_f = zeros(size(r_t));                                      % Mapping Transfer functions: H_{0}(f)

% First Half (1:(fs / 2))
H_0_f(ind0_H_f:ind1_H_f, :) = R_f(ind0_H_f:ind1_H_f, :) ./ E_f(ind0_H_f:ind1_H_f);
% Second Half ((fs / 2):fs)
H_0_f(ind3_H_f:-1:ind2_H_f, :) = real(H_0_f(ind0_H_f:ind1_H_f, :)) - imag(H_0_f(ind0_H_f:ind1_H_f, :)) .* 1i;

% H(f) in the time domain
h_0_t = ifft(H_0_f, 'symmetric');


%% POSITION FOR CORRELATION RANGE

corr_delta = lambda_min / 2;

% Using Polar Coordinates
nb_elem_corr = 4;
nb_elem_pos_corr = nb_elem_corr + 1;
corr_ds = linspace(0, (2 * corr_delta), nb_elem_pos_corr);
corr_th = linspace(0, 2 * pi, nb_elem_pos_corr);

pos_corr = zeros(nb_elem_pos_corr, nb_elem_pos_corr, 2);

for ind_pos_corr = 1:size(pos_corr, 1)
    pos_corr(ind_pos_corr, :, 1) = pos_emt(1) + corr_ds .* cos(corr_th(ind_pos_corr)); 
    pos_corr(ind_pos_corr, :, 2) = pos_emt(2) + corr_ds .* sin(corr_th(ind_pos_corr)); 
end

pos_corr(:, 1, :) = [];

dx_pos_corr = abs(pos_corr(2, 1, 1) - pos_corr(1, 1, 1));
dy_pos_corr = abs(pos_corr(2, 1, 2) - pos_corr(1, 1, 2));

x_surf = pos_corr(:, :, 1);
y_surf = pos_corr(:, :, 2);

plot(x_surf, y_surf)


% % Using Grids ...
% nb_elem_corr = 20;
% corr_ds = linspace(0, 2 * corr_delta, nb_elem_corr + 1);
% corr_th = linspace(0, 2 * pi, nb_elem_corr + 1);
% 
% pos_corr = zeros(nb_elem_corr + 1, nb_elem_corr + 1, 2);
% 
% for ind_corr = 1:(nb_elem_corr + 1)
%     pos_corr(ind_corr, :, 1) = pos_emt(1) + corr_ds .* cos(corr_th(ind_corr)); 
%     pos_corr(ind_corr, :, 2) = pos_emt(2) + corr_ds .* sin(corr_th(ind_corr)); 
% end
% 
% close all
% 
% plot(pos_corr(:, :, 1), pos_corr(:, :, 2), 'k*')

%% CALCULATION THE CORRELATION RANGE

corr_range = zeros(size(pos_corr, 1:2));

for ind_corr_i = 1:nb_elem_corr
    for ind_corr_j = 1:nb_elem_corr
        % Received signal from correlation range
        pos_emt_corr = [pos_corr(ind_corr_i, ind_corr_j, 1); pos_corr(ind_corr_i, ind_corr_j, 2)];
        [r_1_t, Ez] = function_impulse_response(param_simu, pos_emt_corr, pos_rec, Display_rep_imp, time, e_t);
        R_1_f = fft(r_1_t);
        
        % Allocating memory for mapping functions
        H_1_f = zeros(size(r_1_t));                                        % Mapping Transfer functions: H_{0}(f)

        % First Half (1:(fs / 2))
        H_1_f(ind0_H_f:ind1_H_f, :) = R_1_f(ind0_H_f:ind1_H_f, :) ./ E_f(ind0_H_f:ind1_H_f);
        % Second Half ((fs / 2):fs)
        H_1_f(ind3_H_f:-1:ind2_H_f, :) = real(H_1_f(ind0_H_f:ind1_H_f, :)) - imag(H_1_f(ind0_H_f:ind1_H_f, :)) .* 1i;

        % H(f) in the time domain
        h_1_t = ifft(H_1_f, 'symmetric');
        
        corr_rho = corrcoef(h_0_t, h_1_t);
        corr_range(ind_corr_i, ind_corr_j) = corr_rho(1, 2)
    end
end

% corr_range = [corr_range; zeros(1, size(corr_range, 2))]

%%


theta = 45; % to rotate 90 counterclockwise
rot_matrix = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];

pos_corr_plot = zeros(nb_elem_corr + 1, nb_elem_corr + 1, 2);

for ind_corr_i = 1:nb_elem_corr
    for ind_corr_j = 1:nb_elem_corr
        rot_pts = [pos_corr(ind_corr_i, ind_corr_j, 1) - pos_emt(1), ...
            pos_corr(ind_corr_i, ind_corr_j, 2) - pos_emt(2)];
        
        rot_pts = rot_matrix * rot_pts';
        
        ind_xy = 1;
        pos_corr_plot(ind_corr_i, ind_corr_j, ind_xy) = rot_pts(ind_xy) + pos_emt(ind_xy);
        
        ind_xy = 2;
        pos_corr_plot(ind_corr_i, ind_corr_j, ind_xy) = rot_pts(ind_xy) + pos_emt(ind_xy);
    end
end

pos_corr_plot(:, end, :) = [];
pos_corr_plot(end, :, :) = pos_corr_plot(1, :, :);

close all

figure(1)

x_surf = pos_corr(:, :, 1);
y_surf = pos_corr(:, :, 2);

stem3(x_surf, y_surf, corr_range)

% hold on;
% surf(x_surf, y_surf, abs(corr_range))
% plot3(x_surf, y_surf, ones(size(pos_corr, 1:2)), 'k*')
% hold off;
% 
% view(0, 90)
% colorbar
% colormap jet
% 
% figure(2)
% 
% x_surf = pos_corr_plot(:, :, 1);
% y_surf = pos_corr_plot(:, :, 2);
% 
% hold on;
% surf(x_surf, y_surf, abs(corr_range))
% plot3(pos_corr(:, :, 1), pos_corr(:, :, 2), ones(size(pos_corr, 1:2)), 'k*')
% hold off;
% 
% view(0, 90)
% colorbar
% colormap jet

