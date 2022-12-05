
%% ------------------------------------------------------------------------
%                       ESTIMATION OF EMITTER POSITION
%  ------------------------------------------------------------------------

% 1. SIMULATION (VALENTIN'S CODE)
% 2. MAPPING AND RECONSTRUCTION
% 3. ESTIMATING THE POSITION

% i: Refers to a receiver.
% j: Refers to a mapped position.

clear all; close all; clc; format long;
addpath("functions");
addpath("data/H_fs");

new_dir_0 = "images/EstEmtPos"; 
mkdir(new_dir_0);

nb_fig = 0;
n_fontsize = 12;
       

%% ------------------------------------------------------------------------
%                       1. SIMULATION (VALENTIN'S CODE)
%  ------------------------------------------------------------------------

% Definition of simulation of a Reverberant Chamber.

%% PARAMETERS

% Emitter parameters
signal_type = "wavelet";

nb_signals = 1;

if signal_type == "wavelet"
    d_impulse = 8e-9;                                                      % Duration of the impulse [s]
    f0 = 2.4e9;                                                            % Carrier frequency of the impulse [Hz]
elseif signal_type == "chirp_wavelet" || signal_type == "chirp_sinusoidal"
    d_impulse = 200e-9;                                                    % (Chirp) Duration of the impulse [s]
%     d_impulse = 8e-9;                                                    % (Chirp with same BW as wavelet) Duration of the impulse [s]
    f0 = 2.4e9;                                                            % Carrier frequency of the impulse [Hz]
elseif signal_type == "sinc"
    d_impulse = 50e-9;                                                     % (Sinc) Duration of the impulse [s]
    f0 = 8e9;                                                              % Carrier frequency of the impulse [Hz]
end

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
pos_rec = [1.0, 0.8, 1.2; 1.4, 1.4, 1.4];                                  % Alex's location of the receiver [m]
pos_emt = [0.5389; 0.5376];                                                % Alex's location of the emitter [m]

% Display of the impulse respone
Display_rep_imp = 0;                                                       % If=1: display, if=0: no display
[r_t, Ez] = function_impulse_response(param_simu, pos_emt, pos_rec, Display_rep_imp, time, e_t );

% Frequency Range
df = 1 / time(end);
freq = [0:length(time)-1] * df;

new_dir_test = new_dir_0 + "/test-" + num2str(nb_signals) + "_" + signal_type + ...
    "-f0=" + num2str(f0 .* 1e-9) + "GHz_d=" + num2str(d_impulse .* 1e9) + "ns-" + ...
    "emt_x=" + num2str(round(pos_emt(1), 3)) + "m_y=" + num2str(round(pos_emt(2), 3)) + "m";
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
%                  2. MAPPING AND RECONSTRUCTION
%  ------------------------------------------------------------------------

% Mapping of different H_{i,j}(f)s by using parameters from r_{i}(t).
% After creating H_{i,j}(f)s, reconstruct the ~E_{i,j}(f)s.

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


%% PLOTTING OF EMITTED AND RECEIVED SIGNALS

nb_fig = 0;
close all;

nb_fig = nb_fig + 1;
figure(nb_fig);

set(gcf, 'Position', get(0, 'Screensize'));

ax(1) = subplot(2, 2, 1);
hold on;
grid on;
plot(time .* 1e9, e_t);
hold off;

xlabel('time [ns]');
ylabel('Signal Amplitude [V/m]');
title({num2str(size(pos_emt, 2)) + " Emitted Signal(s) (" + signal_type + "): e(t)"});

set(gca,'fontsize',n_fontsize);

ax(2) = subplot(2, 2, 2);
hold on;
grid on;
plot(freq .* 1e-9, E_f_norm);
hold off;

xlabel('Frequency [GHz]');
ylabel('Normalized Amplitude [V/m]');
title({"FFT of " + num2str(size(pos_rec, 2)) + " Emitted Signal(s): E(f)"});

xlim([bw_flo, bw_fhi] .* 1e-9);

set(gca,'fontsize',n_fontsize);

ax(3) = subplot(2, 2, 3);
hold on;
grid on;
plot(time .* 1e9, r_t);
hold off;

xlabel('time [ns]');
ylabel('Signal Amplitude [V/m]');
title({num2str(size(pos_rec, 2)) + " Received Signal(s): r_{i}(t)"});

set(gca,'fontsize',n_fontsize);

ax(4) = subplot(2, 2, 4);
hold on;
grid on;
plot(freq .* 1e-9, R_f_norm);
hold off;

xlabel('Frequency [GHz]');
ylabel('Normalized Amplitude [V/m]');
title({"FFT of " + num2str(size(pos_rec, 2)) + " Received Signal(s): R_{i}(f)"});

xlim([bw_flo, bw_fhi] .* 1e-9);

set(gca,'fontsize',n_fontsize);

figure_title = new_dir_test + "/signal(s)-" + num2str(size(pos_emt, 2)) + "_emt(s)-" + ...
    num2str(size(pos_rec, 2)) + "_rec(s).png";
saveas(gcf, figure_title)


%% POSITIONS FOR TRANSFER FUNCTIONS

map_init = [0.53; 0.53];
map_rows = 10;
map_cols = 10;
map_size = map_rows * map_cols;
map_dist = (lambda_min / 2) + 0.01;

% Position of the emitters for mapping
pos_emt_map = zeros(2, map_size);
ind_map = 1;

for ind_cols = 1:map_cols
    for ind_rows = 1:map_rows
        pos_emt_map(2, ind_map) = map_init(1, 1) + ((ind_rows - 1) * map_dist);
        pos_emt_map(1, ind_map) = map_init(2, 1) + ((ind_cols - 1) * map_dist); 
        ind_map = ind_map + 1;
    end
end


%% PLOTTING OF MAPPING

close all

openfig("data/CavitySpace/cavity_space_empty.fig");
figure(1);

close all;

n_fontsize = 12;
time_cvs = 0;

% openfig("data/CavitySpace/cavity_space_" + num2str(time_cvs) + "ns.fig");
openfig("data/CavitySpace/cavity_space_empty.fig");
figure(1);

% For plotting out a mapping section:
% plot(position_receiver_x, position_receiver_y, 'rx', 'MarkerSize', 30, 'LineWidth', 3)
% legend('', 'Receivers', 'Predefined Area', 'Location', 'northeast');

pos_emt_map_x = pos_emt_map(1, :);
pos_emt_map_y = pos_emt_map(2, :);

pos_rec_x = pos_rec(1, :);
pos_rec_y = pos_rec(2, :);

hold on;
plot(pos_emt_map_x, pos_emt_map_y, 'k.', 'MarkerSize', 20, 'LineWidth', 2)
plot(pos_rec_x, pos_rec_y, 'rx', 'MarkerSize', 20, 'LineWidth', 2)
hold off;

set(gca,'fontsize',n_fontsize);

xlabel("x-axis [m]");
ylabel("y-axis [m]");
title({"Cavity with " + num2str(map_size) + ...
    " Mapping Emitter(s) and " + num2str(size(pos_rec, 2)) + " Receiver(s)"});

new_dir_map = new_dir_test + "/mapping-init(x=" + num2str(round(map_init(1), 3)) + ...
    "m_y=" + num2str(round(map_init(2), 3)) + "m)-rows=" + num2str(map_rows) + ...
    "_cols=" + num2str(map_cols) + "-dist=" + num2str(round(map_dist, 3)) + "m";
mkdir(new_dir_map);

figure_title = new_dir_map + "/cavity-" + num2str(map_size) + "_emt(s)_map-" + ...
    num2str(size(pos_rec, 2)) + "_rec(s).png";
saveas(gcf, figure_title)

% openfig("data/CavitySpace/cavity_space_empty.fig");
% figure(2);

pos_emt_map_x = pos_emt_map(1, :);
pos_emt_map_y = pos_emt_map(2, :);

pos_rec_x = pos_rec(1, :);
pos_rec_y = pos_rec(2, :);

text_dist = map_dist / 10;

hold on;
% plot(pos_rec_x, pos_rec_y, 'rx', 'MarkerSize', 20, 'LineWidth', 2)
for ind_emt_map = 1:map_size
    text(pos_emt_map_x(ind_emt_map) - text_dist, pos_emt_map_y(ind_emt_map), ...
        num2str(ind_emt_map), 'FontSize', n_fontsize, 'HorizontalAlignment','right')
end
hold off;

set(gca,'fontsize',n_fontsize);

xlabel("x-axis [m]");
ylabel("y-axis [m]");
title({"Zoom-in of Cavity with " + num2str(map_size) + ...
    " Mapping Emitter(s) and " + num2str(size(pos_rec, 2)) + " Receiver(s)"});

xlim([min(pos_emt_map_x) - map_dist, max(pos_emt_map_x) + map_dist]);
ylim([min(pos_emt_map_y) - map_dist, max(pos_emt_map_y) + map_dist]);

figure_title = new_dir_map + "/cavity_zoom_in-" + num2str(map_size) + "_emt(s)_map-" + ...
    num2str(size(pos_rec, 2)) + "_rec(s).png";
saveas(gcf, figure_title)

%% CREATING THE EMITTED SIGNAL FOR MAPPING

% lambda_min_map = (c0) / (freq_car + 2 * bw_sig);

d_impulse_map = 1 / bw_sig;                                                % Bandwidth frequency for sinc

% Emitted Signal for Mapping: sinc
e_t_map = sinc(pi .* bw_sig .* (time - (time(end) / 2))) .* ...
    cos(2 .* pi .* freq_car .* (time - (time(end) / 2)));

e_t_map = e_t_map * A_impulse;

E_f_map = fft(e_t_map);                                                    % FFT of e_{map}(t)
E_f_map_norm = abs(E_f_map) ./ nnz(e_t_map);                               % Normalized FFT of e_{map}(t)


%% PLOTTING OF EMITTED SIGNAL FOR MAPPING

nb_fig = 0;
close all;

nb_fig = nb_fig + 1;
figure(nb_fig);

zoom_type = "zoom-in";
% zoom_type = "zoom-out";

set(gcf, 'Position', get(0, 'Screensize'));

ax(1) = subplot(1, 2, 1);
hold on;
grid on;
plot(time .* 1e9, e_t_map);
hold off;

xlabel('time [ns]');
ylabel('Signal Amplitude [V/m]');
title({num2str(size(pos_emt, 2)) + " Mapping Emitted Signal(s) (" + signal_type + "): e_{map}(t)"});

if zoom_type == "zoom-in"
    xlim([(time(end) / 2) - (d_impulse), (time(end) / 2) + (d_impulse)] .* 1e9);
elseif zoom_type == "zoom-out"
    xlim([0, time(end)] .* 1e9);
end

set(gca,'fontsize',n_fontsize);

ax(2) = subplot(1, 2, 2);
hold on;
grid on;
plot(freq .* 1e-9, E_f_map_norm);
xline(bw_flo * 1e-9, 'k--');
xline(bw_fhi * 1e-9, 'k--');
hold off;

xlabel('Frequency [GHz]');
ylabel('Normalized Amplitude [V/m]');
title({"FFT of " + num2str(size(pos_rec, 2)) + " Mapping Emitted Signal(s): E_{map}(f)"});

if zoom_type == "zoom-in"
    xlim([bw_flo, bw_fhi] .* 1e-9);
elseif zoom_type == "zoom-out"
    xlim([bw_flo - (2 * bw_sig), bw_fhi + (2 * bw_sig)] .* 1e-9);
end

set(gca,'fontsize',n_fontsize);

figure_title = new_dir_map + "/signal-emt_map-" + zoom_type + ".png";
saveas(gcf, figure_title)


%% DEFINING H_{i,j}(f)s AND ~E_{i,j}(f)s

[xxx, ind_bw_flo] = min(abs(freq - bw_flo));
[xxx, ind_bw_fhi] = min(abs(freq - bw_fhi));

% Indexs for H(f):
ind0_H_f = ind_bw_flo; 
ind1_H_f = ind_bw_fhi;
ind2_H_f = size(r_t, 1) - ind1_H_f;    
ind3_H_f = size(r_t, 1) - ind0_H_f; 

% Allocating memory for mapping functions
E_f_map = fft(e_t_map);                                                    % FFT of e_{map}(t): E_{map}(f)
R_f_map = zeros([size(r_t), map_size]);                                    % FFT of r_{map}(t): R_{map}(f)
H_f_map = zeros([size(r_t), map_size]);                                    % Mapping Transfer functions: H_{i,j}(f)
E_f_recon = zeros([size(r_t), map_size]);                                  % Reconstructed FFT: ~E(f)
e_t_recon = zeros([size(r_t), map_size]);                                  % Reconstructed Signal: ~e(t) 

% Display of the impulse response
Display_rep_imp = 0;                                                       % If=1: display, if=0: no display

for ind_map_j = 1:map_size 
    [r_t_map, Ez] = function_impulse_response(param_simu, pos_emt_map(:, ind_map_j), pos_rec, Display_rep_imp, time, e_t_map);
    R_f_map(:, :, ind_map_j) = fft(r_t_map);                               % FFT of r_{map}(t)
    
    % First Half (1:(fs / 2))
    H_f_map(ind0_H_f:ind1_H_f, :, ind_map_j) = R_f_map(ind0_H_f:ind1_H_f, :, ind_map_j) ./ E_f_map(ind0_H_f:ind1_H_f);
    % Second Half ((fs / 2):fs)
    H_f_map(ind3_H_f:-1:ind2_H_f, :, ind_map_j) = real(H_f_map(ind0_H_f:ind1_H_f, :, ind_map_j)) - imag(H_f_map(ind0_H_f:ind1_H_f, :, ind_map_j)) .* 1i;

    % Reconstructing ~E(f) from H_{map}(f)
    E_f_recon(ind0_H_f:ind1_H_f, :, ind_map_j) = (R_f(ind0_H_f:ind1_H_f, :) ./ H_f_map(ind0_H_f:ind1_H_f, :, ind_map_j));
    E_f_recon(ind3_H_f:-1:ind2_H_f, :, ind_map_j) = E_f_recon(ind0_H_f:ind1_H_f, :, ind_map_j);

    e_t_recon(:, :, ind_map_j) = real(ifft(E_f_recon(:, :, ind_map_j)));
    counter = num2str(ind_map_j) + " / " + num2str(map_size)
end


%% ------------------------------------------------------------------------
%                     3. ESTIMATING THE POSITION
%  ------------------------------------------------------------------------

% Estimating the emitter's positions, by using the following:
% - Correlation Maximization
% - Entropy Minimization

%% CALCULATIONS

corr_rho_sum = zeros(size(r_t, 2), map_size);
ent_sum = zeros(size(r_t, 2), map_size);

for ind_r_t_i = 1:size(r_t, 2)
    for ind_map_j = 1:map_size
        ind_unequal = (1:map_size ~= ind_map_j);
        for ind_map_k = 1:map_size
            if ind_unequal(ind_map_k) == 1
                corr_rho = corrcoef(e_t_recon(:, ind_r_t_i, ind_map_j), e_t_recon(:, ind_r_t_i, ind_map_k));
                corr_rho = corr_rho(1, 2);
                corr_rho_sum(ind_r_t_i, ind_map_j) = corr_rho_sum(ind_r_t_i, ind_map_j) + (corr_rho);
            end
        end
        
        % Obtaining probabilities
        probs_ent = histcounts(e_t_recon(:, ind_r_t_i, ind_map_j), size(r_t, 1));
        % Removing any 0s
        probs_ent(probs_ent == 0) = [];
        % Normalizing to 1
        probs_ent = probs_ent ./ size(e_t_recon, 1);
        % Calculating the entropy
        ent_sum(ind_r_t_i, ind_map_j) = -sum(probs_ent .* log10(probs_ent));
    end
end

corr_rho_sum_max = max(corr_rho_sum);
[xxx, ind_corr_max] = max(corr_rho_sum_max);

ent_sum = sum(ent_sum);
[xxx, ind_ent_min] = min(ent_sum);


%% PLOTTING THE CORRELATION MAPPING

min_corr_map_X = min(pos_emt_map(1, :)) - (map_dist / 2);
max_corr_map_X = max(pos_emt_map(1, :)) + (map_dist / 2);

min_corr_map_Y = min(pos_emt_map(2, :)) - (map_dist / 2);
max_corr_map_Y = max(pos_emt_map(2, :)) + (map_dist / 2);

corr_map_X = linspace(min_corr_map_X, max_corr_map_X, map_rows + 1);
corr_map_Y = linspace(min_corr_map_Y, max_corr_map_Y, map_cols + 1);
[corr_map_X, corr_map_Y] = meshgrid(corr_map_X, corr_map_Y);

matrix_corr_max = reshape(corr_rho_sum_max, map_rows, map_cols);
surf_corr_max = ones(size(corr_map_X)) .* min(corr_rho_sum_max);
surf_corr_max(1:map_rows, 1:map_cols) = matrix_corr_max;

matrix_ent_min = reshape(ent_sum, map_rows, map_cols); 
surf_ent_min = ones(size(corr_map_X)) .* min(ent_sum);
surf_ent_min(1:map_rows, 1:map_cols) = matrix_ent_min;

nb_fig = 0;

close all

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
axis equal;
plot3(pos_emt_map(1, :), pos_emt_map(2, :),  repmat(max(corr_rho_sum_max), 1, map_size), 'k.', 'MarkerSize', marker_size, 'LineWidth', line_width)
for ind_emt_map = 1:map_size
    text(pos_emt_map_x(ind_emt_map) - text_dist, pos_emt_map_y(ind_emt_map), ...
        max(corr_rho_sum_max), num2str(ind_emt_map), ...
        'FontSize', n_fontsize, 'HorizontalAlignment','right')
end
surf(corr_map_X, corr_map_Y, surf_corr_max)
plot3(pos_emt(1), pos_emt(2), max(corr_rho_sum_max), 'g*', 'MarkerSize', marker_size, 'LineWidth', line_width)
hold off;
view(0, 90)

xlim([min_corr_map_X, max_corr_map_X]);
ylim([min_corr_map_Y, max_corr_map_Y]);

colorbar
colormap("jet")

title("Correlation Maximization Mapping")
xlabel("x-axis [m]");
ylabel("y-axis [m]");

figure_title = new_dir_map + "/corr_max.png";
saveas(gcf, figure_title)

nb_fig = nb_fig + 1;
figure(nb_fig)

hold on;
grid on;
axis equal;
plot3(pos_emt_map(1, :), pos_emt_map(2, :),  repmat(max(ent_sum), 1, map_size), 'k.', 'MarkerSize', marker_size, 'LineWidth', line_width)
for ind_emt_map = 1:map_size
    text(pos_emt_map_x(ind_emt_map) - text_dist, pos_emt_map_y(ind_emt_map), ...
        max(ent_sum), num2str(ind_emt_map), ...
        'FontSize', n_fontsize, 'HorizontalAlignment','right')
end
surf(corr_map_X, corr_map_Y, surf_ent_min)
plot3(pos_emt(1), pos_emt(2), max(ent_sum), 'g*', 'MarkerSize', marker_size, 'LineWidth', line_width)
hold off;
view(0, 90)

xlim([min_corr_map_X, max_corr_map_X]);
ylim([min_corr_map_Y, max_corr_map_Y]);

colorbar
colormap("jet")

title("Entropy Minimization Mapping");
xlabel("x-axis [m]");
ylabel("y-axis [m]");

figure_title = new_dir_map + "/ent_min.png";
saveas(gcf, figure_title)


%% PLOTTING OF MAPPING FUNCTIONS AND RECONSTRUCTED SIGNAL

nb_fig = 0;
close all;

alg_text = ["corr_max", "ent_min"];
ind_alg = [ind_corr_max, ind_ent_min];

for ind_plot = 1:2
    nb_fig = nb_fig + 1;
    figure(nb_fig);

    ind_plot_map = ind_alg(ind_plot);

    xlim_flo = bw_flo - (bw_sig / 2);
    xlim_fhi = bw_fhi + (bw_sig / 2);

    set(gcf, 'Position', get(0, 'Screensize'));

    ax(1) = subplot(2, 3, 1);
    hold on;
    grid on;
    plot(freq .* 1e-9, abs(R_f_map(:, :, ind_plot_map)));
    xline(bw_flo * 1e-9, 'k--');
    xline(bw_fhi * 1e-9, 'k--');
    hold off;

    xlabel('Frequency [GHz]');
    ylabel('FFT Amplitude [V/m]');
    title({"FFT of r_{i," + num2str(ind_plot_map) + "}(t): R_{i," + num2str(ind_plot_map) + "}(f)"});

    xlim([xlim_flo, xlim_fhi] .* 1e-9);

    set(gca,'fontsize',n_fontsize);

    ax(2) = subplot(2, 3, 4);
    hold on;
    grid on;
    plot(freq .* 1e-9, abs(H_f_map(:, :, ind_plot_map)));
    xline(bw_flo * 1e-9, 'k--');
    xline(bw_fhi * 1e-9, 'k--');
    hold off;

    xlabel('Frequency [GHz]');
    ylabel('FFT Amplitude [V/m]');
    title({"Mapped Transfer Function: H_{i," + num2str(ind_plot_map) + "}(f)"});

    xlim([xlim_flo, xlim_fhi] .* 1e-9);

    set(gca,'fontsize',n_fontsize);

    ax(3) = subplot(2, 3, [2, 5]);
    hold on;
    grid on;
    plot(freq .* 1e-9, abs(E_f_recon(:, :, ind_plot_map)));
    hold off;

    xlabel('Frequency [GHz]');
    ylabel('FFT Amplitude [V/m]');
    title({"Reconstruction FFT: E_{i," + num2str(ind_plot_map) + "}(f))"});

    xlim([xlim_flo, xlim_fhi] .* 1e-9);

    set(gca,'fontsize',n_fontsize);

    ax(4) = subplot(2, 3, [3, 6]);
    hold on;
    grid on;
    plot(time .* 1e9, e_t_recon(:, :, ind_plot_map));
    hold off;

    xlabel('Time [ns]');
    ylabel('Signal Amplitude [V/m]');
    title({"Reconstructed Signal: ifft(E_{i," + num2str(ind_plot_map) + "}(f)) -> e_{i," + num2str(ind_plot_map) + "}(t)"});

    xlim([0, time(end)] .* 1e9);

    set(gca,'fontsize',n_fontsize);

    figure_title = new_dir_map + "/signal-mapping_recon-j=" + alg_text(ind_plot) + ".png";
    saveas(gcf, figure_title)
end
