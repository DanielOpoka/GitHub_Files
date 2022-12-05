
%% ------------------------------------------------------------------------
%                       TIME FREQUENCY ANALYSIS
%  ------------------------------------------------------------------------

% Extracting any useful information from the time frequency domain from 
% reverberated signals.

clear all; close all; clc; format long;
addpath("functions");
addpath("data");

% To avoid showing figures
set(0,'DefaultFigureVisible','off')
        

%% ------------------------------------------------------------------------
%                         SIMULATION (VALENTIN'S CODE)
%  ------------------------------------------------------------------------

% Definition of simulation of a Reverberant Chamber.

%% PARAMETERS

% Impulse parameters (For TFA)

signal_type = "wavelet";

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

d_before = d_impulse * 2;                                                  % delay before impulse [s]
d_after = 500e-9;                                                          % Delay after impulse [s]
A_impulse = 1;                                                             % Amplitude of the impulse [V/m]

% Cavity parameters
Lx=1.9; Ly=1.9;                                                            % Cavity dimensions [m,m]
sigma=(6.5)*1e2;                                                           % Conductivity of the walls [S/m]
coef_losses_ok=1;                                                          % If=1: losses, if=0: PEC condition 
coef_regular_chao=2;                                                       % If=1: regular cavity, if=2: chaotic cavity 

% WARNING: lossless chaotic cavity not coded (possible: lossless or with losses reglar cavity, or with losses chaotic cavity)
[time, e_t, param_simu] = function_parameters_simu(Lx,Ly,d_impulse,A_impulse,f0,d_before,d_after,coef_regular_chao,coef_losses_ok,sigma, signal_type);


%% FOR SAMPLING TIME

% Wave parameters
c0 = 3e8;                                                                  % Speed of light [m/s]
f_impulse = 1 / d_impulse;                                                 % Frequency of impulse [Hz]
lambda_min = c0 / (f0 + 2 * f_impulse);                                    % Minimum wavelength of signal [m]
fact_delta=10;                                                             % Discretization []
delta=lambda_min/fact_delta;                                               % Discretized space step (dx, dy) [m]
dt=delta/c0/sqrt(2);                                                       % Discretized time step [s]
fs = 1 / dt;                                                               % Sampling Frequency [Hz]


%% IMPULSE RESPONSE

% Receiver and emitter locations
pos_rec = [1.0, 0.8, 1.2; 1.4, 1.4, 1.4];                                  % Position of receiver(s) [m]
pos_emt = [0.5375; 0.5375];                                                % Position of the emitter(s) [m]

% Display of the impulse respone
Display_rep_imp = 0;                                                       % If=1: display, if=0: no display
% Simulation of revereberated signal
[r_t, Ez] = function_impulse_response(param_simu, pos_emt, pos_rec, Display_rep_imp, time, e_t );

% Frequency Range
df = 1 / time(end);                                                        % Frequency Step [Hz]
freq = [0:length(time)-1] * df;                                            % Frequency [Hz]


%% ------------------------------------------------------------------------
%                          TIME-FREQUENCY ANALYSIS
%  ------------------------------------------------------------------------

% Analysis of the time-frequency domain of signals through FFT and STFT.

% % Reverberation Time
load('data/ReverberationTime/tau_RT_data_1-100CR_Pos4_CRC2-750.mat');
tau_RT = mean(tau_RT_end);                                                 % Reverberation Time [s] 

%% FOURIER TRANSFORM

% For Zero Padding:
cond_r_t_ZP = 0;                                                           % Conditioner for zero-padding or not [0/1]

if cond_r_t_ZP == 1
    n_ZP = 2;                                                              % Number of total length of the zero-padding [#]
    
    % Time range for zero-padded r_t signals
    number_of_time_steps = round(((d_before + d_after) * n_ZP) / dt);
    time_r_t_ZP = dt * [0:number_of_time_steps - 1].';
    
    % Frequency range for zero-padded r_t signals
    df_r_t_ZP = 1 / time_r_t_ZP(end);
    freq_r_t_ZP = [0:length(time_r_t_ZP) - 1] * df_r_t_ZP;         
   
    e_t_ZP = zeros((size(e_t, 2) * n_ZP) + 1, size(e_t, 1));               % Emitted with zero padding
    e_t_ZP(1:size(e_t, 2)) = e_t(:);
    e_t = e_t_ZP';
    
    r_t_ZP = zeros((size(r_t, 1) * n_ZP) + 1, size(r_t, 2));               % Received signal with zero padding
    for r_t_i = 1:size(r_t, 2)
        r_t_ZP(1:size(r_t, 1) , r_t_i) = r_t(:, r_t_i);
    end
    r_t = r_t_ZP;

    time = time_r_t_ZP;
    freq = freq_r_t_ZP;
end

% % Fourier Transform: e(t) -> E(f)
% For Normalization:     
N_e_t = nnz(e_t);                                                          % Number of elements != 0 in e_t
% No zero-padding needed as it already has it.
E_f = fft(e_t)';                                                           % FFT of e(t)
ph_E_f = unwrap(angle(E_f));                                               % Phase angle of e(t)
E_f_norm = abs(E_f) ./ N_e_t;                                              % Normalized E(f)

% PSD of e_t
[P_e_t, f_e_t] = periodogram(e_t, hann(length(e_t), 'periodic'), [], fs);
% 99% Occupied Bandwidth of e_t
[bw_e_t, flo_e_t, fhi_e_t, pow_e_t] = obw(P_e_t, f_e_t);

% [xxx, ind_flo_e_t] = min(abs(f - flo_e_t));
% [xxx, ind_fhi_e_t] = min(abs(f - fhi_e_t));

% % Fourier Transform: r(t) -> R(f)
% For Normalization
abs_r_t_cond = abs(r_t) > 0;                                               % Conditioner if abs(r_t) is greater than 0
N_r_t = sum(abs_r_t_cond == 1);                                            % Number of elements != 0 in r_t

R_f = fft(r_t);                                                            % FFT of r(t)
ph_R_f = unwrap(angle(R_f));                                               % Phase angle of r(t)
R_f_norm = abs(R_f) ./ N_r_t;                                              % Normalized R(f)

% PSD of r_t
sum_flo_r_t = 0;
sum_fhi_r_t = 0;
num_r_t = size(pos_rec, 2);

for ind_r_t = 1:num_r_t
    [P_r_t, f_r_t] = periodogram(r_t(:, ind_r_t), hann(length(r_t(:, ind_r_t)), 'periodic'), [], fs);
    % 99% Occupied Bandwidth of r_t(i)
    [bw_r_t, flo_r_t, fhi_r_t, pow_r_t] = obw(P_r_t, f_r_t);

    sum_flo_r_t = sum_flo_r_t + flo_r_t;
    sum_fhi_r_t = sum_fhi_r_t + fhi_r_t;
end
% Obtained the lowest and highest frequency range index by using the
% average from the received signals 99% occupied bandwidth.
avg_flo_r_t = sum_flo_r_t / num_r_t;
avg_fhi_r_t = sum_fhi_r_t / num_r_t;

% [xxx, ind_flo_r_t] = min(abs(f - sum_flo_r_t));
% [xxx, ind_fhi_r_t] = min(abs(f - sum_fhi_r_t));

% Deconvolution
% Bandwidth of the Signal
% Using 99% occupied bandwidth
bw_flo = flo_e_t;                                                          % Lower frequency of bandwidth [Hz]
bw_fhi = fhi_e_t;                                                          % Higher frequency of bandwidth [Hz]
bw_sig = bw_fhi - bw_flo;                                                  % Bandwidth of signal [Hz]
freq_car = (bw_fhi + bw_flo) / 2;                                          % Carrier frequency [Hz]

% % Using wavelet's bandwidth
% bw_flo = f0 - (2 / d_impulse);
% bw_fhi = f0 + (2 / d_impulse);

% % Using chirp's bandwidth
% bw_flo = 2e9;
% bw_fhi = 3e9;

% Bandwidth index of the Signal
[xxx, ind_bw_flo] = min(abs(freq - bw_flo));
[xxx, ind_bw_fhi] = min(abs(freq - bw_fhi));
[xxx, ind_freq_car] = (min(abs(freq - freq_car)));

% Room Impulse caused by Reverberation in VC
H_f = zeros(size(r_t));
E_f_recon = zeros(size(r_t));

% % For H_f-database:
% H_f_data = [freq', H_f];
% freq_chirp = freq(ind_bw_flo:ind_bw_fhi);
% E_f_chirp = abs(E_f);
% E_f_chirp = E_f_chirp(ind_bw_flo:ind_bw_fhi);
% % For saving the data
% data_name = "data/E_f_chirp_Alex/E_f_chirp" + num2str(d_impulse) + ".txt";
% writematrix(E_f_chirp, data_name);
% type data_name;

% Indexs for H(f):
ind0_H_f = ind_bw_flo; 
ind1_H_f = ind_bw_fhi;
ind2_H_f = length(r_t) - ind1_H_f;    
ind3_H_f = length(r_t) - ind0_H_f; 

for r_t_i = 1:size(r_t, 2)
    % First Half (1:(fs / 2))
    H_f(ind0_H_f:ind1_H_f, r_t_i) = R_f(ind0_H_f:ind1_H_f, r_t_i) ./ E_f(ind0_H_f:ind1_H_f);
    % Second Half ((fs / 2):fs)
    H_f(ind3_H_f:-1:ind2_H_f, r_t_i) = real(H_f(ind0_H_f:ind1_H_f, r_t_i)) - imag(H_f(ind0_H_f:ind1_H_f, r_t_i)) .* 1i;

    figure;
    plot(abs(H_f,:,1))
    bfbffb


    % Reconstructing E_f from H_f
    E_f_recon(ind0_H_f:ind1_H_f, r_t_i) = (R_f(ind0_H_f:ind1_H_f, r_t_i) ./ H_f(ind0_H_f:ind1_H_f, r_t_i));
    E_f_recon(ind3_H_f:-1:ind2_H_f, r_t_i) = E_f_recon(ind0_H_f:ind1_H_f, r_t_i);
end

ph_H_f = unwrap(angle(H_f));                                               % Phase angle of H(f)
h_t = real(ifft(H_f));                                                     % H(f) in time          
N_h_t = nnz(h_t) ./ size(h_t, 2);                                          % Number of useful elements in h(t)
H_f_norm = abs(H_f) ./ N_h_t;                                              % Normalized H(f)

e_t_recon = real(ifft(E_f_recon));
e_t_recon = flip(e_t_recon, 1);

%% INDEXING FOR TFA

% Short Time Fourier Transform Parameters
tau_stft = d_impulse;                                                      % Time window size [s]
% Time win_stftdow in term of elements with respect to the given time data:
ind_stft = find(min(abs(time - tau_stft)) == abs(time - tau_stft));
over_per = 0.95;                                                           % Overlap percentage [%]
over_len = round(ind_stft * over_per);                                     % Overlap Length [#]
win_stft = hann(ind_stft, 'periodic');                                     % Time window [s]

% Wigner-Ville Distribution Parameters
if signal_type == "chirp_wavelet"
    tau_wvd = 50e-9;
    ind_wvd = find(min(abs(time - tau_wvd)) == abs(time - tau_wvd));
    % To make sure WVD windows are odd
    if mod(ind_wvd, 2) == 0
        ind_wvd = ind_wvd + 1;
    end
    win_wvd = hann(ind_wvd, 'periodic');

    twin_wvd = win_wvd;                                                    % Time window [s]
    fwin_wvd = win_wvd;                                                    % Frequency window [Hz]

    freq_pts_wvd = ind_wvd + 2;                                            % Frequency points [#]
    time_pts_wvd = ind_wvd + 1;                                            % Time points [#]
elseif signal_type == "wavelet"
    freq_pts_wvd = round(length(time) * (1 / 2)) + 1;                      % Frequency points [#]
    time_pts_wvd = round(length(time) * (1 / 2));                          % Time points [#]
    % To make sure time points is even
    if mod(time_pts_wvd, 2) == 1
        time_pts_wvd = time_pts_wvd + 1;
    end
    tau_wvd = tau_stft;
    ind_wvd = find(min(abs(time - tau_wvd)) == abs(time - tau_wvd));
    
    twin_wvd = win_stft;                                                   % Time window [s]
    fwin_wvd = win_stft;                                                   % Frequency window [Hz]
end


%% ------------------------------------------------------------------------
%                           PLOTTING OF TFA
%  ------------------------------------------------------------------------

% Plotting of all relevant data.

close all;

n_fontsize = 14;
n_fig = 0; % For an ever increasing number of plots
signal_input_label = signal_type;

new_dir_signal = "images/TFA/" + signal_input_label;
mkdir(new_dir_signal);

new_dir_TFA = new_dir_signal + "/t_end=" + num2str(round(time(end), 9) * 1e9) + "ns";
mkdir(new_dir_TFA);

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

figure_title = new_dir_TFA + "/cavity-" + num2str(size(pos_emt, 2)) + "_emt(s)-" + ...
    num2str(size(pos_rec, 2)) + "_rec(s).png";
saveas(gcf, figure_title)


%% TFA OF e(t)

% Input for TFA plot of e_t
tlim_t0 = d_before - (d_impulse / 2);
tlim_t1 = d_before + (1.5 * d_impulse);
tlim = [tlim_t0, tlim_t1];

ind_flim = 2;
flim_f0 = ((bw_flo + bw_fhi) / 2) - ind_flim * (bw_fhi - bw_flo) / 2;
flim_f1 = ((bw_flo + bw_fhi) / 2) + ind_flim * (bw_fhi - bw_flo) / 2;
flim = [flim_f0, flim_f1];

title_label = ' e(t)';
ind_win_stft = [min(find(e_t)) + 1, min(find(e_t)) + ind_stft];
signal = e_t;
FFT_norm = E_f_norm;

n_fig = n_fig + 1;
figure(n_fig);

STFT_plot(time, freq, fs, ind_win_stft, win_stft, signal, FFT_norm, over_len, ...
    tlim, flim, n_fontsize, title_label)

figure_title = new_dir_TFA + "/STFT-e_t-t_win=" + num2str(tau_stft * 1e9) + "ns.png";
saveas(gcf, figure_title)

% !!! WARNING: DO NOT USE, TOO COMPUTATIONALLY EXPENSIVE !!!
% n_fig = n_fig + 1;
% figure(n_fig);
% 
% ind_win_wvd = [min(find(e_t)) + 1, min(find(e_t)) + ind_wvd];
% 
% WVD_plot(time, freq, fs,  ind_win_wvd, signal, FFT_norm, ...
%     twin_wvd, tlim, flim, n_fontsize, title_label)
% 
% figure_title = new_dir_TFA + "/WVD-e_t-t_win=" + num2str(tau_wvd * 1e9) + "ns.png";
% saveas(gcf, figure_title)


%% TFA OF EVERY r(t)
% !!! CAUTION !!!
% Do not use this if more than 3 receivers, it will crash MATLAB

% Input for TFA plot of r_t & h_t
tlim = [0, time(end)];
ind_win_stft = [1, ind_stft];
ind_win_wvd = [1, ind_wvd];

for r_t_i = 1:size(r_t, 2)
    % For r(t)
    signal = normalize(r_t(:, r_t_i), 'norm', Inf);
    FFT_norm = R_f_norm(:, r_t_i);
    title_label = "r_" + num2str(r_t_i) + "(t)";
    
    n_fig = n_fig + 1;
    figure(n_fig);
    STFT_plot(time, freq, fs, ind_win_stft, win_stft, signal, FFT_norm, ...
        over_len, tlim, flim, n_fontsize, title_label)
    
    figure_title = new_dir_TFA + "/STFT-" + ...
        title_label + "-t_win=" + num2str(tau_stft * 1e9) + "ns.png";
    saveas(gcf, figure_title)
    
    
%     !!! WARNING: DO NOT USE, TOO COMPUTATIONALLY EXPENSIVE !!!
%     n_fig = n_fig + 1;
%     figure(n_fig);
% 
%     WVD_plot(time, freq, fs,  ind_win_wvd, signal, FFT_norm, ...
%         twin_wvd, fwin_wvd, freq_pts_wvd, time_pts_wvd, tlim, flim, n_fontsize, title_label)
% 
%     figure_title = new_dir_TFA + "/WVD-" + ...
%         title_label + "-t_win=" + num2str(tau_wvd * 1e9) + "ns.png";
%     saveas(gcf, figure_title)
    
    % For h(t)
    signal = normalize(h_t(:, r_t_i), 'norm', Inf);
    FFT_norm = H_f_norm(:, r_t_i);
    title_label = " h_" + num2str(r_t_i) + "(t)";
    
    n_fig = n_fig + 1;
    figure(n_fig);
    STFT_plot(time, freq, fs, ind_win_stft, win_stft, signal, FFT_norm, ...
        over_len, tlim, flim, n_fontsize, title_label)

    figure_title = new_dir_TFA + "/STFT-" + title_label + "-t_win=" + num2str(tau_stft * 1e9) + "ns.png";
    saveas(gcf, figure_title)
    
        n_fig = n_fig + 1;
    figure(n_fig);

    WVD_plot(time, freq, fs,  ind_win_wvd, signal, FFT_norm, ...
        twin_wvd, fwin_wvd, freq_pts_wvd, time_pts_wvd, tlim, flim, n_fontsize, title_label)

    figure_title = new_dir_TFA + "/WVD-" + ...
        title_label + "-t_win=" + num2str(tau_wvd * 1e9) + "ns.png";
    saveas(gcf, figure_title)
end

%% PLOTS OF PHASES

close all

n_fig = n_fig + 1;
figure(n_fig);
set(gcf, 'Position', get(0, 'Screensize'));
tiledlayout(3,1);

nexttile

hold on;
grid on;
plot(freq, ph_E_f, 'r');
hold off;

xlim([flim_f0, flim_f1]);

% Labelling
set(gca,'fontsize',n_fontsize);
xlabel("Frequency [Hz]");
ylabel("Phase Angle [rad]");
title("Unwrapped Phase of E(f)");

nexttile

hold on;
grid on;
plot(freq, ph_R_f(:, 1), 'g');
hold off;

xlim([flim_f0, flim_f1]);

% Labelling
set(gca,'fontsize',n_fontsize);
xlabel("Frequency [Hz]");
ylabel("Phase Angle [rad]");
title("Unwrapped Phase of R_{1}(f)");

nexttile

hold on;
grid on;
plot(freq, ph_H_f(:, 1), 'b');
hold off;

xlim([flim_f0, flim_f1]);

% Labelling
set(gca,'fontsize',n_fontsize);
xlabel("Frequency [Hz]");
ylabel("Phase Angle [rad]");
title("Unwrapped Phase of H_{1}(f)");

figure_title = new_dir_TFA + "/" + signal_type + "-unwrapped phase.png";
saveas(gcf, figure_title)

