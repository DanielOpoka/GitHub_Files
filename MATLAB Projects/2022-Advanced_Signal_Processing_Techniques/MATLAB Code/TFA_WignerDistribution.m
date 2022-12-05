
%% ------------------------------------------------------------------------
%                            WIGNER DISTRIBUTION
%  ------------------------------------------------------------------------

% Analytical representation of an impulse response.

clear all; close all; clc; format long;
addpath("functions");
addpath("data");

%% ------------------------------------------------------------------------
%                            ANALYTICAL TESTS
%  ------------------------------------------------------------------------

%% PARAMETERS

n_steps = 1000;
time = linspace(0, 1, n_steps);
freq = linspace(-120, 120, n_steps);


%% Wigner Distribution alpha

alpha = -20;
W_alpha = WDF_alpha(alpha, time, freq);


%% Wigner Distribution h
miu = -5;
omega_c = 60;
W_miu1 = WDF_alpha(miu, time, freq - omega_c);
W_miu2 = WDF_alpha(miu, time, freq + omega_c);
W_miu3 =  WDF_alpha(miu, time, freq);
W_h = ((1 / (4 * omega_c^2)) .* W_miu1 + ...
    (1 / (4 * omega_c^2)) .* W_miu2 - ...
    (1 / (2 * omega_c^2)) .* W_miu3 .* cos(2 .* omega_c .* time));


%% Plots

close all;

n_fontsize = 14;
n_fig = 0; % For an ever increasing number of plots

n_fig = n_fig + 1;
figure(n_fig);

surf(time, freq, W_alpha, 'edgecolor', 'none')
colormap jet
view(0, 90)

xlim([time(1), time(end)])
ylim([freq(1), freq(end)])

xlabel('t', 'interpreter', 'latex')
ylabel('$\omega$', 'interpreter', 'latex')
zlabel('$W_{\alpha}(t, \omega)$', 'interpreter', 'latex')
set(gca,'fontsize',n_fontsize);


%% ------------------------------------------------------------------------
%                         SIMULATION (VALENTIN'S CODE)
%  ------------------------------------------------------------------------

% Definition of simulation of a Reverberant Chamber.

%% PARAMETERS

% Impulse parameters (For TFA)

signal_type = "wavelet";
d_impulse = 8e-9;                                                          % Duration of the impulse [s]
f0 = 2.4e9;                                                                % Carrier frequency of the impulse [Hz]

d_before = 0;                                                              % Delay before impulse [s]
d_after = 1000e-9;                                                         % Delay after impulse [s]
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

% Receiver and emitter locations:
position_receiver=[1.0, 0.8, 1.2; 1.4, 1.4, 1.4];  % Alex's location of the receiver [m]
position_emitter=[0.5375; 0.5375];                 % Alex's location of the emitter [m]

% Display of the impulse respone
Display_rep_imp = 0;                               % If=1: display, if=0: no display
[r_t, Ez] = function_impulse_response(param_simu, position_emitter, position_receiver, Display_rep_imp, time, e_t );

% Frequency Range
df = 1 / time(end);
freq = [0:length(time)-1] * df;

%% BANDWIDTH OF THE SIGNAL

% PSD of e_t
[P_e_t, f_e_t] = periodogram(e_t, hann(length(e_t), 'periodic'), [], fs);
% 99% Occupied Bandwidth of e_t
[bw_e_t, flo_e_t, fhi_e_t, pow_e_t] = obw(P_e_t, f_e_t);

% Using 99% occupied bandwidth
bw_flo = flo_e_t;                                       % Lower frequency of bandwidth [Hz]
bw_fhi = fhi_e_t;                                       % Higher frequency of bandwidth [Hz]
bw_sig = bw_fhi - bw_flo;                               % Bandwidth of signal [Hz]
freq_car = (bw_fhi + bw_flo) / 2;                       % *Carrier frequency [Hz]

% Bandwidth index of the Signal
[xxx, ind_bw_flo] = min(abs(freq - bw_flo));
[xxx, ind_bw_fhi] = min(abs(freq - bw_fhi));
[xxx, ind_freq_car] = (min(abs(freq - freq_car)));

%% WIGNER DISTRIBUTION

% % Reverberation Time
load('data/ReverberationTime/tau_RT_data_1-100CR_Pos4_CRC2-750.mat');
tau_RT = mean(tau_RT_end);

% % Weyl Formula:
cavity_area = nnz(abs(Ez) > 0) * delta^2;               % Cavity Surface Area [m^2]
% Number of Modes for a 2D cavity [#]
N_TE_2D = round(((2 * pi * cavity_area) / (c0^2)) * freq_car * bw_sig);

WD_freq_int = bw_sig / N_TE_2D;                      % Frequency interval according to # of Modes [Hz]
WD_freq_modes = (bw_flo + (WD_freq_int / 2)):WD_freq_int:(bw_fhi - (WD_freq_int / 2));

% % Wigner Distribution

% alpha_WD = (-1 / tau_RT);                               % Constant for WD [Hz]
% freq_WD = freq(ind_bw_flo:ind_bw_fhi) - freq_car;       % Frequency range for WD [Hz]
% W_alpha = WignerDist(alpha_WD, time, freq_WD);       % WD [*amplitude]

%% !!! WARNING: DO NOT USE, TOO COMPUTATIONALLY EXPENSIVE !!!
% % Frequency modes for sum of Wigner Distributions
% nb_elem_bw = ind_bw_fhi - ind_bw_flo + 1;
% 
% sum_method = "(+-)";
% 
% % Frequency range for sum of Wigner Distributions
% [xxx, ind_WD_sum_flo] = min(abs(freq - (WD_freq_modes(1) - (bw_sig / 2)) ));
% [xxx, ind_WD_sum_fhi] = min(abs(freq - (WD_freq_modes(end) + (bw_sig / 2)) ));
% WD_sum_freq = freq(ind_WD_sum_flo:ind_WD_sum_fhi);
% 
% % Allocating space for the sum of Wigner Distributions
% WD_alpha_sum = zeros(length(time), length(WD_sum_freq));
% 
% ind_WD_sum_freq = 1;
% 
% for ind_modes = 1:(N_TE_2D - 1)
%     [xxx, ind_WD_sum_freq] = min(abs(WD_sum_freq - ...
%         (WD_freq_modes(ind_modes) - (bw_sig / 2)) ))
%     
%     for elem_time = 1:length(time)
%         for elem_freq = 1:nb_elem_bw
%             if sum_method == "(+-)"
%                 WD_alpha_sum(elem_time, elem_freq + ind_WD_sum_freq) = ...
%                 WD_alpha_sum(elem_time, elem_freq + ind_WD_sum_freq) + ...
%                 W_alpha(elem_time, elem_freq);
%             elseif sum_method == "(abs)"
%                 WD_alpha_sum(elem_time, elem_freq + ind_WD_sum_freq) = ...
%                 WD_alpha_sum(elem_time, elem_freq + ind_WD_sum_freq) + ...
%                 abs(W_alpha(elem_time, elem_freq));
%             end
%         end
%     end
% end
% 
% new_dir = "images/WignerDistribution/cavity";
% mkdir(new_dir);
% 
% close all;
% nb_fig = 0;
% 
% nb_fig = nb_fig + 1;
% figure(nb_fig)
% figure('visible','off');
% 
% [surf_Y, surf_X] = meshgrid(WD_sum_freq, time);
% surf_Z = abs(WD_alpha_sum);
% 
% surf(surf_X, surf_Y, surf_Z, 'edgecolor', 'none')
% colormap jet;
% colorbar;
% 
% title('Wigner Distribution of Room Impulse Response', 'interpreter', 'latex')
% xlabel('time', 'interpreter', 'latex')
% ylabel('Frequency', 'interpreter', 'latex')
% 
% tlim_t0 = 0;
% tlim_t1 = tau_RT;
% tlim = [tlim_t0, tlim_t1];
% 
% flim_f0 = bw_flo;
% flim_f1 = bw_fhi;
% flim = [flim_f0, flim_f1];
% 
% xlim([tlim(1), tlim(2)]);
% ylim([flim(1), flim(2)]);
% 
% view(0, 90)
% 
% figure_title = new_dir + "/sum_W_alpha-" + num2str(N_TE_2D) + "_modes" + ...
%     sum_method + "_omega.png";
% saveas(gcf, figure_title)

%% TEST 1

alpha = -5;
fc_wvd = 60;
nb_elem = 5000;
time_wvd = linspace(0, 1, nb_elem);
dt_wvd = time_wvd(2) - time_wvd(1);
fs_wvd = 1 / dt_wvd; 
df_wvd = 1 / time_wvd(end);
freq_wvd = (0:length(time_wvd) - 1) * df_wvd;

wvd_exp = exp(2 * alpha * time_wvd) .* (sin(2 * pi * fc_wvd * time_wvd) / (pi * 2 * pi * fc_wvd));

signal = wvd_exp;
FFT_norm = abs(fft(signal)) ./ nnz(signal);

tau_wvds = [0.2, 0.4, 1];

for ind_tau_wvd = 1:length(tau_wvds)
    tau_wvd = tau_wvds(ind_tau_wvd);
    ind_twin_wvd = find(min(abs(time_wvd - tau_wvd)) == abs(time_wvd - tau_wvd));
    if (rem(ind_twin_wvd, 2) == 0)
        ind_twin_wvd = ind_twin_wvd - 1;
    end
    if tau_wvd == time_wvd(end)
        twin_wvd = ones(1, ind_twin_wvd);
    else
        twin_wvd = hamming(ind_twin_wvd, 'periodic'); 
    end

    ind_fwin_wvd = length(signal);
    if (rem(ind_twin_wvd, 2) == 1)
        ind_fwin_wvd = ind_fwin_wvd - 1;
    end
    fwin_wvd = ones(1, ind_fwin_wvd);

    tlim = [time_wvd(1), time_wvd(end)];
    freq_zoom = 50;
    flim = [fc_wvd - freq_zoom, fc_wvd + freq_zoom];

    new_dir_WVD = "images/WignerDistribution/troubleshooting";
    title_label = ' test 1';
    n_fontsize = 14;

    ind_win_wvd = [1, ind_twin_wvd];

    WVD_plot(time_wvd, freq_wvd, fs_wvd,  ind_win_wvd, signal, FFT_norm, ...
        twin_wvd, fwin_wvd, tlim, flim, n_fontsize, title_label)

    figure_title = new_dir_WVD + "/WVD-test_1-t_win=" + num2str(tau_wvd) + "s.png";
    saveas(gcf, figure_title)
end


%% TEST 2

alpha = -5;
fc_wvd = 60;
nb_elem = 5000;
time_wvd = linspace(0, 1, nb_elem);
dt_wvd = time_wvd(2) - time_wvd(1);
fs_wvd = 1 / dt_wvd; 
df_wvd = 1 / time_wvd(end);
freq_wvd = (0:length(time_wvd) - 1) * df_wvd;

wvd_exp = exp(2 * alpha * time_wvd) .* (...
    (sin(2 * pi * fc_wvd * time_wvd) / (pi * 2 * pi * fc_wvd)) + ...
    (sin(2 * pi * 3 * fc_wvd * time_wvd) / (pi * 2 * pi * fc_wvd)) + ...
    (sin(2 * pi * 2 * fc_wvd * time_wvd) / (pi * 2 * pi * fc_wvd)));

signal = wvd_exp;
FFT_norm = abs(fft(signal)) ./ nnz(signal);

tau_wvds = [0.2, 0.4, 1];

for ind_tau_wvd = 1:length(tau_wvds)
    tau_wvd = tau_wvds(ind_tau_wvd);
    ind_twin_wvd = find(min(abs(time_wvd - tau_wvd)) == abs(time_wvd - tau_wvd));
    if (rem(ind_twin_wvd, 2) == 0)
        ind_twin_wvd = ind_twin_wvd - 1;
    end
    if tau_wvd == time_wvd(end)
        twin_wvd = ones(1, ind_twin_wvd);
    else
        twin_wvd = hamming(ind_twin_wvd, 'periodic'); 
    end

    ind_fwin_wvd = length(signal);
    if (rem(ind_twin_wvd, 2) == 1)
        ind_fwin_wvd = ind_fwin_wvd - 1;
    end
    fwin_wvd = ones(1, ind_fwin_wvd);

    tlim = [time_wvd(1), time_wvd(end)];
    freq_zoom = 50;
    flim = [fc_wvd - freq_zoom, 3 * fc_wvd + freq_zoom];

    new_dir_WVD = "images/WignerDistribution/troubleshooting";
    title_label = ' test 2';
    n_fontsize = 14;

    ind_win_wvd = [1, ind_twin_wvd];

    WVD_plot(time_wvd, freq_wvd, fs_wvd,  ind_win_wvd, signal, FFT_norm, ...
        twin_wvd, fwin_wvd, tlim, flim, n_fontsize, title_label)

    figure_title = new_dir_WVD + "/WVD-test_2-t_win=" + num2str(tau_wvd) + "s.png";
    saveas(gcf, figure_title)
end

%%

alpha = (-1 / tau_RT);
h_ana_t = zeros(size(time));

for ind_N_TE_2D = 1:N_TE_2D
    h_ana_t = h_ana_t + ...
        exp(2 * alpha * time) .* ...
        (sin(2 * pi * WD_freq_modes(ind_N_TE_2D) * time) ./ (pi * 2 * pi * WD_freq_modes(ind_N_TE_2D))); 
end

close all

ax(1) = subplot(1, 2, 1);
plot(time, h_ana_t)

ax(2) = subplot(1, 2, 2);
plot(freq, abs(fft(h_ana_t)))


%% ------------------------------------------------------------------------
%                               FUNCTIONS
%  ------------------------------------------------------------------------

% All relevant functions.

% Wigner Distribution Function with alpha
function [W_alpha] = WDF_alpha(alpha, time, freq)
    W_alpha = zeros(length(time));

%     for ind_t = 1:length(time)
%         for ind_f = 1:length(time)
%             W_alpha(ind_f, ind_t) = (exp(2 * alpha * time(ind_t)) * ...
%                 (sin(2 * freq(ind_f) * time(ind_t)) / (pi * freq(ind_f)))); 
%         end
%     end

    % Using omega = 2 * pi * freq:
    % No difference when not using omega.
    for ind_t = 1:length(time)
        for ind_f = 1:length(freq)
            W_alpha(ind_f, ind_t) = (exp(2 * alpha * time(ind_t)) * ...
                (sin(2 * 2 * pi * freq(ind_f) * time(ind_t)) / (pi * 2 * pi * freq(ind_f)))); 
        end
    end
end

