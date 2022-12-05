
%% ------------------------------------------------------------------------
%                             DATA SMOOTHING
%  ------------------------------------------------------------------------

% Studying the effects of how many receiver antennas are needed to have an
% accurate enough representation of the original emitted signal in spite of
% the reverberations of the chamber.

clear all; close all; clc; format long;
addpath("functions");
addpath("data");

rng('default');


%% INPUT PARAMETERS

% Impulse parameters (For TFA)

signal_type = "chirp_wavelet";

N_recs = 500;                                                              % Number of receiver antennas [#]

new_dir_0 = "images/Statistics/" + signal_type;
mkdir(new_dir_0);


%% ------------------------------------------------------------------------
%                         SIMULATION (VALENTIN'S CODE)
%  ------------------------------------------------------------------------

% Definition of simulation of a Reverberant Chamber.

if signal_type == "wavelet"
    d_impulse = 8e-9;                                                      % Duration of the impulse [s]
    f0 = 2.4e9;                                                            % Carrier frequency of the impulse [Hz]
elseif signal_type == "chirp_wavelet" || signal_type == "chirp_sinusoidal"
%     d_impulse = 200e-9;                                                  % (Chirp) Duration of the impulse [s]
    d_impulse = 50e-9;                                                     % (Chirp) Duration of the impulse [s]
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
f_impulse = 1 / d_impulse;                     
lambda_min = c0 /(f0 + 2 * f_impulse);
fact_delta = 10;
delta = lambda_min / fact_delta;
dt= delta / c0 / sqrt(2);
fs = 1 / dt;                                                               % Sampling Frequency [Hz]


%% POSITIONING OF EMITTER(S) AND RECEIVER(S) (RANDOM)

corr_delta = lambda_min / 2;                                               % Distance at which correlation changes [m]
position_receiver_all = zeros(2, N_recs + 1);                              % Location of receivers [m]

pos_emt = [0.5375; 0.5375];                                                % Alex's location of the emitter [m]

low_X = 0;                                                                 % Lower boundary for random X point [m]
upp_X = 2.2;                                                               % Upper boundary for random X point [m]

low_Y = corr_delta;                                                        % Lower boundary for random Y point [m]
upp_Y = 1.8;                                                               % Upper boundary for random Y point [m]

position_receiver_all(1, 1) = (upp_X-low_X).*rand(1) + low_X; 
position_receiver_all(2, 1) = (upp_Y-low_Y).*rand(1) + low_Y;

% For conditioners
cavity_center = [1.125, 0.75];
cavity_radius = 1.125 - corr_delta;

emit_center = [pos_emt(1), pos_emt(2)];
emit_radius = corr_delta * 10;

i = 2;

while (i < N_recs + 2)
    rand_X = (upp_X-low_X).*rand(1) + low_X;                               % Random X point [m]
    rand_Y = (upp_Y-low_Y).*rand(1) + low_Y;                               % Random Y point [m]
    
    incorrect = false;
    for j = 1:(i - 1) 
        % 1. Is there a corr_delta distance between random points?
        dist_corr_delta = sqrt((rand_X - position_receiver_all(1, j))^2 + ...
            (rand_Y - position_receiver_all(2, j))^2);
        if (dist_corr_delta <= corr_delta)
            incorrect = true; 
            break;
        end
        
        % 2. Is the random point in the cavity?
        dist_cavity = ((rand_X - cavity_center(1))^2 + ...
            (rand_Y - cavity_center(2))^2) > (cavity_radius^2);
        if dist_cavity
            incorrect = true; 
            break;
        end
        
%         % 3. Is the random point far away enough from the emitter?
%         dist_emit = ((rand_X - emit_center(1))^2 + ...
%             (rand_Y - emit_center(2))^2) < (emit_radius^2);
%         if dist_emit
%             incorrect = true; 
%             break;
%         end
    end
    
    if incorrect == false
        position_receiver_all(1, i) = rand_X; 
        position_receiver_all(2, i) = rand_Y;
        
        i = i + 1;
    end
end

pos_rec = position_receiver_all(:, 2:(N_recs + 1));

close all

% To check positioning
openfig('data/CavitySpace/cavity_space_empty.fig');
figure(1);
hold on;
viscircles(cavity_center, cavity_radius, 'Color', 'red');
% viscircles(emit_center, emit_radius, 'Color', 'green');
plot(pos_rec(1, :), pos_rec(2, :), 'rx')
plot(pos_emt(1), pos_emt(2), 'g*')
hold off;

figure_title = new_dir_0 + "/Cavity-" + num2str(N_recs) + "Recs.png";
saveas(gcf, figure_title)


%% RECEIVED SIGNAL

% Display of the impulse respone
Display_rep_imp = 0;                                                       % If=1: display, if=0: no display
[r_t, Ez] = function_impulse_response(param_simu, pos_emt, pos_rec, Display_rep_imp, time, e_t);
% [r_t_est, Ez] = function_impulse_response(param_simu, pos_pos, pos_rec, Display_rep_imp, time, e_t);

% Frequency Range
df = 1 / time(end);
freq = [0:length(time) - 1] .* df;


%% FOURIER TRANSFORM

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

% Bandwidth index of the Signal
[xxx, ind_bw_flo] = min(abs(freq - bw_flo));
[xxx, ind_bw_fhi] = min(abs(freq - bw_fhi));

%% Room Impulse caused by Reverberation in VC
H_f = zeros(size(r_t));
% E_f_recon = zeros(size(r_t));

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

%     % Reconstructing E_f from H_f
%     E_f_recon(ind0_H_f:ind1_H_f, r_t_i) = (R_f(ind0_H_f:ind1_H_f, r_t_i) ./ H_f(ind0_H_f:ind1_H_f, r_t_i));
%     E_f_recon(ind3_H_f:-1:ind2_H_f, r_t_i) = E_f_recon(ind0_H_f:ind1_H_f, r_t_i);
end

% Symmetric in the IFFT, treats H(f) as if it were conjugate symmetric by 
% ignoring the second half of its elements (that are in the negative 
% frequency spectrum).
h_t = real(ifft(H_f));                                                     % H(f) in time          

% e_t_recon = real(ifft(E_f_recon));
% e_t_recon = flip(e_t_recon, 1);


%%

close all;

H_f_mean = mean(abs(H_f), 2);

nb_fig = 0;

nb_fig = nb_fig + 1;
figure(nb_fig);

% set(gcf, 'Position', get(0, 'Screensize'));

hold on;
grid on;
plot(freq .* 1e-9, abs(H_f));
plot(freq .* 1e-9, abs(H_f_mean), 'k*-');
hold off;

xlim([bw_flo, bw_fhi] .* 1e-9);

xlabel("Frequency [GHz]");
ylabel("FFT Amplitude []");
title("Mean of " + num2str(N_recs) + " abs(H(f))");

figure_title = new_dir_0 + "/H_f-mean_abs.png";
saveas(gcf, figure_title)

H_f_mean = mean(H_f, 2);

nb_fig = nb_fig + 1;
figure(nb_fig);

% set(gcf, 'Position', get(0, 'Screensize'));

hold on;
grid on;
plot(freq .* 1e-9, abs(H_f));
plot(freq .* 1e-9, abs(H_f_mean), 'k*-');
hold off;

xlim([bw_flo, bw_fhi] .* 1e-9);

xlabel("Frequency [GHz]");
ylabel("FFT Amplitude []");
title("Mean of " + num2str(N_recs) + " H(f) (complex)");

figure_title = new_dir_0 + "/H_f-mean_complex.png";
saveas(gcf, figure_title)


%% ------------------------------------------------------------------------
%                              STATISTICS
%  ------------------------------------------------------------------------

% Analyzing the required number of antennas for obtaining an accurate emmited 
% signal e(t) from an n number of received signal r(t).

%% CORRELATION STUDIES

clc

freqs_smt = [25, 50, 100, 200] .* 1e6;                                     % Smoothing frequencies [Hz]
% ----------[ 1,  2,   3,   4]                                             Indexes for smoothing frequencies 

% IMPORTANT:
% M -> Refers to the smoothing frequency
% N -> Refers to the number of receivers

PE_M_N = zeros(length(freqs_smt), N_recs);
corr_smt_M_N = zeros(length(freqs_smt), N_recs);
corr_alpha_M_N = zeros(length(freqs_smt), N_recs);

R_f_avg_M_N = zeros(length(freqs_smt), N_recs, length(time));
R_f_smt_avg_M_N = zeros(length(freqs_smt), N_recs, length(time));

for ind_freq_smt = 1:length(freqs_smt)
    freq_smt = freqs_smt(ind_freq_smt);                                    % Current smoothing frequencies [Hz]
    [xxx, ind_smoothing] = min(abs(freq - freq_smt));                      % Index for smoothing frequency [#]

    for ind_N_recs = 1:N_recs
        % METHOD: Mean -> Smooth
        % Mean from N abs(R(f)):
        R_f_avg = mean(R_f_norm(:, 1:ind_N_recs), 2);   
        R_f_avg_M_N(ind_freq_smt, ind_N_recs, :) = R_f_avg;
        % Data smoothing using 'movmean':
        R_f_smt_avg = smoothdata(R_f_avg, 'movmean', ind_smoothing);   
        R_f_smt_avg_M_N(ind_freq_smt, ind_N_recs, :) = R_f_smt_avg;

        % E(f) normalized to have a maximum of 1:
        E_f_norm_1 = normalize(E_f_norm(:, 1), 'norm', Inf);
        % R(f)_{smoothed} normalized to have a maximum of 1:
        R_f_smt_avg_norm_1 = normalize(R_f_smt_avg(:, 1), 'norm', Inf);

        % Calculating the Percentage Error between the normalized values to 1 of
        % smooth(mean(R_f)) and E(f) in the bandwidth of the signal:
        approx = R_f_smt_avg_norm_1(ind_bw_flo:ind_bw_fhi);
        exact = E_f_norm_1(ind_bw_flo:ind_bw_fhi);
        per_error = (abs(approx - exact) ./ exact) * 100;
        PE_rec = mean(per_error);

        PE_M_N(ind_freq_smt, ind_N_recs) = PE_rec;
        
        corr_A = E_f_norm;                                                 % Data A for Correlation []
        corr_B = R_f_smt_avg;                                              % Data B for Correlation []
        corr_elem = 100;                                                   % Number of elements for correlation line fitting [#] 

        corr_smt = corrcoef(corr_A, corr_B);                                   % Correlation coefficient between A and B []
        corr_smt = max(abs(xcorr(corr_A, corr_B,'normalized')));                                   % Correlation coefficient between A and B []


        corr_smt_M_N(ind_freq_smt, ind_N_recs) = corr_smt;
        
        % corr_smt_eq = sum(((corr_A - mean(corr_A)) / std(corr_A)) .* ...
        %     ((corr_B - mean(corr_B)) / std(corr_B))) .* ...
        %     ((1 / (length(corr_A) - 1)));

        % Correlation line according to MATLAB function:
        corr_line_Y_1 = linspace(min(corr_B), max(corr_B), corr_elem) .* corr_smt;
        % X-axis correlation line:
        corr_line_X = linspace(min(corr_A), max(corr_A), corr_elem);
        % Fitted correlation line according to data:
        corr_polyfit = polyfit(corr_A, corr_B, 1);
        corr_linefit_Y = polyval(corr_polyfit, corr_line_X);

        % Checking if whether A or B should in the X or Y axes:
        if corr_linefit_Y(end) > corr_line_Y_1(end)
            corr_switch = 1;
            corr_A = R_f_smt_avg;
            corr_A_label = "smt(<|R(f)|>)";
            corr_B = E_f_norm;
            corr_B_label = "|E(f)|";

            corr_line_Y_1 = linspace(min(corr_B), max(corr_B), corr_elem) .* corr_smt;
            corr_line_X = linspace(min(corr_A), max(corr_A), corr_elem);
            corr_polyfit = polyfit(corr_A, corr_B, 1);
            corr_linefit_Y = polyval(corr_polyfit, corr_line_X);
        else
            corr_switch = 0;
            corr_A_label = "|E(f)|";
            corr_B_label = "smt(<|R(f)|>)";
        end

        % Percentage error between correlation from MATLAB function and fitted line [%]:
        % Percentage error between smoothing and E(f)_{norm_1}:
        approx = corr_line_Y_1;
        exact = corr_linefit_Y;
        per_error = (abs(approx - exact) ./ exact) * 100;
        corr_PE = mean(per_error);

        % For calculating the correct correlation coefficient:
        corr_acc = 5;
        corr_alpha = 1;
        corr_line_Y_alpha = linspace(min(corr_B), max(corr_B), corr_elem) .* corr_alpha;

        while (corr_PE > corr_acc)  
            corr_alpha = corr_alpha - 0.0001;

            if corr_switch == 0
                corr_line_Y_alpha = linspace(min(corr_B), max(corr_B), corr_elem) .* corr_alpha;
            elseif corr_switch == 1
                corr_line_Y_alpha = linspace(min(corr_A), max(corr_A), corr_elem) .* corr_alpha;
            end
            
            approx = corr_line_Y_alpha;
            exact = corr_linefit_Y;
            per_error = (abs(approx - exact) ./ exact) * 100;
            corr_PE = mean(per_error);
        end

        corr_alpha_M_N(ind_freq_smt, ind_N_recs) = corr_alpha;
    end
end


%% PLOTTING

% For an ever increasing number of plots
n_fontsize = 14;
n_markersize = 18;

% ind_fsmt = 2;
inds_rec = [1, 2, 4, 5, 10, 100, 500];

for ind_fsmt = 1:length(freqs_smt)
    new_dir_1 = new_dir_0 + "/f_smt=" + num2str(freqs_smt(ind_fsmt) * 1e-6) + "MHz";
    mkdir(new_dir_1)
    
    for ind_plot = 1:length(inds_rec)
        close all
        
        ind_rec = inds_rec(ind_plot);
        % To check positioning
        openfig('data/CavitySpace/cavity_space_empty.fig');
        cavity_fig = gca;
        cavity_fig_plot = get(cavity_fig, 'children');

        nb_fig = 1; 
        nb_fig = nb_fig + 1;
        figure(nb_fig);
        set(gcf, 'Position', get(0, 'Screensize'));

        ax_1 = subplot(6, 6, [13, 14, 19, 20, 25, 26]);

        copyobj(cavity_fig_plot, ax_1);

        hold on;
        grid on;
        viscircles(cavity_center, cavity_radius, 'Color', 'red');
        plot(pos_rec(1, 1:end), pos_rec(2, 1:end), 'rx')
        plot(pos_emt(1), pos_emt(2), 'g*', 'MarkerSize', n_markersize)
        hold off;

        xlim([0, Lx + 0.4]);
        ylim([0, Ly]);

        % Labelling
        xlabel('x-axis [m]');
        ylabel('y-axis [m]');
        title({'Cavity', "(1 Emitter, " + num2str(N_recs) + " Receivers)"});

        ax_2 = subplot(6, 6, [3, 4, 5, 6, 9, 10, 11, 12]);

        hold on;
        grid on;
        plot(1:N_recs, PE_M_N(ind_fsmt, :))
        hold off;

        % Labelling
        xlabel('Number of Receivers [#]');
        ylabel('Percentage Error [%]');
        title({'Percentage Error of smt(<|R_{1-N}(f)|>) and |E(f)|', ...
            "using movmean and f_{smt}=" + num2str(freqs_smt(ind_fsmt) * 1e-6) + " MHz"});

        ax_3 = subplot(6, 6, [21, 22, 23, 24, 27, 28, 29, 30]);

        hold on;
        grid on;
        plot(1:N_recs, corr_alpha_M_N(ind_fsmt, :))
        hold off;

        % Labelling
        xlabel('Number of Receivers [#]');
        ylabel('Correlation [0, 1]');
        title({"\rho_{\alpha}(" + corr_A_label + "," + corr_B_label + ")", ...
            "using movmean and f_{smt}=" + num2str(freqs_smt(ind_fsmt) * 1e-6) + " MHz"});

        figure_title = new_dir_1 + "/PE_corr.png";
        saveas(gcf, figure_title)

        nb_fig = nb_fig + 1;
        figure(nb_fig);
        set(gcf, 'Position', get(0, 'Screensize'));

        ax_1 = subplot(4, 6, [7, 8, 13, 14]);

        copyobj(cavity_fig_plot, ax_1);

        hold on;
        grid on;
        viscircles(cavity_center, cavity_radius, 'Color', 'red');
        plot(pos_rec(1, 1:ind_rec), pos_rec(2, 1:ind_rec), 'rx', 'MarkerSize', n_markersize)
        plot(pos_emt(1), pos_emt(2), 'g*', 'MarkerSize', n_markersize)
        hold off;

        xlim([0, Lx + 0.4]);
        ylim([0, Ly]);

        % Labelling
        xlabel('x-axis [m]');
        ylabel('y-axis [m]');
        title({'Cavity', "(1 Emitter, " + num2str(ind_rec) + " Receivers)"});

        ax_2 = subplot(4, 6, [3, 4, 9, 10]);

        E_f_norm_1 = normalize(E_f_norm(:, 1), 'norm', Inf);
        R_f_smt_avg_norm_1 = normalize(R_f_smt_avg_M_N(ind_fsmt, ind_rec, :), 'norm', Inf);

        hold on;
        grid on;
        plot(freq .* 1e-9, E_f_norm_1)
        plot(freq .* 1e-9, squeeze(R_f_smt_avg_norm_1))
        hold off;

        xlim([bw_flo, bw_fhi] .* 1e-9);

        % Labelling
        xlabel('Frequency [GHz]');
        ylabel('Normalized to 1 FFT Amplitude [0, 1]');
        title({'Waveform Comparison between', ...
            "|E(f)| and smt(<|R_{1-" + num2str(ind_rec) + "}(f)|>) with f_{smt}=" + ...
            num2str(freqs_smt(ind_fsmt) * 1e-6) + " MHz"});

        ax_3 = subplot(4, 6, [15, 16, 21, 22]);

        hold on;
        grid on;
        plot(freq .* 1e-9, R_f_norm(:, 1), 'r--')
        plot(freq .* 1e-9, squeeze(R_f_avg_M_N(ind_fsmt, ind_rec, :)), 'g*-')
        plot(freq .* 1e-9, squeeze(R_f_smt_avg_M_N(ind_fsmt, ind_rec, :)), 'b+-')
        hold off;

        xlim([bw_flo, bw_fhi] .* 1e-9);

        % Labelling
        xlabel('Frequency [GHz]');
        ylabel('FFT Amplitude []');

        str_1 = '|R_{1}(f)|';
        str_2 = "<|R_{1-" + num2str(ind_rec) + "}(f)|>";
        str_3 = 'smoothed';
        legend(str_1, str_2, str_3, 'location', 'best')

        ax_4 = subplot(4, 6, [11, 12, 17, 18]);

        corr_A = E_f_norm;                                                       
        corr_B = squeeze(R_f_smt_avg_M_N(ind_fsmt, ind_rec, :));                                        

        % Correlation line according to MATLAB function:
        corr_line_Y_1 = linspace(min(corr_B), max(corr_B), corr_elem) .* corr_smt_M_N(ind_fsmt, ind_rec);
        % X-axis correlation line:
        corr_line_X = linspace(min(corr_A), max(corr_A), corr_elem);
        % Fitted correlation line according to data:
        corr_polyfit = polyfit(corr_A, corr_B, 1);
        corr_linefit_Y = polyval(corr_polyfit, corr_line_X);

        % Checking if whether A or B should in the X or Y axes:
        if corr_linefit_Y(end) > corr_line_Y_1(end)
            corr_switch = 1;
            corr_A = R_f_smt_avg;
            corr_A_label = "smt(<|R(f)|>)";
            corr_B = E_f_norm;
            corr_B_label = "|E(f)|";

            corr_line_Y_1 = linspace(min(corr_B), max(corr_B), corr_elem) .* corr_smt_M_N(ind_fsmt, ind_rec);
            corr_line_X = linspace(min(corr_A), max(corr_A), corr_elem);
            corr_polyfit = polyfit(corr_A, corr_B, 1);
            corr_linefit_Y = polyval(corr_polyfit, corr_line_X);
        else
            corr_switch = 0;
            corr_A_label = "|E(f)|";
            corr_B_label = "smt(<|R(f)|>)";
        end

        corr_line_Y_alpha = linspace(min(corr_B), max(corr_B), corr_elem) .* corr_alpha_M_N(ind_fsmt, ind_rec);

        hold on;
        grid on;
        plot(corr_A, corr_B, 'b*')
        plot(corr_line_X, corr_line_Y_1, 'r')
        plot(corr_line_X, corr_line_Y_alpha, 'go')
        plot(corr_line_X, corr_linefit_Y, 'g--')
        hold off;

        % Labelling
        xlabel(corr_A_label);
        ylabel(corr_B_label);
        title({"\rho(" + corr_A_label + "," + corr_B_label + ") = " + num2str(corr_smt_M_N(ind_fsmt, ind_rec)), ...
            "\rho_{\alpha}(" + corr_A_label + "," + corr_B_label + ") = " + num2str(corr_alpha_M_N(ind_fsmt, ind_rec))})

        figure_title = new_dir_1 + "/N_recs=" + num2str(ind_rec) + ".png";
        saveas(gcf, figure_title)
    end
end


%% ACCURACY VS. # OF RECEIVED SIGNALS 

% freqs_smt = [25, 50, 100, 200] .* 1e6;
% freqs_smt = [freqs_smt, round(bw_sig, -6)];
% 
% for ind_fs = 1:length(freqs_smt)
%     freq_smt = freqs_smt(ind_fs);
%     % Using x MHz for smoothing window
%     [xxx, ind_smoothing] = min(abs(freq - freq_smt));
% 
%     new_dir_1 = new_dir_0 + "/smoothing-f_smt=" + num2str(freq_smt .* 1e-6) + "MHz";
%     mkdir(new_dir_1);
% 
%     % smt_method_cell = {'movmean', 'movmedian', 'gaussian', 'lowess', 'loess', ...
%     %     'rlowess', 'rloess', 'sgolay'};
%     smt_method_cell = {'movmean', 'movmedian', 'gaussian', 'lowess', 'loess', 'sgolay'};
%     n_smt_methods = size(smt_method_cell, 2);
% 
%     % Percentage error averages according to smoothing
%     per_error_smt_avg = zeros(N_recs, 1);
%     % Correlation between E(f) and smooth(mean(R(f)))
%     corr_smt_avg = zeros(N_recs, 1);
% 
%     % Mean percentage error averages from smoothed methods
%     mean_per_error_smt_avg_method = zeros(n_smt_methods, 1);
%     % Percentage error averages according to smoothed methods
%     per_error_smt_avg_method = zeros(N_recs, n_smt_methods);
% 
%     % Mean correlation from smoothed methods
%     mean_corr_smt_avg_method = zeros(n_smt_methods, 1);
%     % Correlation according to smoothed methods
%     corr_smt_avg_method = zeros(N_recs, n_smt_methods);
% 
%     for ind_smt = 1:n_smt_methods
%         smt_method = smt_method_cell{ind_smt};
% 
%         for ind_rec = 1:N_recs
%     %         % METHOD 1: Smooth -> Mean
%     %         % Smoothed R(f) according to the window made up of ind_smoothing
%     %         R_f_smt = smoothdata(R_f_norm(:, 1:ind_rec), smt_method, ind_smoothing);
%     %         % Averaging the R(f) smoothed
%     %         R_f_smt_avg = sum(R_f_smt, 2) / size(R_f_smt, 2);
% 
%             % METHOD 2: Mean -> Smooth
%             % Average from the magnitude: |R(f)|s
%             R_f_avg = mean(R_f_norm(:, 1:ind_rec), 2);
%             R_f_smt_avg = smoothdata(R_f_avg, smt_method, ind_smoothing);
% 
%             % E(f) normalized to have a maximum of 1:
%             E_f_norm_1 = normalize(E_f_norm(:, 1), 'norm', Inf);
%             % R(f)_{smoothed} normalized to have a maximum of 1:
%             R_f_smt_avg_norm = normalize(R_f_smt_avg(:, 1), 'norm', Inf);
%             % Percentage error between smoothing and E(f)_{norm_1}:
%             approx = R_f_smt_avg_norm(ind_bw_flo:ind_bw_fhi);
%             exact = E_f_norm_1(ind_bw_flo:ind_bw_fhi);
%             per_error_smt = (abs(approx - exact) ./ exact);
%             per_error_smt_avg(ind_rec) = mean(per_error_smt);
% 
%             % Correlation
%             corr_rho = corrcoef(E_f_norm(ind_bw_flo:ind_bw_fhi), R_f_smt_avg(ind_bw_flo:ind_bw_fhi));
%             corr_smt_avg(ind_rec) = corr_rho(1, 2);
%         end
% 
%         mean_per_error_smt_avg_method(ind_smt) = mean(per_error_smt_avg);
%         per_error_smt_avg_method(:, ind_smt) = per_error_smt_avg;
% 
%         mean_corr_smt_avg_method(ind_smt) = mean(corr_smt_avg);
%         corr_smt_avg_method(:, ind_smt) = corr_smt_avg;
% 
%         % % Plotting 
%         close all;
% 
%         n_fig = n_fig + 1;
%         figure(n_fig);
% 
%         set(gcf, 'Position', get(0, 'Screensize'));
%         ax(1) = subplot(2, 2, [1, 3]);
%         hold on;
%         grid on;
%         plot(freq .* 1e-9, E_f_norm_1);
%         plot(freq .* 1e-9, R_f_smt_avg_norm);
%         hold off;
% 
%         xlim([bw_flo, bw_fhi] .* 1e-9);
% 
%         xlabel('Frequency [GHz]');
%         ylabel('Normalized Amplitude [0, 1] ');
%         title({"Smoothing of " + num2str(N_recs) + " R(f)s using " + smt_method, ...
%             "with a smoothing window of " +  num2str(freq_smt * 1e-6) + " MHz"});
%         legend('E(f)_{norm}', 'R(f)_{smt avg, norm}', 'location', 'best')
% 
%         set(gca,'fontsize',n_fontsize);
% 
%         ax(2) = subplot(2, 2, 2);
%         hold on;
%         grid on;
%         plot(1:N_recs, per_error_smt_avg);
%         hold off;
% 
%         xlabel('Number of Receivers [#]');
%         ylabel('Percentage Error [%] ');
%         title({"Percentage Error (%) of " + num2str(N_recs) + " R(f)s using " + smt_method, ...
%             "with a smoothing window of " +  num2str(freq_smt * 1e-6) + " MHz"});
% 
%         set(gca,'fontsize',n_fontsize);
% 
%         ax(3) = subplot(2, 2, 4);
%         hold on;
%         grid on;
%         plot(1:N_recs, corr_smt_avg);
%         hold off;
% 
%         xlabel('Number of Receivers [#]');
%         ylabel('Correlation [-1,1] ');
%         title({"Correlation of " + num2str(N_recs) + " R(f)s using " + smt_method, ...
%             "with a smoothing window of " +  num2str(freq_smt * 1e-6) + " MHz"});
% 
%         set(gca,'fontsize',n_fontsize);
% 
%         figure_title = new_dir_1 + "/" + num2str(ind_smt) + "-" + smt_method + ".png";
%         saveas(gcf, figure_title)
%     end
% 
%     n_fig = n_fig + 1;
%     figure(n_fig);
% 
%     set(gcf, 'Position', get(0, 'Screensize'));
%     ax(1) = subplot(1, 2, 1);
%     hold on;
%     grid on;
%     plot(mean_per_error_smt_avg_method, 'k*-');
%     hold off;
% 
%     xticks([1:n_smt_methods])
%     xticklabels(smt_method_cell);
%     xlabel("Smoothing Method");
%     ylabel("Mean of Percentage Error Averages [%]");
%     title({"Mean of PE of Different Methods", ...
%         "with a smoothing window of " +  num2str(freq_smt * 1e-6) + " MHz"});
% 
%     ax(2) = subplot(1, 2, 2);
%     hold on;
%     grid on;
%     plot(mean_corr_smt_avg_method, 'k+-');
%     hold off;
% 
%     xticks([1:n_smt_methods])
%     xticklabels(smt_method_cell);
%     xlabel("Smoothing Method");
%     ylabel("Mean of Correlations [-1,1]");
%     title({"Mean of Correlations of Different Methods", ...
%         "with a smoothing window of " +  num2str(freq_smt * 1e-6) + " MHz"});
% 
%     figure_title = new_dir_1 + "/method-per_error_mean.png";
%     saveas(gcf, figure_title)
% 
%     n_fig = n_fig + 1;
%     figure(n_fig);
% 
%     set(gcf, 'Position', get(0, 'Screensize'));
%     ax(1) = subplot(1, 2, 1);
%     hold on;
%     grid on;
%     plot(per_error_smt_avg_method);
%     hold off;
% 
%     % xticklabels(smt_method_cell);
%     xlabel('Number of Receivers [#]');
%     ylabel("Percentage Error of Averages [%]");
%     title({"Analysis (PE) of Different Methods", ...
%         "with a smoothing window of " +  num2str(freq_smt * 1e-6) + " MHz"});
%     legend(smt_method_cell, 'location', 'best');
% 
%     ax(2) = subplot(1, 2, 2);
%     hold on;
%     grid on;
%     plot(corr_smt_avg_method);
%     hold off;
% 
%     % xticklabels(smt_method_cell);
%     xlabel('Number of Receivers [#]');
%     ylabel("Correlation [-1,1]");
%     title({"Analysis (corr) of Different Methods", ...
%         "with a smoothing window of " +  num2str(freq_smt * 1e-6) + " MHz"});
%     legend(smt_method_cell, 'location', 'best');
% 
%     figure_title = new_dir_1 + "/method-per_error_method.png";
%     saveas(gcf, figure_title)
% end


%% --------------------------- FUNCTIONS ----------------------------------

function [x1, y1] = coord_circle(xc, yc, r)
    % xc: x-axis center of circle
    % yc: y-axis center of circle
    % r: radius
    
    ang = 0:0.01:2*pi; 
    x1 = r * cos(ang) + xc;
    y1 = r * sin(ang) + yc;
end

