
clear all; close all; clc;
addpath("functions")
 
%% ------------------------------------------------------------------------
%                               FOURIER TRANSFORM
%  ------------------------------------------------------------------------

% Comparing the analytical solution of the Fourier Transform with the 
% results of the MATLAB fft function.

%% Initial Parameters
N_samp = 1000;                                                             % Number of useful samples
fs = 1000;                                                                 % Sampling frequency [Hz]
dt = 1 / fs;                                                               % Sampling time [s]
time = (0:(N_samp - 1)) .* dt;                                             % Time domain [s]
N_samp = length(time);                                                     % Number of useful samples [#]
df = fs / N_samp;
freq = (0:(N_samp - 1)) .* df;                                             % Frquency range [Hz]

f0 = 25;                                                                   % Carrier frequency [hZ]
a = 2 * pi * f0;                                                           % Constant for signal
f_t0 = @(t) (1 - cos((a ./ 2) .* t)) .* (sin(a .* t));                     % Input Function Handle

%% Computation of Magnitude
F_w__Int_Re = []; F_w_Int_Im = []; Mg_F_w = [];
for zeta = 1:length(freq)
    F_w_Re = f_t0(time) .* cos(2 .* pi .* freq(zeta) .* time);
    F_w_Int_Re(zeta) = trapz(time, F_w_Re);

    F_w_Im = f_t0(time) .* sin(2 .* pi .* freq(zeta) .* time);
    F_w_Int_Im(zeta) = trapz(time, F_w_Im);

    Mg_F_w(zeta) = sqrt(F_w_Int_Re(zeta).^2 + F_w_Int_Im(zeta).^2); 
end

F_w_fft = fft(f_t0(time));
F_w_fft_scaled = F_w_fft ./ N_samp;                                        % FFT of Input Function

% Computational Inverse Fourier Transform
f_t_Mg_ifft = zeros(1, length(freq));
for zeta = 1:length(freq)
    f_t_Mg_ifft = f_t_Mg_ifft + ...
        (F_w_Int_Re(zeta) .* cos(2 .* pi .* freq(zeta) .* time)) + ...
        (F_w_Int_Im(zeta) .* sin(2 .* pi .* freq(zeta) .* time));
end
f_t_Mg_ifft = f_t_Mg_ifft ./ time(end);

f_t_ifft = ifft(F_w_fft);


%% ------------------------------------------------------------------------
%                         SHORT TIME FOURIER TRANSFORM
%  ------------------------------------------------------------------------ 

f_t_time = f_t0(time);                                                     % Original signal [...]

%% Varying Time Index

figure(1);
plot(time, f_t_time);
n_fontsize = 14;
set(gca,'fontsize',n_fontsize);
title('Original Signal', 'interpreter', 'latex'); 
xlabel('Time [s] ', 'interpreter', 'latex'); 
ylabel('Amplitude [] ', 'interpreter', 'latex');

tau_stft = 0.1:0.1:1;
ind_stft = zeros(length(tau_stft), 1);
time_stft = {};
freq_stft = {};
f_t_stft = {};                                                             % Segmented signals [...]
F_w_stft = {};                                                             % FFT segmented signals [...]

for i = 1:length(tau_stft)
    ind_stft(i) = find(min(abs(time - tau_stft(i))) == abs(time - tau_stft(i)));

    f_t_stft{i} = [f_t_time(1:(ind_stft(i))), ...
        zeros(1, (length(f_t_time) - ind_stft(i)))];
    F_w_stft{i} = fft(f_t_stft{i});
    F_w_stft{i} = (F_w_stft{i} ./ ind_stft(i)).^2;
    
    figure(2);
    hold on;
    grid on;
%     plot3(ones(ind_stft(i) - 1, 1) .* tau_stft(i), freq_stft{i}, abs(F_w_stft{i}));
    plot3(ones(length(f_t_time), 1) .* tau_stft(i), freq, abs(F_w_stft{i}));
end

view(3);
n_fontsize = 14;
set(gca,'fontsize',n_fontsize);
title('Fourier Transform with Varying Time Index $\tau$', 'interpreter', 'latex'); 
xlabel('Time Index $\tau$ [s] ', 'interpreter', 'latex'); 
ylabel('Frequency $\omega$ [Hz] ', 'interpreter', 'latex');
zlabel('Power of Magnitude $|F(\omega)|^2$', 'interpreter', 'latex');
ylim([0, (freq(end) / 2)]);

%% Short Time Fourier Transform Computations

tau_stft = 0.083;                                                          % STFT Time Index [s]
% tau_stft = 100e-9;
ind_stft = find(min(abs(time - tau_stft)) == abs(time - tau_stft));

time_stft = time((ind_stft + 1):end);                                      % Time for segmented signals [s]
freq_stft = linspace(freq(1), freq(end), ind_stft);
f_t_stft = zeros(length(time), (length(time) - ind_stft));                 % Segmented signals [...]
F_w_stft = zeros(length(time), (length(time) - ind_stft));                 % FFT segmented signals [...]

for i = 1:(length(time) - ind_stft)
    win = hann(ind_stft)';                                                 % Window used
    f_t_stft(1:ind_stft, i) = f_t_time(i:(ind_stft + i - 1)) .* win;
    
    F_w_stft(:, i) = fft(f_t_stft(:, i)) ./ ind_stft;
%     F_w_stft(:, i) = (F_w_stft(:, i) ./ ind_stft);
end

n_fig = 0;

% n_fig = n_fig + 1;
% figure(n_fig);
% stft(f_t_time, fs,'FrequencyRange','onesided','FFTLength',length(time));
% colormap jet;
% colorbar;

n_fig = n_fig + 1;
figure(n_fig);
ylim_1 = 0;
ylim_2 = (freq(end) / 2);
% set(gcf, 'Position', get(0, 'Screensize'));
ax(1) = subplot(3, 3, [1, 2]);
hold on;
grid on;
plot(time(1:ind_stft), win, 'r')
plot(time, f_t_time, 'b')
hold off;
% Labelling


ax(2) = subplot(3, 3, [6, 9]);
hold on;
grid on;
plot(abs(fft(f_t_time)), freq)
hold off;
ylim([ylim_1, ylim_2]);
% Labelling


% Time vs. Frequency
ax(3) = subplot(3, 3, [4, 8]);
[X_time_stft, Y_freq_stft] = meshgrid(time_stft, freq);
surf_stft = surf(X_time_stft, Y_freq_stft, abs(F_w_stft), 'FaceColor', 'interp');
view(0, 90)
set(surf_stft,'LineStyle','none');
n_fontsize = 10;
set(gca,'fontsize',n_fontsize);
% Colorbar definition
% cbar = colorbar('westoutside');
% cbar.Label.String = '\fontsize{12} []';
% cbar_limits = get(cbar,'Limits');
% set(cbar,'fontsize',n_fontsize,'Ticks',linspace(cbar_limits(1), cbar_limits(2), 10));
colormap jet;
% Labelling
xlabel('Time [s] ', 'interpreter', 'latex'); 
ylabel('Frequency $\omega$ [Hz] ', 'interpreter', 'latex');
% zlabel('Power of Magnitude $|F(\omega)|$', 'interpreter', 'latex');
ylim([ylim_1, ylim_2]);
grid on;

sgtitle('Short-Time Fourier Transform')

%% ------------------------------------------------------------------------
%                               VALENTIN's CODE
%  ------------------------------------------------------------------------

% PARAMETERS

% Impulse parameters
d_impulse=10e-9;                                                           % Duration of the impulse [s]
f0=2.4e9;                                                                  % Carrier frequency of the impulse [Hz]
d_before=0e-9;                                                             % delay before impulse [s]
d_after=500e-9;                                                            % delay after impulse [s]
A_impulse=1;                                                               % Amplitude of the impulse [V/m]

% Cavity parameters
Lx=1.9; Ly=1.9;                                                            % Cavity dimensions [m,m]
sigma=(6.5)*1e2;                                                           % Conductivity of the walls [S/m]
coef_losses_ok=1;                                                          % If=1: losses, if=0: PEC condition 
coef_regular_chao=1;                                                       % If=1: regular cavity, if=2: chaotic cavity 
% % WARNING: lossless chaotic cavity not coded (possible: lossless or with losses reglar cavity, or with losses chaotic cavity)
[time, e_t, param_simu] = function_parameters_simu(Lx,Ly,d_impulse,A_impulse,f0,d_before,d_after,coef_regular_chao,coef_losses_ok,sigma);

% IMPULSE RESPONSE

% Receiver and emitter locations
position_receiver=[0.6;0.2];                                               % Location of the receiver [m]
position_emitter=[1.2;1.3];                                                % Location of the emitter [m]

% Display of the impulse respone
Display_rep_imp=0;                                                         % If=1: display, if=0: no display
[r_t, Ez] = function_impulse_response(param_simu,position_emitter,position_receiver,Display_rep_imp,time,e_t);

%% ------------------------------------------------------------------------
%                        ANALYTICAL INITIAL IMPULSE
%  ------------------------------------------------------------------------

% For time duration
c0=3e8;                                                                    % Speed of light [m/s]
f_impulse=1/d_impulse;
lambda_min=c0/(f0+2*f_impulse);
fact_delta=10;
delta=lambda_min/fact_delta;

dt=delta/c0/sqrt(2);
fs = 1 / dt;

% Frequency Range
df=1/time(end);
freq = [0:length(time)-1].*df;

t_f=d_before+d_after;
number_of_time_steps=round(t_f/dt);
time_befr = 0:dt:d_before;                                                 % Before emission
time_emit = (d_before + dt):dt:(d_before + d_impulse);                     % During emission
time_aftr = (d_before + d_impulse + dt):dt:d_after;                        % After emission

time_ana = [time_befr, time_emit, time_aftr];

f_impulse= 1 / d_impulse;                                                  % f_impulse = 100 MHz
e_impulse0 = -0.5*(1-cos(2*pi*f_impulse*time_emit));
e_impulse1 = cos(2*pi*f0*time_emit);
e_impulse2 = e_impulse0 .* e_impulse1;
e_impulse2 = e_impulse2 * A_impulse;                                       % Analytical emitter impulse [V/m]

% Zero padding
zeros_befr = zeros(1, length(time_befr));
zeros_aftr = zeros(1, length(time_aftr));
e_impulse = [zeros_befr, e_impulse2, zeros_aftr];

% Scaling Factors
abs_e_t_cond = abs(e_t) > 0;                                               % Conditioner if abs(e_t) is greater than 0
N_e_t = sum(abs_e_t_cond == 1);                                            % Number of elements != 0 in e_t

abs_e_impulse_cond = abs(e_impulse) > 0;                                   % Conditioner if abs(r_t) is greater than 0
N_e_impulse = sum(abs_e_impulse_cond == 1);                                % Number of elements != 0 in r_t

E_f = fft(e_t);
E_f_scaled = E_f ./ N_e_t;
e_fft = fft(e_impulse) ./ N_e_impulse;

% Magnitude Verfication
f_t0 = e_impulse;
F_w_Int_Re = []; F_w_Int_Im = []; Mg_F_w = [];
time0 = reshape(time, 1, []);
for zeta = 1:length(freq)
    F_w_Re = f_t0 .* cos(2 .* pi .* freq(zeta) .* time0);
    F_w_Int_Re(zeta) = trapz(time0, F_w_Re);

    F_w_Im = f_t0 .* sin(2 .* pi .* freq(zeta) .* time0);
    F_w_Int_Im(zeta) = trapz(time0, F_w_Im);

    Mg_F_w(zeta) = sqrt(F_w_Int_Re(zeta).^2 + F_w_Int_Im(zeta).^2);
end

% Inverse Fourier Transform Verification
f_t_Mg_ifft = zeros(1, length(freq));
for zeta = 1:length(freq)
    f_t_Mg_ifft = f_t_Mg_ifft + ...
        (F_w_Int_Re(zeta) .* cos(2 .* pi .* freq(zeta) .* time0)) + ...
        (F_w_Int_Im(zeta) .* sin(2 .* pi .* freq(zeta) .* time0));
end
f_t_Mg_ifft = f_t_Mg_ifft ./ time0(end);

%% Analytical Fourier Transform

N = length(time) * 100;                                                    % Number of points []
fs = 1 / time(2);                                                          % Sampling frequency [Hz]
df = fs / N;                                                               % Spectral Resolution [Hz]
freq_ana = (0:(N - 1)) .* df;                                              % Frequency Range [Hz]
f0 = 2.4e9;                                                                % Carrier Frequency [Hz]
fi = 100e6;                                                                % Impulse Frequency [Hz]

% Analytical Solution of E(f)
A = 1;                                                                     % For analytical solutions
freq_w_t = freq_ana - (fs / 2);                                            % -(fs / 2) to +(fs / 2)

F_f_ana_inf = zeros(1, N);
% For 1 period of e(t) -> N_tau_w = 24: tau_w == d_impulse
N_period = 0.5;
N_tau_w = d_impulse * f0 * N_period;

% For f0
% Diracs
[min_F, F_Sf_ind1] = min(abs(freq_w_t - f0));
[min_F, F_Sf_ind2] = min(abs(freq_w_t + f0));
F_Sf_ana = zeros(1, N);
F_Sf_ana(F_Sf_ind1) = 1 / 2; F_Sf_ana(F_Sf_ind2) = 1 / 2;
F_f_ana_inf = F_f_ana_inf + F_Sf_ana;

% For (f0 + fi)
% Diracs
[min_F, F_Sf_ind1] = min(abs(freq_w_t - (f0 + fi)));
[min_F, F_Sf_ind2] = min(abs(freq_w_t + (f0 + fi)));
F_Sf_ana = zeros(1, N);
F_Sf_ana(F_Sf_ind1) = 1 / 4; F_Sf_ana(F_Sf_ind2) = 1 / 4;
F_f_ana_inf = F_f_ana_inf + F_Sf_ana;

% For (f0 - fi)
% Diracs
[min_F, F_Sf_ind1] = min(abs(freq_w_t - (f0 - fi)));
[min_F, F_Sf_ind2] = min(abs(freq_w_t + (f0 - fi)));
F_Sf_ana = zeros(1, N);
F_Sf_ana(F_Sf_ind1) = 1 / 4; F_Sf_ana(F_Sf_ind2) = 1 / 4;
F_f_ana_inf = F_f_ana_inf + F_Sf_ana;

% Window function
% tau_w = N_tau_w / (f0);
tau_w = mean([N_tau_w / (f0), N_tau_w / (f0 + fi), N_tau_w / (f0 - fi)])
[min_time, N_impulse] = min(abs(time - tau_w));
F_Wf_ana = tau_w .* sinc(freq_w_t .* (tau_w));

F_f_ana_inf = -(1 ./ 2) .* fftshift(F_f_ana_inf);
F_f_ana_fin = conv(F_f_ana_inf, F_Wf_ana, 'same') ./ tau_w;

% close all
% plot(freq_ana .* 1e-9, abs(F_Wf_ana), 'k')
% plot(freq_ana .* 1e-9, abs(F_f_ana_fin), 'k')

%% ------------------------------------------------------------------------
%                               PLOTTING
%  ------------------------------------------------------------------------

close all;

n_fig = 0;
n_fontsize = 14;

% Plot of the initial impulse
n_fig = n_fig + 1;
figure(n_fig);
ax(1) = subplot(2, 2, [1, 3]);
hold on;
grid on;
plot(time*1e9, e_t, 'b')
hold off;
set(gca,'xlim',[0 d_impulse].*1e9)
xlabel('Time [ns]','interpreter','latex')
ylabel('$E$ [V/m]','interpreter','latex')
title('Initial Impulse $e(t)$','interpreter','latex')

ax(2) = subplot(2, 2, [2, 4]);
hold on;
grid on;
% plot(time*1e9, e_impulse0, 'r')
% plot(time*1e9, e_impulse1, 'g')
plot(time_ana*1e9, e_impulse, 'r')
hold off;
set(gca,'xlim',[0 d_impulse].*1e9)
xlabel('Time [ns]','interpreter','latex')
ylabel('$E$ [V/m]','interpreter','latex')
title('Analytical Initial Impulse $s(t)$','interpreter','latex')

n_fig = n_fig + 1;
figure(n_fig);
hold on;
grid on;
plot(freq .* 1e-9, abs(E_f_scaled), 'b')
plot(freq .* 1e-9, Mg_F_w .* 1e8, 'r')
plot(freq_ana .* 1e-9, abs(F_f_ana_inf), 'g--')
plot(freq_ana .* 1e-9, abs(F_f_ana_fin), 'k')
hold off;
set(gca,'fontsize',n_fontsize,'xlim',[f0 - 4/d_impulse, f0 + 4/d_impulse]*1e-9)
% set(gca,'fontsize',n_fontsize,'xlim',[0 fs/2]*1e-9)
title('Fourier Transform of $e(t)$','interpreter','latex')
xlabel('Frequency [GHz]','interpreter','latex')
ylabel('$|F(f)|$','interpreter','latex')
legend('$E(f)$','$E(f) \quad Computational$','$E(f) \quad  Inf$', ...
    '$E(f) \quad  Fin$',...
    'interpreter','latex')

% n_fig = n_fig + 1;
% figure(n_fig);
% hold on;
% grid on;
% plot(time*1e9, ifft(E_f), 'b')
% plot(time*1e9, f_t_Mg_ifft, 'r')
% hold off;
% set(gca,'fontsize',n_fontsize)
% title('Inverse Fourier Transform')
% xlabel('Time [ns]','interpreter','latex')
% ylabel('$E$ [V/m]','interpreter','latex')
% legend('$ifft(E(f))$','$ifft(E(f)) \quad Computational$','interpreter','latex')


% close all;

% fig = figure(1);
% subplot(2, 1, 1)
% hold on;
% grid on;
% plot(time, f_t0(time));
% title('Signal'); xlabel('Time [s]'); ylabel('Amplitude []');
% 
% subplot(2, 1, 2)
% hold on;
% grid on;
% % xline(f0);
% plot(freq, abs(F_w_fft_scaled));
% plot(freq, Mg_F_w);
% hold off;
% set(gca,'xlim',[0 fs / 2])
% title('Fourier Transform'); xlabel('Frequency [Hz]'); ylabel('Magnitude []');
% legend('fft','Mg')
% 
% % curr_name = 'example5';
% % saveas(fig, strcat(pwd,'\images\analytics\',curr_name,'.jpg')) 
% 
% fig = figure(2);
% hold on;
% grid on;
% plot(time, f_t_ifft);
% plot(time, f_t_Mg_ifft);
% title('Signal from ifft');
% xlabel('Time [s]');
% ylabel('Amplitude []');
% legend('ifft','Mg ifft');
% 
% % curr_name = 'example5_ifft';
% % saveas(fig, strcat(pwd,'\images\analytics\',curr_name,'.jpg')) 

%% ------------------------------------------------------------------------
%                               FUNCTIONS
%  ------------------------------------------------------------------------

function d1 = dirac_1(x)
    d1 = dirac(x);
    idx = d1 == Inf;
    d1(idx) = 1;
end

