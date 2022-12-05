
%% ------------------------------------------------------------------------
%                             REVERBERATION TIME
%  ------------------------------------------------------------------------

% By simulating a reverberant chamber, its reverberation time can be
% calculated. Different approaches are used to calculated it, such a
% uniform distribution or a random distribution of receivers. 

% The equation used is from:
% Lerosey, G. (2006). Retournement temporel d'ondes électromagnétiques et application à la télécommunication en milieux complexes (Doctoral dissertation, ESPCI ParisTECH).


clear all; close all; clc; format long;
addpath("functions")
addpath("data");

new_dir = "images/ReverberationTime";
n_fontsize = 14;

%% SIMULATION PARAMETERS

% Impulse parameters
signal_type = "wavelet";
d_impulse=8e-9;                                                            % Duration of the impulse [s]
f0=2.4e9;                                                                  % Carrier frequency of the impulse [Hz]
d_before=0e-9;                                                             % delay before impulse [s]
d_after = 750e-9;                                                          % delay array after impulse [s]
A_impulse=1;                                                               % Amplitude of the impulse [V/m]

% Cavity parameters
Lx=1.9; Ly=1.9;                                                            % Cavity dimensions
sigma=(6.5)*1e2;                                                           % Conductivity of the walls [S/m]
coef_losses_ok=1;                                                          % If=1: losses, if=0: PEC condition 
coef_regular_chao=2;                                                       % If=1: regular cavity, if=2: chaotic cavity 

% WARNING: lossless chaotic cavity not coded (possible: lossless or with losses reglar cavity, or with losses chaotic cavity)
[time, e_t, param_simu] = function_parameters_simu(Lx,Ly,d_impulse,A_impulse,f0,d_before,d_after,coef_regular_chao,coef_losses_ok,sigma, signal_type);

% Wave parameters
c0=3e8;                                                                    % Speed of light [m/s]
f_impulse=1/d_impulse;
lambda_min=c0/(f0+2*f_impulse);
fact_delta=10;
delta=lambda_min/fact_delta;
dt=delta/c0/sqrt(2);
fs = 1 / dt;                                                               % Sampling Frequency [Hz]

%% For tau_RT vs. different emitter positions

N_recs = 500;                                                              % Number of receivers [#]
corr_delta = lambda_min / 2;                                               % Distance at which correlation changes [m]
position_receiver_all = zeros(2, N_recs + 1);                              % Location of receivers [m]

low_X = 0;                                                                 % Lower boundary for random X point [m]
upp_X = 2.2;                                                               % Upper boundary for random X point [m]

low_Y = corr_delta;                                                        % Lower boundary for random Y point [m]
upp_Y = 1.8;                                                               % Upper boundary for random Y point [m]

emt_d = 1.75;
pos_emt = [(emt_d / 4), (emt_d / 2), ((3 * emt_d) / 4), ...
    (emt_d / 4), (emt_d / 2), ((3 * emt_d) / 4), ...
    (emt_d / 4), (emt_d / 2), ((3 * emt_d) / 4); ...
    (emt_d / 4), (emt_d / 4), (emt_d / 4), ...
    (emt_d / 2), (emt_d / 2), (emt_d / 2), ...
    ((3 * emt_d) / 4), ((3 * emt_d) / 4), ((3 * emt_d) / 4)];

n_emt = size(pos_emt, 2);

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

nb_fig = 0;

close all

% To check positioning
nb_fig = nb_fig + 1;
figure(nb_fig);
openfig('data/CavitySpace/cavity_space_empty.fig');

hold on;
viscircles(cavity_center, cavity_radius, 'Color', 'red');
% viscircles(emit_center, emit_radius, 'Color', 'green');
plot(pos_rec(1, :), pos_rec(2, :), 'rx')
for ind_emt = 1:n_emt
    plot(pos_emt(1, ind_emt), pos_emt(2, ind_emt), 'g*')
    text(pos_emt(1, ind_emt), pos_emt(2, ind_emt), " "+ num2str(ind_emt), ...
        'Color', 'green', 'FontSize', n_fontsize);
end
hold off;

figure_title = new_dir + "/Cavity-" + num2str(N_recs) + "Recs.png";
saveas(gcf, figure_title)

tau_RT_emt = zeros(1, n_emt);

% Display of the impulse respone
Display_rep_imp = 0;                                                       % If=1: display, if=0: no display
    
for ind_emt = 1:n_emt
    % Received Signal
    [r_t, Ez] = function_impulse_response(param_simu, pos_emt(:, ind_emt), ...
        pos_rec, Display_rep_imp, time, e_t);
    tau_RT_emt(ind_emt) = ReverberationTime(time, r_t)
end

tau_RT_emt_plot = round(tau_RT_emt, 10) .* 1e9;

% To check positioning
% nb_fig = nb_fig + 1;
% figure(nb_fig);
openfig('data/CavitySpace/cavity_space_alex.fig');

hold on;
viscircles(cavity_center, cavity_radius, 'Color', 'red');
% viscircles(emit_center, emit_radius, 'Color', 'green');
plot(pos_rec(1, :), pos_rec(2, :), 'rx')
for ind_emt = 1:n_emt
    plot(pos_emt(1, ind_emt), pos_emt(2, ind_emt), 'g*')
    text(pos_emt(1, ind_emt), pos_emt(2, ind_emt), " " + num2str(tau_RT_emt_plot(ind_emt)) + "ns", ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'Color', 'green', 'FontSize', n_fontsize);
end
hold off;

title("Reverberation Time of Cavity with " + num2str(n_emt) + " emitters", 'interpreter', 'latex');

figure_title = new_dir + "/Cavity-" + num2str(N_recs) + "Recs-tau_RT.png";
saveas(gcf, figure_title)

%% For tau_RT vs. Sim. Time

tau_RT_array = [];                              % Array for Reverberation Time [s]

d_after_array = (10:10:1000) .* 1e-9;           % delay array after impulse [s]

for i = 1:length(d_after_array)
    d_after = d_after_array(j);
    [time, e_t, param_simu] = function_parameters_simu(Lx,Ly,d_impulse,A_impulse,f0,d_before,d_after,coef_regular_chao,coef_losses_ok,sigma);

    % Receiver and emitter locations
    position_receiver=[0.6;0.2];                % Location of the receiver [m]
    position_emitter=[1.2;1.3];                 % Location of the emitter [m]

    % Varying position for receiver
    % % Location of the central emitter [step]
    position_receiver_N_save=round(position_receiver/delta); 
    Coef_receiver=1;                            % Number of receivers in both directions (+ and - x) from the central receiver
    m=1;
    for j = -Coef_receiver:1:Coef_receiver
        position_receiver_N=position_receiver_N_save+[j;0];
        position_receiver(:,m)=position_receiver_N*delta;
        m=m+1;
    end

    % Display of the impulse respone
    Display_rep_imp = 0;                         % If=1: display, if=0: no display
    [r_t,Ez] = function_impulse_response(param_simu, position_emitter, position_receiver, Display_rep_imp, time, e_t );

    tau_RT_array = [tau_RT_array, ReverberationTime(time, r_t)]
end

%% Receiver positions - Uniform Distribution

CR_rows = 10;
CR_cols = 20;

% Location of receivers
position_receiver_all = zeros(2, CR_rows * CR_cols);
position_receiver_all(1, 1) = corr_delta * 10;
position_receiver_all(2, 1) = corr_delta * 5;

position_emitter=[1.2;1.3];                 % Location of the emitter [m]

i = 0;
for j = 1:CR_rows
    for k = 1:CR_cols
        i = i + 1;
        position_receiver_all(1, i) = position_receiver_all(1, 1) + (k - 1) * corr_delta;
        position_receiver_all(2, i) = position_receiver_all(2, 1) + (j - 1) * corr_delta;
    end
end

%% tau_RT_end vs CR

N_CR = 100;                                 % Number of Receivers
CR_array = 1:N_CR;
tau_RT_end = zeros(1, N_CR);

for i = 1:N_CR
    [time, e_t, param_simu] = function_parameters_simu(Lx,Ly,d_impulse,A_impulse,f0,d_before,d_after,coef_regular_chao,coef_losses_ok,sigma);
    
    position_receiver = [];
    
    position_receiver(1, 1:i) = position_receiver_all(1, 1:i);
    position_receiver(2, 1:i) = position_receiver_all(2, 1:i);
    
    % Display of the impulse respone
    Display_rep_imp = 0;                         % If=1: display, if=0: no display
    [r_t,Ez] = function_impulse_response(param_simu, position_emitter, position_receiver, Display_rep_imp, time, e_t );

    tau_RT_end(i) = ReverberationTime(time, r_t)
end


%% ------------------------------------------------------------------------
%                       REVERBERATION TIME POST-PROCESSING
%  ------------------------------------------------------------------------

% Parameter: Coef_receiver = 5
% Parameter: coef_regular_chao=2

% d_after_array = [0, d_after_array];
% tau_RT_array = [0, tau_RT_array];

% Naming convention for saved files:
% CR -> Coef_receiver, CRC -> coef_regular_chao
% Pos 1 = position_receiver=[0.6;0.2];
% Pos 2 = position_receiver=[1.2;0.2];
% Pos 3 = Uniform distribution at intervals of 5*delta;
% Pos 4 = Random distribution at intervals of 5*delta;
% Last numbers are the final simulation time

% % For tau_RT vs. Sim. Time
% file_name_exp1 = 'data/tau_RT_data_CR01_Pos1_CRC2-1e3.mat';  
% % save(file_name_exp1,'d_after_array','tau_RT_array')
% load(file_name_exp1);
% 
% % For tau_RT(end) vs CR
% file_name_exp2 = 'data/tau_RT_data_1-100CR_Pos4_CRC2-1e3.mat';    
% save(file_name_exp2,'CR_array','tau_RT_end')
% load(file_name_exp2);

%% Plotting
% 
% close all;
% 
% n_fontsize = 14;
% n_fig = 0; % For an ever increasing number of plots
% 
% n_fig = n_fig + 1;
% figure(n_fig);
% figure(1);
% hold on;
% grid on;
% plot(d_after_array .* 1e9, tau_RT_array .* 1e9, 'LineWidth', 2)
% 
% plot(d_after_array(20) .* 1e9, tau_RT_array(20) .* 1e9, 'r*', 'MarkerSize', 12)
% plot(d_after_array(50) .* 1e9, tau_RT_array(50) .* 1e9, 'g*', 'MarkerSize', 12)
% plot(d_after_array(100) .* 1e9, tau_RT_array(100) .* 1e9, 'b*', 'MarkerSize', 12)
% 
% plot([200 , 200], [0, tau_RT_array(20) .* 1e9], 'r--')
% plot([500 , 500], [0, tau_RT_array(50) .* 1e9], 'g--')
% plot([1000 , 1000], [0, tau_RT_array(100) .* 1e9], 'b--')
% 
% yline(tau_RT_array(end) * 1e9, '--');
% % plot(d_after_array .* 1e9, tau_RT_fit_smooth .* 1e9);
% hold off;
% xlim([0, d_after_array(end)] .* 1e9);
% set(gca, 'fontsize', n_fontsize);
% % title('Reverberation Time $\tau_{RT}$ vs. Simulation Time','interpreter','latex');
% xlabel('Simulation Time [ns]', 'interpreter', 'latex');
% ylabel('Reverberation Time $\tau_{RT}$ [ns]', 'interpreter', 'latex');
% % legend_label = strcat('Coeff Receiver=',strcat(' ',file_name_exp1(20:21)));
% legend('', '200 ns' , '500 ns', '1000 ns', 'Location', 'best')
% 
% 
% n_fig = n_fig + 1;
% figure(n_fig);
% hold on;
% grid on;
% plot(CR_array, tau_RT_end .* 1e9, 'LineWidth', 2);
% % xlim([0, 163]);
% set(gca,'fontsize',n_fontsize);
% % title('Reverberation Time $\tau_{RT}$ vs. Number of Receivers','interpreter','latex');
% xlabel('Number of Receivers []','interpreter','latex');
% ylabel('$\tau_{RT}(end)$ [ns]','interpreter','latex');
% legend(file_name_exp2(26:29), 'Location', 'best')
% % legend('random', 'uniform', 'Location', 'best')
