clear; clc; close all; format SHORTENG;

%%   *** GEAR GEOMETRIES ***
% % LEGO Gear of 21 teeth geometry 
% % ** Same gear for pinion and wheel **
z = [21 21]; % Number of gear teeth
m0 = 0.5; % Module [mm]
alfa0 = 20; % Pressure angle [degrees]
b = [8 8]; % Face width [mm]
x = [0 0]; % Profile shift coefficients

Ca = [1 1]; % Addendum coefficient
ha = [(Ca(1)+x(1))*m0 (Ca(2)+x(2))*m0]; % Tooth addendum
Cf = [1.25 1.25]; % Dedendum coefficient
hf = [(Cf(1)-x(1))*m0 (Cf(2)-x(2))*m0]; % Tooth dedendum

r = (m0.*z)/2; % Pitch radius
rf = r-hf; % Base radius
ra = r+ha; % Outside radius
pb = 2*pi*(rf(1)./z(1)); % Base pitch

%%   *** TORSIONAL STIFFNESS MODELING ***
% Lego Shaft dimensions
D = 4.72e-3;
t1 = 1.83e-3;
t2 = (D-t1)/2;
L = 95.2e-3;

% Lego shaft material properties
E = 2.588e9; % Young's Modulus [GPa]
G = 875e6; % Shear Modulus [MPa]
rho = [1030 1070]; % Density [kg/m^3]

% Experimental Analysis
% Lego Lever
ml = [3.83 .52 .22 .28].*1e-3; % Mass of lego components [kg]
dl = [59.65 15.5 24 113].*1e-3; % CoG of lego components
g = 9.81; % Constant acceleration of gravity [m/s^2]
Tl = (sum(ml.*g.*dl)); % Sum of torques of lego components [Nm]

d = .1113; % Length of lever [m]
mx = [.15 .3 .45 .6 .7 .75 .9]; % Mass of Water [kg]
m = (.170+mx); % mass of water bottle [kg]
W = m*9.81; 
T = W*d+Tl;
phi_deg = [57 67 73 74 78 78 80];
phi_rad = phi_deg.*(pi/180);
k_exp1 = T./phi_rad; % Experimental torsional stifness [Nm]
k_exp = mean(k_exp1);

% Analytical analysis
% Analytical polar moment of inertia
J1 = (pi*(D/2)^4)/2; % Version 1: Circular [m^4]
J2 = (1/6)*(D^4-16*t2^4); % Version 2: Rectangular minus spaces [m^4]
J3 = (pi/32)*(D^4-4*t2^4); % Version 3: Circular minus spaces [m^4]
J = [J1 J2 J3];
k_ana = J.*G./L; % Analytical torsional stiffness [Nm/rad]
per_error_tor = ((k_exp-k_ana)./k_exp).*100; % Percentage error [%]
[~,ind] = min(per_error_tor);
J_ana = J(ind); % ·Best" analytical polar moment of inertia
k_ana1 = k_ana(ind);

%   *** EV3 DATA PROCESSING ***
FN = dir(fullfile('LegoEV3_Data', '\*.csv'));
names = {FN.name}; n = numel(names);
names_omega = {'omega_O_M','omega_0_T','omega_4_5_T','omega_9_0_T','omega_1_3_5_T','omega_1_8_0_T'};

% Assumed torque from https://www.philohome.com/motors/motorcomp.htm
T1 = 6.64e-2; % Torque of Lego motor [Nm]

%% *** DATA PROCESSING ***
figure(1);
hold on;
xlim([0 60]); ylim([160 210]); grid;
title('Motor Rotational Speed vs. Time');
xlabel('Time [s]'); ylabel('Motor Rotational Speed [RPM]');
for i = 1:n
    data{i} = dlmread(['LegoEV3_Data/', names{i}],'R');
    data1 = data{i};
    time{i} = data1(:,1);
    omega11{i} = data1(:,2);
    omega1{i} = data1(:,2).*(pi/30); % Rotational Speed [rad/s]
    OM{i} = T1.*omega1{i};
    plot(time{i},omega11{i},'LineWidth',2);
end
legend(names_omega{:},'Location','bestoutside');
hold off;

time_in = time{1};
omega_in = omega1{1};

%% *** POWER LOSSES MODEL ***
% Gear geometry that affects power loss due to friction:
gf = sqrt(ra(1)^2-(r(1)*cos(alfa0))^2)-r(1)*sin(alfa0); % Length of approach path
ga = sqrt(ra(2)^2-(r(2)*cos(alfa0))^2)-r(2)*sin(alfa0); % Length of recess path

n_fr = (numel(OM{1})-numel(OM{2})-1); Pin_fr = OM{1}; Pfr = OM{2};
if n_fr > 0
    Pin_fr(end-n_fr:end) = [];
    i2 = 2;
elseif n_fr < 0
    n_fr = abs(n_fr);
    Pfr(end-n_fr+2:end) = [];
    i2 = 1;
end
% The experimental power lost from the motor when engaged with the test rig:
PLfr_exp = Pin_fr-Pfr;

% Range of values for mean friction coefficient of ABS plastic:
% Obtained from: https://plastics.ulprospector.com/generics/1/c/t/acrylonitrile-butadiene-styrene-abs-properties-processing
r_miu_m = [0.11 0.46];
miu_m = r_miu_m(1); % Mean friction coefficient
n1 = 0.01; % Accuracy for friction coefficient

per_error_PLfr = 100; 

% Programming construct to obtain the best mean friction coefficient from the given range:
while (per_error_PLfr >= 1)
    % The analytical friction power loss:
    PLfr_ana = pi*((1/z(1))+(1/z(2)))*(1-((gf+ga)/pb)+(gf/pb)^2+(ga/pb)^2).*Pin_fr.*miu_m;
    
    per_error_PLfr = ((PLfr_exp-PLfr_ana)./(PLfr_exp))*100;
    per_error_PLfr = mean(per_error_PLfr,'omitnan');
    
    if miu_m >= r_miu_m(2)
        break;
    end
    miu_m = miu_m+n1;
end

figure(2);
hold on; xlim([0 60]); grid;
title('Friction Power Loss vs. Time');
xlabel('Time [s]'); ylabel('Power Loss [W]');
plot(time{i2},PLfr_exp,'LineWidth',2);
plot(time{i2},PLfr_ana,'--','LineWidth',2);
namesfr = {'PLfr_e_x_p','PLfr_a_n_a'}; legend(namesfr);
hold off;

Tfr_exp = mean(PLfr_exp)/mean(omega1{1});
Tfr_ana = mean(PLfr_ana)/mean(omega1{1});

% System Power Losses
% Test rig dimensions:
L1 = 43e-3; L1d = 48e-3; L2 = 98e-3; L2d = 122e-3;  % Effective length of shafts
LTDS = [130e-3 130e-3 130e-3 131e-3]; % Effective length of TDS for different configurations
phi = z(1)/z(2); % Gear ratio 
j1 = 0;

figure(3); 
hold on; namesS = []; xlim([0 60]); ylim([0 0.2]); grid;
colors = {'b','r','g','m'};
title('System Power Losses vs. Time');
xlabel('Time [s]'); ylabel('Power Loss [W]');
for i1 = 3:(numel(names))
    j1 = j1+1;
    Pin_T = OM{1}; % Input mechanical power [W]
    P_T = OM{i1};
    n_T = (numel(Pin_T)-numel(P_T)-1);
    if n_T > 0
        Pin_T(end-n_T:end) = [];
        i3 = i1;
    elseif n_T < 0
        n_T = abs(n_T);
        P_T(end-n_T+2:end) = [];
        i3 = 1;
    end
    PLS_exp{j1} = Pin_T-P_T;
    TS_exp(j1) = mean(PLS_exp{j1})/mean(omega1{1})-Tfr_exp; 
    thetac = (pi/4).*j1; % thetac = [35, 80, 125, 200].*(pi/180);
    keq(j1) = (3*J_ana*G)/(3*L2d+3*phi*L1d+phi*LTDS(j1)+3*phi*L1+3*L2); % Equivalent torsional stiffness of model test rig [Nm/rad]
    phic(j1) = ((L2d*(LTDS(j1)-3*L1))/(L2*(3*L1d+LTDS(j1))))*thetac; % Total angle of twist [rad]
    TS_ana(j1) = keq(j1)*phic(j1);
    PLS_ana{j1} = OM{1}-((T1-Tfr_ana-TS_ana(j1)).*omega1{1});
    per_error_PLS(j1) = ((mean(PLS_exp{j1})-mean(PLS_ana{j1}))/(mean(PLS_exp{j1})))*100; 
    
    plot(time{i3},PLS_exp{j1},colors{j1},'LineWidth',2);
    plot(time{1},PLS_ana{j1},['--',colors{j1}],'LineWidth',2);
    
    PLS_m_exp(j1) = mean(PLS_exp{j1},'omitnan'); 
    PLS_m_ana(j1) = mean(PLS_ana{j1},'omitnan');
    
    thetac_d = 45*j1;
    namesS = [namesS,"PLS_e_x_p-"+thetac_d+"T","PLS_a_n_a-"+thetac_d+"T"];
end
legend(namesS{:},'Location','eastoutside');
hold off;

%% *** RESULTS ***

Tfr_exp
miu_m
Tfr_ana

TS_exp
TS_ana

PLS_m_exp
PLS_m_ana

per_error_PLfr
per_error_PLS
mean_per_error_PL = mean([per_error_PLfr,per_error_PLS])

