
%% Kinematics of box-pushing mechanism

% Knowing the dimensions of the mechanism, equations describing the
% angular velocities are developped such as to obtained the speeds found
% at the end of the mechanism that help in understanding the overall
% kinematics of the system.
% This is accomplished by solving the equations using Euler's forward 
% method as it could be solved by simple linear algebra.

clear; clc; close all;

hold on;

A=0.464;

% Displacement
Pa2=-pi:pi/50:pi;
X=(-1.2*sin(Pa2)+sqrt(((1.2*sin(Pa2)).^2-(4)*(0.36-0.48*cos(Pa2)))))/(2); 
X1=(-1.2*sin(Pa2)-sqrt(((1.2*sin(Pa2)).^2-(4)*(0.36-0.48*cos(Pa2)))))/(2); 

figure(1);
hold on;
grid;
plot(Pa2,X);
plot(Pa2,X1);
xlabel('Psi 2 [radians]');
ylabel('Displacement [meters]');
plot(0,0.3464,'or');

legend('Displacement with +sqrt','Displacement with -sqrt','Displacement when Psi 2 = 0');

%Euler's forward method
h=0.1;
P1(1)=-pi;
P2(1)=-0.23914;
SR(1)=Speed_Ratio(P1(1),P2(1));
Psi_1_dot=1; %1 rad/s
Psi_2_dot(1)=SR(1)*Psi_1_dot;
V(1)=abs(Velocity(P2(1))*Psi_2_dot(1));

N=length([P1(1):h:pi]);

for n=1:N
    P2(n+1)=P2(n)+h*(Output_Input_Relation(P1(n),P2(n)));
    P1(n+1)=P1(1)+n*h;
    SR(n+1)=Speed_Ratio(P1(n+1),P2(n+1));
    Psi_2_dot(n+1)=SR(n+1)*Psi_1_dot;
    V(n+1)=(Velocity(P2(n+1))*Psi_2_dot(n+1));
end

axis equal;

figure(2);
hold on; grid;
plot(P1,P2,'b');
plot(0,0,'or');
ylim_P2=get(gca);
xlabel('Psi 1 (Input angle) [radians]'); ylabel('Psi 2 (Output angle) [radians]');

% Displacement:
X=(-1.2*sin(P2)+sqrt(((1.2*sin(P2)).^2-(4)*(0.36-0.48*cos(P2)))))/(2);

figure(3);
hold on; grid;
plot(P2,X,'b');
ylim_P2=get(gca);
xlabel('Psi 2 [radians]'); ylabel('Displacement [meters]');
plot(0,0.3464,'or');

% Angular Speed Ratio:
figure(4);
hold on; grid; 
plot(P1,SR);
xlabel('Psi 1 [radians]'); ylabel('Speed ratio: Psi 2 dot / Psi 1 dot');

figure(5);
hold on; grid;
plot(P1,V);
xlabel('Psi 1 [radians]'); ylabel('Velocity [meters/seconds]');

function f = Output_Input_Relation(P1,P2)
A=0.464; 
f = (0.02*cos(P1)+0.01*sin(P1)+0.0224*sin(A+P2-P1))/(0.0224*sin(A+P2-P1)+0.0896*cos(A+P2)+0.0448*sin(A+P2));
end

function SR = Speed_Ratio(P1,P2)
A=0.464;
SR=(0.0224*sin(A+P2-P1)+0.1*(2*cos(P1)+sin(P1)))/(0.0896*cos(A+P2)+0.0448*sin(A+P2));
end

function V1 = Velocity(P2)
V1=(0.5.*(-1.2.*cos(P2)+((1.44.*sin(2*P2)-1.92.*sin(P2))/2.*sqrt(1.44.*(sin(P2)).^2-4.*(0.36-0.48.*cos(P2))))));
end
