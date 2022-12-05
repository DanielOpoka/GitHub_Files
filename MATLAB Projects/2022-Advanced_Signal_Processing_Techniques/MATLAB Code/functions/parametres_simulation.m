function [ delta,dt,t_H,time,f_max,f_0,lambda_min ] = parametres_simulation(  c0,t_f,d_impulse,f0,Lx,Ly,modulation_ok,fact_delta )

f_impulse=1/d_impulse;

if modulation_ok==0
    lambda_min=c0/(2*f_impulse);
    f_min=0;
    f_max=2*f_impulse;
    BW=f_max-f_min;
end

if modulation_ok==1
    lambda_min=c0/(f0+2*f_impulse);
    f_min=f0-2*f_impulse;
    f_max=f0+2*f_impulse;
    BW=f_max-f_min;    
end

% Deltas
delta=lambda_min/fact_delta;
dt=delta/c0/sqrt(2);

% Vecteur temps
number_of_time_steps=round(t_f/dt);
time=dt*[0:number_of_time_steps-1].' ;


% Temps d'Heisenberg
f_0=(f_min+f_max)/2;
% t_H=((pi^2)*Lx*Ly/(c0^2))*f_0; % ANCIEN
% t_H=(16/(c0^2))*(Lx*Ly)*f_0;
t_H=4*pi*Lx*Ly/(c0^2)*f0;

% L=sqrt((Lx^2)+(Ly^2));
% S=Lx*Ly;
% t_H=((2*pi*S*f_0)/(c0^2))*(((2*f0*S+c0*L)/(2*f0*S+c0*L/2))^2)+(pi*S/c0*(f0)^2)*((4*(S^2)*c0*L*f0+2*(c0^2)*(L^2)*S)/((2*S*f0+c0*L/2)^3))



end