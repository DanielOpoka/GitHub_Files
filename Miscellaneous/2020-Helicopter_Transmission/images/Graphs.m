clear; clc; close all;

WI=[5:50];
W1=2*(WI./(2.13));

R1=2.83; R2=17.16;
W2=(R1/R2).*W1;

R3=2.36; R5=4.25;
W4=(R3/2*(R3+R5))*W2;

plot(W1,W4);