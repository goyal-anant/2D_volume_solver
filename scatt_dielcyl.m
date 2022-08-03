clc; clear all; close all;

%% defining constants
% epsilon0 = 8.853e-12;
epsilonr = 4;           %relative permitivity of the dielectric shell
epsilon  = epsilonr;     
theta0   = pi/2;        %angle between z-axis and propagation axis
phi0     = 0;           %angle between x-axis and projection of k on xy plane
lambda   = 1;           %wavelength 
E0       = 1;           %electric field amplitude 
k        = 2*pi/lambda; %wave vector
outer    = 0.30;        %outer radius
inner    = 0.25;        %inner radius
N        = 50;          %no. of pixels in the space

%% building the scatterer
%function [X1, Y1] = dielectric_shell(outer,inner,lambda,N,plot_flag)
[X1, Y1] = dielectric_shell(0.30,0.25,lambda,50,1);

%% formulating and solving the problem
a = lambda / 20; %side of the square patch
r = a/sqrt(pi);    %radius of equivalent circle with same cross section

%forming matrix C in the equation [C][E] = Ei
for n = 1:length(X1)
    for m = 1:length(X1)
        if n == m
            C(n,m) = 1 + ((epsilon - 1) * (1i/2) * (pi * k * r * besselh(1,2,(k * r)) - 2i));
        else
            rho_mn = sqrt((X1(n) - X1(m))^2 + (Y1(n) - Y1(m))^2);
            C(n,m) = (1i * pi * k * r / 2) * (epsilon - 1) * (besselj(1, (k * r))) * (besselh(0,2,(k * rho_mn)));
        end
    end
end

%for Ei vector
alpha = (X1 * sin(theta0) * cos(phi0)) + (Y1 * sin(theta0) * sin(phi0));
Ei = E0 * exp(-1i * k * alpha);     %incident electric field 

%solving matrix equation
E = C\Ei';      %total electric field on the dielectric volume

%% getting the plot of E in terms of phi
phi = atan(Y1./X1);
phi = sort(phi);
plot(phi*180/pi,abs(E));
% xlim([0 180])


% %% calculating the total field 
% phi = linspace(0,pi,228);   %for circular variation around the shell
% %rho0 = (0.3*lambda - 0.25*lambda) + 0.25*lambda;
% rho0 = 0.5*lambda;          %distance where we want to check the field
% 
% %temporary variable 'temp' to store the amplitude of the field
% temp = -1i * (pi * k / 2) * sqrt(2i / (pi * k * rho0)) ...
%     * exp(-1i * k * rho0) * (epsilon-1) * r * besselj(1,(k * r)) ;
% temp1 = cos(phi.')*X1 + sin(phi.')*Y1;  %to hold angular variation
% Es = temp * exp(1i * k * temp1) * E;    %required field
% Es = flip(abs(Es));                     %absolute value of the complex field
% 
% %% plotting the field w.r.t. phi variation
% AxesH = axes('XTick',0:10:180, 'NextPlot', 'add');
% plot(phi*180/pi,Es,'linewidth',3);
% hold on; grid on; set(gca,'fontsize',20)
% xlabel('\phi(degrees)');
% ylabel('|E|');
% 
% %% extracted from richmond figure 3 (reference plot)
% import_rr = load('extrd_fig3_rich.mat');
% rr = import_rr.rr;
% plot(rr(:,1),rr(:,2),linewidth = 2);
% 
% %% calculating the error between reference and implemented data
% error_vec = Es - rr(:,2);
% error = norm(error_vec);
% % plot(phi*180/pi,error_vec,linewidth = 2);