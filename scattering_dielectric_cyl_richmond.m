clc; clear variables; close all;

%% defining constants
% epsilon0 = 8.853e-12;
epsilonr = 4;           %relative permitivity of the dielectric shell
epsilon = epsilonr;     
theta0 = pi/2;          %angle between z-axis and propagation axis
phi0 = 0;               %angle between x-axis and projection of k on xy plane
lambda = 1;             %wavelength 
E0 = 1;                 %electric field amplitude 
k = 2*pi/lambda;        %wave vector

%% defining the scatterer
shellsize = lambda;     %temporary space to carve out the dielectric shell
N = 150 ;                 %number of pixels in which space is broken
x1    = linspace(-shellsize,shellsize,N);   %x-space of square space
y1    = linspace(-shellsize,shellsize,N);   %y-space of square space
[X Y] = meshgrid(x1,y1);
X     = reshape(X,[(N)^2,1]);           %to generate the coordinate system
Y     = reshape(Y,[(N)^2,1]);           %to generate the coordinate system

%plotting the square space
 scatter(X,Y,'k','filled');
 hold on; grid on; axis('equal');

%carving out dielectric shell from the square space generated by X and Y
t = 0;      %temporary variable to hold index 
for m = 1:(N)^2
    if ((X(m)^2 + Y(m)^2) <= (0.3*lambda)^2) && ((X(m)^2 + Y(m)^2) >= (0.25*lambda)^2)
       t = t + 1; 
       X1(t) = X(m);
       Y1(t) = Y(m);
    end
end

%plotting the dielectric shell
 scatter(X1,Y1,'w','filled');
 grid on; axis('equal'); hold off;
set(gca,'fontsize',20)
%% formulating and solving the problem
a = lambda / 55; %side of the square patch
r = 0.56 * a;    %radius of equivalent circle with same cross section

%forming matrix C in the equation [C][E] = Ei
for m = 1:length(X1)
    for n = 1:length(X1)
        if m == n
            C(m,n) = 1 + ((epsilon - 1) * (1i/2) * (pi * k * r * besselh(1,2,(k * r)) - 2i));
        else
            rho = sqrt((X1(m) - X1(n))^2 + (Y1(m) - Y1(n))^2);
            C(m,n) = (1i * pi * k * r / 2) * (epsilon - 1) * (besselj(1, (k * r))) * (besselh(0,2,(k * rho)));
        end
    end
end

%for Ei vector
alpha = X1 * sin(theta0) * cos(phi0) + Y1 * sin(theta0) * sin(phi0);
Ei = E0 * exp(-1i * k * alpha);     %incident electric field 

%solving matrix equation
E = C\Ei';      %total electric field on the dielectric volume

%% plotting the total field 
phi = linspace(0,pi,256);   %for circular variation around the shell
% rho0 = ((0.3*lambda - 0.25*lambda)/2) + (0.25*lambda);
rho0 = 0.3*lambda;          %distance where we want to check the field

%temporary variable 'temp' to store the amplitude of the field
temp = -1i * (pi * k / 2) * sqrt(2i / (pi * k * rho0)) ...
    * exp(-1i * k * rho0) * (epsilon-1) * r * besselj(1,(k * r)) ;
temp1 = cos(phi.')*X1 + sin(phi.')*Y1;  %to hold angular variation
Es = temp * exp(1i * k * temp1) * E;    %required field
Es = abs(Es);                           %absolute value of the complex field

%plotting the field w.r.t. phi variation
AxesH = axes('XTick',0:10:180,'YTick', 0:0.25:2.5, 'NextPlot', 'add');
plot(phi*180/pi,flip(Es),'linewidth',3);
grid on; set(gca,'fontsize',20)
xlabel('\phi(degrees)');
ylabel('|E|');