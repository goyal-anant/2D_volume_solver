clc; clear all; close all;

%% plotting echowidth and E for different values of patch size

N = [50 100 150 200];
select_structure = 0; %0 for full shell and 1 for half-shell

AxesH = axes('XTick',0:10:200, 'NextPlot', 'add');
for i = 1:length(N)
    [phi, Wphi] = richmond1965(N(i),select_structure);
    kk = num2str(N(i));
    plot(phi*180/pi,Wphi,'linewidth',8,'DisplayName',strcat('patch size = \lambda/',num2str(N(i)/2)));
    hold on; 
end
grid on; set(gca,'fontsize',30);
xlabel('\phi(degrees)');
ylabel('Echo width/ Wavelength');
legend('FontSize',30)

% import_rr = load('extrd_fig4_rich.mat');
% rr = import_rr.fullcirc;
% plot(rr(:,1),rr(:,2),'--m','linewidth',10,'DisplayName',"published result");