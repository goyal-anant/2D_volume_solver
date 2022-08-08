clc; clear all; close all;

%% for plotting dielectric half shell for different values of N on same figure


N = [25 50 75 100 125 150 175 200 250];
for i = 1:length(N)
    [X1 ,Y1 ,r ,X ,Y ] = dielectric_half_shell(0.30,0.25,1,N(i),0);

    subplot(3,3,i)
    %plotting the square space
    scatter(X,Y,'k','filled');
    hold on; grid on; axis('equal','tight');

    %plotting the dielectric shell
    scatter(X1,Y1,'red','filled');
    axis('equal'); hold off;
    title("N = " + N(i) + ", patch size = \lambda/" + N(i)/2)
    set(gca,'fontsize',20)
end

