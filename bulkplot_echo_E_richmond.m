clc; clear all; close all;

%% plotting echowidth and E for different values of patch size

N = [25 50 75 100 125 150 175 200 250];

for i = 1:length(N)
    richmond1965(N(i));
    hold on;
end

