clear all
close all
clc
addpath('./examples/')
%% Define the MJLS
% Eaxmple 1: scaler not MSS
% Example 2: scalar MSS
% Example 3: Unstable+Unstable=Stable  
% Example 4: Stable+Stable=Unstable
MJLS=example_7_network_stabilization_toy();      
%% Analyze the defined MJLS
[B_cal,T_cal]=get_B_T_matrices(MJLS);
rho_B=max(abs(eig(B_cal)));
rho_T=max(abs(eig(T_cal)));