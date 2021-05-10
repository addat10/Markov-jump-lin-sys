clear all
close all
clc
addpath('../lib')
addpath('./connectivity_data_mahdi')
%% Define the MJLS
MJLS=example_9_mahdi_data();  
% p_ic happens to be the stationary distribution
%% Analyze the defined MJLS
% [B_cal,T_cal]=get_B_T_matrices(MJLS);
% tic
% rho_B=max(abs(eig(B_cal)));
% toc
% tic
% rho_T=max(abs(eig(T_cal)));
% toc
%% Analyze probability dynamics
 eigs=eig(MJLS.T);
 figure
 plot(eigs,'o')
%% Simulate MJLS
samples=100;
steps=30;
%ref=[5;4;3;2;1]*ones(1,steps)+1;
ref=1*[5;4;3;2;1]+1*(1:steps);
[x,p] = Simulate_MJLS(MJLS,steps,samples,ref);
%% Extract positions
N=MJLS.nx/2;
pos=zeros(N,steps,samples);
vel=zeros(N,steps,samples);
for i=1:samples
    pos(:,:,i)=x(1:N,:,i);
    vel(:,:,i)=x(N+1:end,:,i);
end
%% Plots
plot_ensemble_trajs(pos,samples)
%% Plot positions
plot_samples=samples;
figure()
for i=1:plot_samples
    for agent=1:N
        plot(1:steps,pos(agent,:,i))
        hold on
    end        
end
title('Positions')
% Plot velocities
figure()
for i=1:plot_samples
    for agent=1:N
        plot(1:steps,vel(agent,:,i))
        hold on
    end        
end
title('Velocities')
%% Plot probabilities
figure()
for i=1:steps
    plot(1:MJLS.N,p(:,i))
    hold on
end
