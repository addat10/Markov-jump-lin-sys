clear all
close all
clc
addpath('../lib')
addpath(genpath('./connectivity_data_mahdi'))

dir1='2458.CC';% 888.Base
dir2='300-400';

%TPM_path='./connectivity_data_mahdi/Transition_Probability.mat'; % First data
TPM_path=['./connectivity_data_mahdi/MC_probability/',dir1,'/',dir2,'/Transition_Probability.mat']; % First data

%% Define the MJLS
MJLS=example_9_mahdi_data(TPM_path);  
% p_ic happens to be the stationary distribution
%% Analyze the defined MJLS
[B_cal,T_cal]=get_B_T_matrices(MJLS);
tic
rho_B=max(abs(eig(B_cal)));
toc
tic
rho_T=max(abs(eig(T_cal)));
toc
%% Analyze probability dynamics
eigs=eig(MJLS.T);
figure
plot(eigs,'o')
title('Spectrum of the transition Matrix')
%% Visualize the Transition probability Matrix
figure()
plot_matrix_data(MJLS.T,'Transition Probability Matrix')
title('Transition Probability Matrix')
 
%% Simulate MJLS

% Generate Reference
%ref=[5;4;3;2;1]*ones(1,steps)+1;
%ref=1*[5;4;3;2;1]+1*(1:steps);
[ref_leader]=gen_ref(MJLS.dt);
ref=1*[5;4;3;2;1]+ref_leader;

% Set appropriate initial conditions
MJLS.x_ic=[ref(:,1);20*ones(5,1)];

% Time-stepping
steps=size(ref,2);
samples=1;
[x,p] = Simulate_MJLS(MJLS,steps,samples,ref);
%% Extract positions
N=MJLS.nx/2;
pos=zeros(N,steps,samples);
vel=zeros(N,steps,samples);
acc=zeros(N,steps,samples);
for i=1:samples
    pos(:,:,i)=x(1:N,:,i);
    vel(:,:,i)=x(N+1:end,:,i);
    acc(:,1,i)=zeros(N,1);
    acc(:,2:end,i)=(1/MJLS.dt)*[x(N+1:end,2:end,i)-x(N+1:end,1:(end-1),i)];    
end
%% Plots
plot_ensemble_trajs(pos,samples)
%% Plot states

% Plot Positions
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

% Plot Accelerations
figure()
for i=1:plot_samples
    for agent=1:N
        plot(1:steps,acc(agent,:,i))
        hold on
    end        
end
title('Acceleration')
%% Plot probabilities
figure()
for i=1:steps
    if i==steps
        plot(1:MJLS.N,p(:,i),'o')
    else
        plot(1:MJLS.N,p(:,i))
    end
    hold on
end
title('Dynamics of probability distributions')
