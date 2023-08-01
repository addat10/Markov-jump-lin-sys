clear all
close all
clc
addpath('../lib')
addpath(genpath('./connectivity_data_mahdi'))
dir1='888.CC';% 888.Base
dir2='100-200';

%TPM_path='./connectivity_data_mahdi/Transition_Probability.mat'; % First data
TPM_path=['./connectivity_data_mahdi/MC_probability/',dir1,'/',dir2,'/Transition_Probability.mat']; % First data
%% Define the MJLS
MJLS=example_single_integrator_mahdi_data(TPM_path);  
% MJLS=example_mahdi_double_int(TPM_path,1,1,0.1);
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
[U,eigs]=eig(MJLS.T);
figure
plot(eigs,'o')
title('Spectrum of the transition Matrix')
%% Compute the bound on the perturbation
pi_s=(1/sum(U(:,1)))*U(:,1);
T0=pi_s*ones(1,size(MJLS.T,1));
Delta_T=MJLS.T-T0;
eps_T=max(max(Delta_T))
Delta_T_1=norm(Delta_T,1)
%% Analysis with LMIs
eps_T=0.2;
cvx_tol=1e-3;
n=size(MJLS.As{1},1);

cvx_begin sdp
variable P(n,n) symmetric
cvx_precision high

exp_AtPA=zeros(n);
for i=1:size(MJLS.T,1)
    exp_AtPA=(pi_s(i)+eps_T)*MJLS.As{i}'*P*MJLS.As{i};
end

minimize 1

subject to:

exp_AtPA + cvx_tol*eye(n) <= P
P >= eye(n)
cvx_end
status=cvx_status; 
%% Check if the found P indeed satisfies the LMI
eig(P)
eig(exp_AtPA-P)