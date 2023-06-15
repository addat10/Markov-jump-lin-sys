% Find the largest robustness margin via bisection 
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
%% Bisection
eps_lims=[0;1];
bisect_tol=1e-2;
% Start the bisection and repeat until the limits are within bisect_tol
while eps_lims(2)-eps_lims(1)>bisect_tol
    eps_mid=mean(eps_lims); % bisect the current limits           
    [status,P]=verify_stab(eps_mid,MJLS,pi_s)
    if status % Choose the upper half as the new limits if feasible
        eps_lims(1)=eps_mid;
        P_ret=P;
    else % Choose the lower half as the new limits if infeasible
        eps_lims(2)=eps_mid;        
    end
end
eps_best=eps_lims(1);


function [status,P]=verify_stab(eps_mid,MJLS,pi_s)
    cvx_tol=1e-3;
    n=size(MJLS.As{1},1);
    
    cvx_begin sdp
    variable P(n,n) symmetric
    cvx_precision high
    
    exp_AtPA=zeros(n);
    for i=1:size(MJLS.T,1)
        exp_AtPA=(pi_s(i)+eps_mid)*MJLS.As{i}'*P*MJLS.As{i};
    end
    
    minimize 1
    
    subject to:
    
    exp_AtPA + cvx_tol*eye(n) <= P
    P >= eye(n)
    cvx_end
    if strcmp('Solved',cvx_status)
        status=true;
    else
        status=false;
    end
end
