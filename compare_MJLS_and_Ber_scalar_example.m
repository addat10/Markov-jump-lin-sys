% This script is used to generate random examples and compare the stability
% of the Markov model and the corresponding iid model with the stationary
% distribution (Bernoulli)
clear all
close all
clc
addpath('./examples/')
addpath('./lib/')

rng(20)
no_samples=1000;
rho_Ber=zeros(no_samples,1);
rho_MJLS=zeros(no_samples,1);
for i=1:no_samples
    %% Define the MJLS
    [MJLS,Ber]=example_10();      
    %% Analyze the defined MJLS
    [B_cal_MJLS,T_cal_MJLS]=get_B_T_matrices(MJLS);
    rho_B_MJLS=max(abs(eig(B_cal_MJLS)));
    rho_MJLS(i)=max(abs(eig(T_cal_MJLS)));

    %% Analyze the defined Bernoulli system
    [B_cal_Ber,T_cal_Ber]=get_B_T_matrices(Ber);
    rho_B_Ber=max(abs(eig(B_cal_Ber)));
    rho_Ber(i)=max(abs(eig(T_cal_Ber)));
end
%% Find examples where Bernoulli is stable but MJLS isn't
id1=find(rho_Ber<1); % Bernouli stable
id2=find(rho_MJLS<1); % MJLS stable

%% Plots
figure()
plot(1:no_samples,rho_MJLS,'r')
hold on
plot(1:no_samples,rho_Ber,'b')
plot(1:no_samples,1*ones(1,no_samples),'k--')
legend('MJLS spectral radius','Ber spectral radius')
title('All examples')

% Plot examples which are Bernoulli stable
figure()
plot(1:size(id1,1),rho_MJLS(id1),'r')
hold on
plot(1:size(id1,1),rho_Ber(id1),'b')
plot(1:size(id1,1),1*ones(1,size(id1,1)),'k--')
legend('MJLS spectral radius','Bernoulli spectral radius')
title('Examples where the Bernouli model is stable')

% Plot examples which are MJLS stable
figure()
plot(1:size(id2,1),rho_MJLS(id2),'r')
hold on
plot(1:size(id2,1),rho_Ber(id2),'b')
plot(1:size(id2,1),1*ones(1,size(id2,1)),'k--')
legend('MJLS spectral radius','Bernoulli spectral radius')
title('Examples where the MJLS model is stable')


%% Helper functions
function [MJLS,Ber]=example_10()
% Create an MJLS system: Example 10: scalar example to compare Bernoulli vs
% Markov modeling

MJLS.N=2;
A1=2*rand();
A2=2*rand();        
MJLS.nx=size(A1,1);
MJLS.As={A1,A2};
%MJLS.T=[0.85,0.6;0.15,0.4];
MJLS.T=rand(2);
MJLS.T=MJLS.T*inv(diag(ones(1,2)*MJLS.T));
MJLS.x_ic=1;
MJLS.p_ic=[0.7;0.3];


% Find the stationary distribution
[eig_vecs,~] = eigs(MJLS.T);
p_eqm=eig_vecs(:,1)/sum(eig_vecs(:,1)); % Normalize it so that the sum is 1

% Return a Bernoulli model with the stationary distribution
Ber=MJLS;
Ber.T=p_eqm*ones(1,size(p_eqm,1),1);
end