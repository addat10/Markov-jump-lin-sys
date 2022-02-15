clear all
close all
clc
%% main
addpath('../lib')
addpath(genpath('./connectivity_data_mahdi'))
dir1='2458.CC';% 888.Base
dir2='0-100';

%TPM_path='./connectivity_data_mahdi/Transition_Probability.mat'; % First data
TPM_path=['./connectivity_data_mahdi/MC_probability/',dir1,'/',dir2,'/Transition_Probability.mat']; % First data

%% Define the MJLS
Kp=0.45:0.01:0.7;
T=1;
rho_T_MJLS=zeros(length(Kp),1);
rho_T_Ber=zeros(length(Kp),1);
lead_flag=0;% Set to zero for consensus problem and 1 if leader present
for i=1:length(Kp)
[MJLS,Ber]=example_ber_vs_markov_first_order_consensus_mahdi_data(TPM_path,Kp(i),T,lead_flag);  
% Analyze the defined MJLS
[B_cal_MJLS,T_cal_MJLS]=get_B_T_matrices(MJLS);
[B_cal_Ber,T_cal_Ber]=get_B_T_matrices(Ber);
% rho_B_MJLS=max(abs(eig(B_cal_MJLS)));
% rho_B_Ber=max(abs(eig(B_cal_Ber)));
rho_T_MJLS(i)=max(abs(eig(T_cal_MJLS)));
rho_T_Ber(i)=max(abs(eig(T_cal_Ber)));
end

%% Plots
figure()
plot(Kp*T,rho_T_Ber,'r')
hold on
plot(Kp*T,rho_T_MJLS,'b')
legend('Bernoulli','Markov')


%% Analyze probability dynamics
eigs=eig(MJLS.T);
figure
plot(eigs,'o')
title('Spectrum of the transition Matrix')


%% helper functions
function [MJLS,Ber]=example_ber_vs_markov_first_order_consensus_mahdi_data(TPM_path,Kp,T,lead_flag)
% Create an MJLS system with data from Mahdi
MJLS.N=64;
% Need to create 64 A matrices
MJLS.As={};
MJLS.Bs={};
MJLS.Laps={};
MJLS.dt=T;
%Kp=3;
n=5;
At=diag(ones(n,1))+diag(ones(n-1,1),1); % Always connected
for i=1:MJLS.N            
    At=update_adjacency(At,i);
    D=diag(At'*ones(n));
    if lead_flag==1
        % Don't need the leader for the consensus problem
        D(1)=2;
    end
    Li=diag(D)-At';
    MJLS.Laps{i}=Li;
    [Ai,Bi]=get_AB(Li,Kp,MJLS.dt,lead_flag);
    MJLS.As{i}=Ai;
    MJLS.Bs{i}=Bi;
end
MJLS.nx=size(Ai,1);
MJLS.x_ic=0*rand(n,1);
load(TPM_path,'Transition_Probability');
MJLS.T=Transition_Probability;



% Find the stationary distribution
[eig_vecs,~] = eigs(Transition_Probability);
p_eqm=eig_vecs(:,1)/sum(eig_vecs(:,1)); % Normalize it so that the sum is 1
%Genetare a random initial distribution but is required only for simulation
p_ic=rand(64,1);
MJLS.p_ic=1/sum(p_ic)*p_ic;

Ber=MJLS;
Ber.T=p_eqm*ones(1,size(p_eqm,1),1);
end
function At=update_adjacency(At,i)
    vec=dec2bin(i-1,6);
    At(1,3)=bin2dec(vec(1));
    At(1,4)=bin2dec(vec(2));
    At(1,5)=bin2dec(vec(3));
    At(2,4)=bin2dec(vec(4));
    At(2,5)=bin2dec(vec(5));
    At(3,5)=bin2dec(vec(6));
end
function [A,B]=get_AB(L,Kp,T,lead_flag)
n=size(L,1);
if lead_flag==0
    % Need the projection for the consensus system (marginally stable)
    Proj=(eye(n)-ones(n)/n);
    A=Proj*(eye(n)-Kp*T*L);
else
    A=eye(n)-Kp*T*L;
end
B=T*Kp*L;
end
