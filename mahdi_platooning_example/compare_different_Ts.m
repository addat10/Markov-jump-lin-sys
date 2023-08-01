clear all
close all
clc
addpath('../lib')
addpath(genpath('./connectivity_data_mahdi'))
dir1='888.Base';% 888.Base
dir2='100-200';

%TPM_path='./connectivity_data_mahdi/Transition_Probability.mat'; % First data
TPM_path=['./connectivity_data_mahdi/MC_probability/',dir1,'/',dir2,'/Transition_Probability.mat']; % First data
%% Define the MJLS
TPM_path=['./connectivity_data_mahdi/MC_probability/888.Base/100-200/Transition_Probability.mat']; % First data
MJLS1=example_single_integrator_mahdi_data(TPM_path);
T1=MJLS1.T;
[B_cal1,T_cal1]=get_B_T_matrices(MJLS1);
rho_B1=max(abs(eig(B_cal1)));
rho_T1=max(abs(eig(T_cal1)));
[U1,eigs1]=eig(T_cal1);
figure()
plot(eigs1,'o')
title('Spectrum of the transition Matrix')


TPM_path=['./connectivity_data_mahdi/MC_probability/888.Base/200-300/Transition_Probability.mat']; % First data
MJLS2=example_single_integrator_mahdi_data(TPM_path);
T2=MJLS2.T;
[B_cal2,T_cal2]=get_B_T_matrices(MJLS2);
rho_B2=max(abs(eig(B_cal2)));
rho_T2=max(abs(eig(T_cal2)));
[U2,eigs2]=eig(T_cal2);
figure()
plot(eigs2,'o')
title('Spectrum of the transition Matrix')


TPM_path=['./connectivity_data_mahdi/MC_probability/888.Base/300-400/Transition_Probability.mat']; % First data
MJLS3=example_single_integrator_mahdi_data(TPM_path);
T3=MJLS3.T;
[B_cal3,T_cal3]=get_B_T_matrices(MJLS3);
rho_B3=max(abs(eig(B_cal3)));
rho_T3=max(abs(eig(T_cal3)));
[U3,eigs3]=eig(T_cal3);
figure()
plot(eigs3,'o')
title('Spectrum of the transition Matrix')


% p_ic happens to be the stationary distribution
%% Analyze the defined MJLS