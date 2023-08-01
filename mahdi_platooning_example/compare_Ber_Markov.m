clear
clc
%%
data_1=load('2458_Base_200_300');
eigs_1=sort(abs(eig(data_1.MJLS.T)),'descend');
N=size(data_1.MJLS.T,1);

data_2=load('2458_CC_200_300.mat');
eigs_2=sort(abs(eig(data_2.MJLS.T)),'descend');

MJLS=data_1.MJLS;
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

