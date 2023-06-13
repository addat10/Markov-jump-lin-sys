clear
clc
%%
data_1=load('2458_Base_200_300');
eigs_1=sort(abs(eig(data_1.MJLS.T)),'descend');
N=size(data_1.MJLS.T,1);

data_2=load('2458_CC_200_300.mat');
eigs_2=sort(abs(eig(data_2.MJLS.T)),'descend');

figure()
plot(1:N,eigs_1,'b')
hold on
plot(1:N,eigs_2,'r')
legend('Base','CC')