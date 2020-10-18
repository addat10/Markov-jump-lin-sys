clear all
close all
clc 
%% Initialize parameters
x=[1;2];
% non-symmetric, convergent towards high probab of theta=1
T=[0.95,0.9;0.05,0.1];
% Symmetric, convergent towards average (Recall consensus)
%T=[0.7,0.3;0.3,0.7];
% Non-convergent(symmetric) towards average
%T=[0,1;1,0];
% Another symmetric Example
T=[0.1,0.9;0.9,0.1];

p_ic=[0.3;0.7];
steps=10;
%% Theoretically investigations
[V,lambda]=eig(T);
% Scale the eig vector corresponding to lambad=1 so that sum of its entries
% is 1
pi_vec=(1/(V(:,1)'*[1;1]))*V(:,1);

%% Propogate probabilities and plot
p=zeros(2,steps);
p(:,1)=p_ic;
figure()
plot(x,p(:,1))
xlim([0.5,2.5])
ylim([-0.1,1.1])
ylabel('probability')
xlabel('possible values of theta are 1 and 2')
title('evolution of probabilities')
hold on
pause
for i=2:steps
    p(:,i)=T*p(:,i-1);
    plot(x,p(:,i))
    pause
end

