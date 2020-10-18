clear all
close all
clc 
%% Initialize parameters
n=10;
x=(1:n)';
steps=30;
%p_ic=(1/n)*ones(n,1);
p_ic=zeros(n,1);
p_ic(floor(n/2))=1;


% Non-convergent(symmetric) towards average
density=0.4;
T=sprand(n,n,density);
T=T*inv(diag(ones(1,n)*T));
%% Theoretically investigations
[V,lambda]=eig(T);
% Scale the eig vector corresponding to lambad=1 so that sum of its entries
% is 1
pi_vec=(1/(V(:,1)'*ones(n,1)))*V(:,1);

%% Propogate probabilities and plot
p=zeros(n,steps);
p(:,1)=p_ic;
figure()
plot(x,p(:,1))
hold on
for i=2:steps
    %% 
    p(:,i)=T*p(:,i-1);    
    plot(x,p(:,i))
    pause    
end
xlim([0.5,n+0.5])
ylim([-0.1,1.1])
ylabel('probability')
xlabel('possible values of theta are 1 and 2')
title('evolution of probabilities')

