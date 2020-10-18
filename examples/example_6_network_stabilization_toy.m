function MJLS=example_6_network_stabilization_toy()
% Create an MJLS system: For Example
% MJLS.N=2;
% A1=[-0.5,2;-0.5,0.5];
% A2=[0.5,0;1,0]; 
% MJLS.nx=size(A1,1);
% MJLS.As={A1,A2};
% MJLS.T=[0.5,0.5;0.5,0.5]; 
% MJLS.x_ic=[1;1];
% MJLS.p_ic=[0.5;0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
no_agents=5;        % number of agents
density=[0.8,0.1];  % Probability that a link is present between agents i and j
N=2;                % Number of Markov modes (possible graphs)
k=10;               % (Spring stiffness)gain between informed agent and origin
T=0.1;              % Sampling time
rng(1)              % Seed for the random number generator
for i=1:N    
    L=gen_rand_Laplacian(no_agents,density(i));    
    if i==1
        while min(L(:,1)==zeros(no_agents,1))
            L=gen_rand_Laplacian(no_agents,density(i)); % atleast one agent connected to informed agent
        end
    else
        while max(L(:,1)~=zeros(no_agents,1))
            L=gen_rand_Laplacian(no_agents,density(i)); % atleast one agent connected to informed agent
        end        
    end    
    L(1,:)=zeros(1,no_agents);  % Informed agent unaffected by others
    L(1,1)=k*T;                 % Connect informed agent to origin
    A{i}=eye(no_agents)-T*L;    
end




MJLS.N=N;
MJLS.nx=size(A{1},1);
MJLS.As=A;
T1=[0.01;0.99];
MJLS.T=T1*ones(1,2); 
MJLS.x_ic=rand(no_agents,1);
MJLS.p_ic=(1/N)*ones(N,1);
end