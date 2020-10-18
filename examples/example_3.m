function MJLS=example_3()
% Create an MJLS system: Example 1: Unstable+Unstable=Stable
MJLS=struct;

MJLS.N=2;
A1=[2,-1;0.1,0.1]; 
A2=[0.2,1;-0.1,2];            

MJLS.nx=size(A1,1);
MJLS.As={A1,A2};
MJLS.T=[0.1,0.9;0.9,0.1];  

MJLS.x_ic=[-10;10];
MJLS.p_ic=[0.5;0.5];
end