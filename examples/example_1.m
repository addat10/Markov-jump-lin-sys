function MJLS=example_1()
% Create an MJLS system: Example 3: scalar example not MSS
MJLS.N=2;
A1=1.15;
A2=0.1;        
MJLS.nx=size(A1,1);
MJLS.As={A1,A2};
MJLS.T=[0.85,0.6;0.15,0.4];
MJLS.x_ic=1;
MJLS.p_ic=[0.7;0.3];
end