function MJLS=example_2()
% Create an MJLS system: Example 3: scalar example MSS
MJLS.N=2;
A1=1.1;
A2=0.1;        
MJLS.nx=size(A1,1);
MJLS.As={A1,A2};
MJLS.T=[0.75,0.6;0.25,0.4];
MJLS.x_ic=1;
MJLS.p_ic=[0.7;0.3];
end