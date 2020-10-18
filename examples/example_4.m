function MJLS=example_4()
% Create an MJLS system: Example 2: Stable+Stable=Unstable
MJLS.N=2;
A1=[-0.5,2;-0.5,0.5];
A2=[0.5,0;1,0];

MJLS.nx=size(A1,1);
MJLS.As={A1,A2};
MJLS.T=[0.5,0.5;0.5,0.5]; 
MJLS.x_ic=[1;1];
MJLS.p_ic=[0.5;0.5];
end