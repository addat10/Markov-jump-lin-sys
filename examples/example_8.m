function MJLS=example_8()
% Create an MJLS system: mass spring damper system with rand var k
T=0.1;
k_1=-10;
k_2=1;

MJLS.N=2;

A1=[1,T;-k_1*T,1-T];
A2=[1,T;-k_2*T,1-T];
p1=0.25;
p2=0.75;

MJLS.nx=size(A1,1);
MJLS.As={A1,A2};
MJLS.T=[p1,p1;p2,p2]; 
MJLS.x_ic=[1;1];
MJLS.p_ic=[p1;p2];
end