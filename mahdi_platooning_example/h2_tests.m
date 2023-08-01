clear
clc
close all
%%
% eps=0.1;
lambda=1:10;
eps=linspace(0.0001,2/lambda(end)-0.001,1000);

for i=1:length(eps)
    A=diag(1-eps(i)*lambda);
    B=eye(length(lambda));
    C=diag(lambda);
    G=ss(A,B,C,0,1); 
    h2(i)=norm(G,2);
    hinf(i)=norm(G,inf);
    spec_rad(i)=max(abs(eig(A)));
end

[hinf_min,id_hinf]=min(hinf);
[h2_min,id_h2]=min(h2);
[spec_rad_min,id_spec_rad]=min(spec_rad);


eps_opt_h2=eps(id_hinf)
eps_opt_hinf=eps(id_h2)
eps_opt_spec_rad=eps(id_spec_rad)



figure()
plot(eps,h2)

figure()
plot(eps,hinf)

figure()
plot(eps,spec_rad)
