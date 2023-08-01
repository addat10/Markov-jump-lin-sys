clear
clc
close all
%%
lam=1:10;
n=length(lam);
N=100;
eps_max=2/max(lam);
eps=eps_max/N:eps_max/N:eps_max;

for t=1:length(eps)
    g(t)=0;
    for i=1:n
        f(i)=lam(i)/(eps(t)*(2-lam(i)*eps(t)));
        g(t)=g(t)+f(i);
    end
    
end
figure()
plot(eps,g,'*')
g_min=min(g);
eps_star=eps(find(g==g_min))
1/max(lam)