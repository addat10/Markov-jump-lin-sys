% Test the space of distributions
clear all
close all
clc
%%
n=4; % Choose n=2^m
[U,S,V]=svd(ones(n,1));
corners=U(:,2:end)';
%%
p=0:0.01:1;
p_vec=[];
for i=1:length(p)    
    switch n
        case 2
            p_vec=[p_vec,[p(i);1-p(i)]]  ;      
        case 4
            p_vec=[p_vec,[p(i)^2;p(i)*(1-p(i));(1-p(i))*p(i);(1-p(i))^2]];
    end 
end
p_proj=U(:,2:end)'*p_vec;    
%% Visualize the set of bernoulli distributions

switch n
    case 2
        figure()
        plot(p_proj,zeros(length(p_proj),1),'*')
        hold on
        plot(corners(1),0,'ro')
        plot(corners(2),0,'ro')
    case 4
        figure()
        plot3(p_proj(1,:),p_proj(2,:),p_proj(3,:),'r*')
        hold on
        draw_lines(corners(:,1),corners(:,2))
        draw_lines(corners(:,1),corners(:,3))
        draw_lines(corners(:,1),corners(:,4))
        draw_lines(corners(:,2),corners(:,3))
        draw_lines(corners(:,2),corners(:,4))
        draw_lines(corners(:,3),corners(:,4))

end
%% Functions
function [] = draw_lines(point1,point2)
        xyz=[point1';point2'];
        line(xyz(:,1),xyz(:,2),xyz(:,3))
end