% Test the space of distributions
clear all
close all
clc
%%
n=4; 
[U,S,V]=svd(ones(n,1));
corners=U(:,2:end)'; % Corners of the tetrahedron (shadow of set of distributions in R4)
%%
p=0:0.01:1; % Success probability in edge 1
q=0:0.01:1; % Success probability in edge 2
p_vec=[];
p_test=[];
for i=1:length(p)        
    for j=1:length(q)
        p_vec=[p_vec,[p(i)*q(j);p(i)*(1-q(j));(1-p(i))*q(j);(1-p(i))*(1-q(j))]]; % collect set of all Bernoulli models                
        if (i==1||i==length(p))&& ((j==1||j==length(q)))
            p_test=[p_test,[[p(i)*q(j);p(i)*(1-q(j));(1-p(i))*q(j);(1-p(i))*(1-q(j))]]];
        end
    end
end
p_proj=U(:,2:end)'*p_vec; 
p_test_proj=U(:,2:end)'*p_test;
%% Visualize the set of bernoulli distributions

figure()
plot3(p_proj(1,:),p_proj(2,:),p_proj(3,:),'r.')
hold on
plot3(p_test_proj(1,:),p_test_proj(2,:),p_test_proj(3,:),'b.')
hold on
draw_lines(corners(:,1),corners(:,2))
draw_lines(corners(:,1),corners(:,3))
draw_lines(corners(:,1),corners(:,4))
draw_lines(corners(:,2),corners(:,3))
draw_lines(corners(:,2),corners(:,4))
draw_lines(corners(:,3),corners(:,4))

% figure()
% [k2,av2] = convhull(p_proj(1,:),p_proj(2,:),p_proj(3,:),'Simplify',true);
% 
% trisurf(k2,p_proj(1,:),p_proj(2,:),p_proj(3,:),'FaceColor','cyan')
% axis equal
% hold on
% draw_lines(corners(:,1),corners(:,2))
% draw_lines(corners(:,1),corners(:,3))
% draw_lines(corners(:,1),corners(:,4))
% draw_lines(corners(:,2),corners(:,3))
% draw_lines(corners(:,2),corners(:,4))
% draw_lines(corners(:,3),corners(:,4))

figure()
[k2,av2] = convhull(p_test_proj(1,:),p_test_proj(2,:),p_test_proj(3,:),'Simplify',true);

trisurf(k2,p_test_proj(1,:),p_test_proj(2,:),p_test_proj(3,:),'FaceColor','blue')
axis equal
hold on
draw_lines(corners(:,1),corners(:,2))
draw_lines(corners(:,1),corners(:,3))
draw_lines(corners(:,1),corners(:,4))
draw_lines(corners(:,2),corners(:,3))
draw_lines(corners(:,2),corners(:,4))
draw_lines(corners(:,3),corners(:,4))
%% Function
function [] = draw_lines(point1,point2)
        xyz=[point1';point2'];
        line(xyz(:,1),xyz(:,2),xyz(:,3))
end