% A simulation that uses consensus variables to bring all agents to a
% predefined displacement. Based on the paper "Consensus and Cooperation in
% Networked Multi-Agent Systems" by Olfati-Saber et al in the proceedings
% of the IEEE, 2007

%position based consensus analysis
clear; close all; clc;

%simulation parameters
n = 10;   % number of agents
s = 10;   % initial region size
R = 4;    % communication range
%control variables
r = 3;    % desired relative position
k = 1;    % proportionality constant
%initial state
pos = (rand(n,2) - 0.5)*s;
vel = zeros(size(pos));
%plot initial state
PlotFlock(pos, vel);

%control loop
dt = 0.01; ti = 0; tf = 100*dt;
figure(2); clf;
movegui('northeast');
[G, A] = Adjacency(pos, R);
for t = ti:dt:tf
    v = zeros(size(pos));
    for i = 1:n   %for each agent
        b = [0 0];
        v = [0 0];
        for j = 1:n %for every neighbor
            if i ~= j && A(i,j)
                dp = (pos(j,:) - pos(i,:));
                dv = (vel(j,:) - vel(i,:));
                v = v + k*dv;
                b = b - r*dp/sqrt(sum(dp.^2));
            end
        end
        vel(i,:) = vel(i,:) + (v+b)*dt;
    end
    
    pos = pos + vel.*dt;
    
    [G, A] = Adjacency(pos, R);
    PlotFlock(pos, vel);
    pause(0.001);
    
end
