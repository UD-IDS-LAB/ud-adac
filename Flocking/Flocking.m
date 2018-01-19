%A flocking simulation script based on the paper "Flocking in Fixed and
%Switching Networks" by Tanner et al, IEEE Transations on Automatic
%Control, 2007

clear; close all; clc;
%simulation of flocking behaviors based on Tanner's paer
R = 5;% sensing radius, [m]
Rhat = 1; %desired spacing, [m]
n = 10; %number of agents
%random x, y, heading
position = rand(n, 2);
velocity = rand(n, 2);
%scaling position vector to (-1, 1)
position = (position-0.5)*2; %x
velocity = (velocity-0.5)*2; %y
%scale position and velocity
velocity = velocity * 1.5;
position = position * 5;
%final velocity?
vf = mean(velocity);
vf = vf ./ sqrt(sum(vf.^2));
%plot initial position
PlotFlock(position, velocity);
saveas(gcf, 'initialState.png');
movegui(gca, 'northwest');
[G, A] = Adjacency(position, R);
saveas(gcf, 'initialGraph.png');
movegui(gca, 'northeast');
%%
%run the RK4 solver for the dynamic system
dt = 0.05; %seconds
tf = 5;
r_save = position;
v_save = velocity;
for t = 0:dt:tf
    %values used for velocity matching
    dv = zeros(n, 2); %velocity difference
    for i = 1:n
        for j = i+1:n
            dv(i,:) = dv(i,:) + (velocity(j,:) - velocity(i,:))*A(i,j);
            dv(j,:) = dv(j,:) - (velocity(j,:) - velocity(i,:))*A(i,j);
        end
    end
    %potential function values
    dr = 1e-3; %gradient numeric step size
    gv = zeros(n, 2); %gradient potential
  %  fprintf('\n');
     for i = 1:n
        for j = i+1:n
             rij= sqrt((position(i,1) - position(j,1))^2 + ...
                  (position(i,2) - position(j,2))^2 );
        %    fprintf('r%g,%g: %g ', i,j, rij);
            %evaluate potential function and derivative
            [V, dvdr] = Potential(rij, Rhat, R);
            %calculate partial dr/dx for chain rule
            drdx = (position(i,1) - position(j,1))/rij;
            drdy = (position(i,2) - position(j,2))/rij;
            %update potential gradient sum
            gv(i,1) = gv(i,1) - dvdr*drdx;
            gv(j,1) = gv(j,1) - dvdr*drdy;
            gv(i,2) = gv(i,2) + dvdr*drdx;
            gv(j,2) = gv(j,2) + dvdr*drdy;
        end
     end
    %euler's method for double integrator dynamics
    acceleration = zeros(n, 2) + dv + gv;
    velocity = velocity + acceleration * dt;
    position = position + velocity*dt;
    %save velocity and position
    r_save = [r_save, position];
    v_save = [v_save, velocity];
    %plot position and connectivity graphs
    PlotFlock(position, velocity);
    [G, A] = Adjacency(position, R);
    pause(0.05);
end
%%
vfa = mean(velocity);
vfa = vfa ./ sqrt(sum(vfa.^2));

%%
nt = length(r_save);
figure(3); clf; hold on;
lw = 2;
for i = 1:n %for each agent
    c = rand(1, 3);
    plot(r_save(i, 1:2:nt), r_save(i, 2:2:nt), '--', 'color', c, 'linewidth', lw);
    plot(r_save(i, end-1), r_save(i, end), 'o', 'color', c, 'linewidth', lw);
    quiver(position(i,1), position(i,2), velocity(i,1), velocity(i,2), 0, ...
                                            'color', c, 'linewidth', lw/2);
end
xlabel('X Position [m]');
ylabel('Y Position [m]');
grid on;
box on;
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on')
hold off;
saveas(gcf, 'finalState.png');

[G, A] = Adjacency(position, R);
saveas(gcf, 'finalGraph.png');
%%
figure(10); clf;
movie = VideoWriter('flocking.avi'); % Name it.
vidTime = 10; %video time in seconds
movie.FrameRate = nt/2/vidTime;
open(movie);
for s = 1:nt/2
    clf; hold on;

    
    plot(r_save(:, 2*s-1), r_save(:, 2*s), 'ok');
    quiver(r_save(:, 2*s-1), r_save(:, 2*s),v_save(:, 2*s-1), v_save(:, 2*s), 0);
    
    mip = min(r_save);
    map = max(r_save);
    axis([min(mip(1:2:end))-1 max(map(1:2:end))+1 ...
          min(mip(2:2:end))-1 max(map(2:2:end))+1]);
    xlabel('X Position [m]');
    ylabel('Y Position [m]');
    grid on;
    box on;
    set(gca,'YMinorTick','on')
    set(gca,'XMinorTick','on')
    pause(0.01);
    writeVideo(movie, getframe(gcf));
end
close(movie);







