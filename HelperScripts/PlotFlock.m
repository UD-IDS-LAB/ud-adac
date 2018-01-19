% this is a helper function that plots a flock of agents, using the
% position to draw a point and the velocity to draw an arrow
function PlotFlock(position, velocity)
figure(1); clf; hold on;
plot(position(:,1), position(:,2), 'ok');
quiver(position(:,1), position(:,2), ...
       velocity(:,1), velocity(:,2), 0);

s = size(position);
for i = 1:s(1)
    text(position(i,1) + 0.05, position(i, 2), num2str(i))
end

grid on;
box on;
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on')
hold off;
end