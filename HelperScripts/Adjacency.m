%this is a helper function that calculates an adjacency matrix and builds a
%graph based on an array of positions and a sensing radius
function [ G, A ] = Adjacency(position, R)
s = size(position);
n = s(1);
A = zeros(n);
for i = 1:n
    for j = i+1:n
        rij= sqrt((position(i,1) - position(j,1))^2 + ...
                  (position(i,2) - position(j,2))^2 );
        A(i,j) = (rij <= R); %upper right half
        A(j,i) = A(i,j); %symmetric
    end
end
G = graph(A);
%plot
figure(2); clf; hold on;
plot(G); box on; hold off;
end

