function [A,b,u] = matrix_assemble(dimension, h, source, solution)
num_of_point = 1/h +1;
N = num_of_point^dimension;
b = zeros(N,1);
u = zeros(N,1);
A = zeros(N,N);
x = linspace(0,1,num_of_point);
y = linspace(0,1,num_of_point);
[X,Y] = meshgrid(x,y);
is_boundary = (X == 1) + (X == 0) + (Y == 1) + (Y == 0);
boundary = find(is_boundary);
internal = find(~is_boundary);
Dxx = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,num_of_point));
Dyy = toeplitz(sparse([1,1],[1,2],[-2,1]/h^2,1,num_of_point));
%     Dxx = toeplitz([-2/h^2, 1/h^2, zeros(1, num_of_point-2)]); %no sparse version
%     Dyy = toeplitz([-2/h^2, 1/h^2, zeros(1, num_of_point-2)]);
b(internal) = source(X(internal),Y(internal));
u(internal) = solution(X(internal),Y(internal));
% A = kron(Dxx,eye(num_of_point))+kron(eye(num_of_point), Dyy);
A = kron(Dxx, speye(num_of_point)) + kron(speye(num_of_point), Dyy);
for i = 1:max(size(boundary))
    b(boundary(i)) = 0;
    u(boundary(i)) = 0;
    A(boundary(i),:) = 0;
    A(boundary(i),boundary(i)) = 1;
end
A = sparse(A);
b = sparse(b);
end