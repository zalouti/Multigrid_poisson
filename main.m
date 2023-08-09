
% define parameters and equation
h_coarsest = 1/4;
h_finest = 1/64;
finest_num_1D = 1/h_finest +1;
source = @(x,y) -2.*pi.^2.*sin(pi.*x).*sin(pi.*y);
solution = @(x,y) sin(pi.*x).*sin(pi.*y);
iteration1 = 10;
iteration2 = 10;
x = zeros(finest_num_1D^2,1);
err = zeros(1,5);

% solve and plot
for i =1:10
    [error, x] = multigrid(h_coarsest,h_finest, iteration1, iteration2, source, solution, x);
    err(i) = error;
end
lgerr = log10(err);
k = (lgerr(10) - lgerr(1))/(10-1);
speed = 10^(abs(k));
plot(1:1:10,err,'LineWidth',2);
title("V-cycle convergence speed")
xlabel("iteration number")
ylabel("infty norm")

% multigrid implementation
function [error, numeric_solution] = multigrid(h_coarsest,h_finest, iteration1, iteration2, source, solution, initialguess)
    gridnum = log2(h_coarsest/h_finest) + 1;
    Grid = grid_generator(h_finest,gridnum);
    for i = 1 : gridnum - 1
        if i == 1
            [A, Grid(i).f, Grid(i).u] = matrix_assemble(2,h_finest * 2^(i-1),source,solution);
            Grid(i).v = initialguess;
        else
            [A, ~, ~] = matrix_assemble(2,h_finest * 2^(i-1),source,solution);
        end
        Grid(i).v = relax(Grid(i).v, A, Grid(i).f,iteration1);
        Grid(i+1).f = restrict(Grid(i).f - A*Grid(i).v);
    end
    [A, ~, ~] = matrix_assemble(2,h_coarsest,source,solution);
    Grid(gridnum).v = gauss_seidel(A,Grid(gridnum).f,Grid(gridnum).v,iteration1);
    for j = gridnum - 1:-1:1
        [A,~,~] = matrix_assemble(2,h_finest * 2^(j-1),source,solution);
        Grid(j).v = Grid(j).v + grid_interpolate(Grid(j+1).v);
        Grid(j).v = relax(Grid(j).v, A, Grid(j).f, iteration2);
    end
    error = norm(Grid(1).v-Grid(1).u,inf);
    numeric_solution = Grid(1).v;
    end