h_coarsest = 1/4;
h_finest = 1/64;
finest_num_1D = 1/h_finest +1;
source = @(x,y) -2.*pi.^2.*sin(pi.*x).*sin(pi.*y);
solution = @(x,y) sin(pi.*x).*sin(pi.*y);
iteration1 = 10;
iteration2 = 10;
x = zeros(finest_num_1D^2,1);
err = zeros(1,5);
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
