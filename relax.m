function current_approximation = relax(current_approximation, A, rhs, num_sweeps_down)
% [A, rhs] = matrix_assemble(2,sqrt(max(size(rhs))),source);
current_approximation = gauss_seidel(A,rhs,current_approximation,num_sweeps_down);
end