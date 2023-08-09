function x = gauss_seidel(A,b,x,step)
% solve equation Ax= b using Gauss Seidel method
% parameters:
% A: coefficient matrix
% b: coefficient vector
% x: solution to be solved
% step: iteration steps
    N = max(size(b));
    [B,C] = deal(zeros(N,N));
    C = inv(tril(A));
    B = - C * triu(A,1);
    f = tril(A)\b;
    for i = 1:step
        x_new = B * x + f;
        if x_new == x
            break;
        else
            x = x_new;
        end
    end
end