function rhs = restrict(residual)
fine_num_1D = sqrt(max(size(residual)));
coarse_num_1D = (fine_num_1D - 1)/2 + 1;
rhs = zeros(coarse_num_1D^2,1);
fine_mesh = reshape(residual,[fine_num_1D,fine_num_1D]);
coarse_mesh = reshape(rhs,[coarse_num_1D,coarse_num_1D]);

for i = 1: coarse_num_1D
    for j = 1: coarse_num_1D
        if i > 1 && i <coarse_num_1D && j > 1 && j <coarse_num_1D
            coarse_mesh(i,j) = (fine_mesh(2*(i-1),2*(j-1)) +fine_mesh(2*(i-1),2*j)+fine_mesh(2*i,2*(j-1))+...
                fine_mesh(2*i,2*j)+ 2*(fine_mesh(2*i-1,2*(j-1))+fine_mesh(2*i-1,2*j)+ ...
                fine_mesh(2*(i-1),2*j-1)+fine_mesh(2*i,2*j-1)+4*fine_mesh(2*i-1,2*j-1)))/16;
        else
            coarse_mesh(i,j) = fine_mesh(2*i-1,2*j-1);
        end
    end
end
rhs = reshape(coarse_mesh,[coarse_num_1D^2,1]);
end