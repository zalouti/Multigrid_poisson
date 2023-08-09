function vector_fine = grid_interpolate(vector_coarse)
coarse_num_1D = sqrt(max(size(vector_coarse)));
fine_num_1D = coarse_num_1D*2 - 1;
vector_fine = zeros(fine_num_1D^2,1);
fine_mesh = reshape(vector_fine,[fine_num_1D,fine_num_1D]);
coarse_mesh = reshape(vector_coarse,[coarse_num_1D,coarse_num_1D]);

for i = 1: coarse_num_1D
    for j = 1: coarse_num_1D
        if i >= 1 && i < coarse_num_1D && j >= 1 && j < coarse_num_1D
            fine_mesh(2*i-1,2*j-1) =  coarse_mesh(i,j);
            fine_mesh(2*i,2*j-1) = (coarse_mesh(i,j)+coarse_mesh(i+1,j))/2;
            fine_mesh(2*i-1,2*j) = (coarse_mesh(i,j)+coarse_mesh(i,j+1))/2;
            fine_mesh(2*i,2*j) = (coarse_mesh(i,j)+coarse_mesh(i+1,j)+coarse_mesh(i,j+1)+coarse_mesh(i+1,j+1))/4;
        else
            fine_mesh(2*i-1,2*j-1) =  coarse_mesh(i,j);
            if i == coarse_num_1D && j ~= coarse_num_1D
                fine_mesh(2*i-1,2*j) = (coarse_mesh(i,j)+coarse_mesh(i,j+1))/2;
            elseif i ~= coarse_num_1D && j == coarse_num_1D
                fine_mesh(2*i,2*j-1) = (coarse_mesh(i,j)+coarse_mesh(i+1,j))/2;
            end
        end
    end
end
vector_fine = reshape(fine_mesh,[fine_num_1D^2,1]);
end