function Grid = grid_generator(h_finest,gridnum)
    field1 = 'f';  value1 = zeros(2,1);
    field2 = 'v';  value2 = zeros(2,1);
    field3 = 'u';  value3 = zeros(2,1);
    grid = struct(field1,value1,field2,value2,field3,value3);
    Grid = repmat(grid,1,gridnum);
    interval_num_finest = 1/h_finest;
    for i = 1: gridnum
        Grid(i).f = zeros((interval_num_finest/ 2^(i-1) + 1)^2,1);
        Grid(i).v = Grid(i).f;
        Grid(i).u = Grid(i).f;
    end
end