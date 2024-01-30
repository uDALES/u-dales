function TR_out = operateUnitCube(scale, shift, divisions)
    % Shift defines where the centre of the cubes are
    % Scale defines how large they are

    TR_in = stlread('cube.stl'); 
    % this cube is centred on the origin and has side length 1
    % scaling is done first, then shifting
    % So for a regular array of cube size(s) size Hx, Hy, Hz
    % scale = [Hx Hy Hz], shift(n,:) = [2*n*Hx 2*n*Hy Hz/2];

    TR = TR_in;
    for n=1:divisions
        TR = divideFaces(TR);
    end

    cube1_vert = TR.Points;
    cube1_conn = TR.ConnectivityList;

    ncubes = size(shift,1);
    cubes_vert = [];
    cubes_conn = [];
    for n=1:ncubes
        cubes_vert = [cubes_vert; cube1_vert.*scale + shift(n,:)];
        cubes_conn = [cubes_conn; cube1_conn + size(TR.Points,1) * (n-1)];
    end

    TR_out = triangulation(cubes_conn, cubes_vert);

end