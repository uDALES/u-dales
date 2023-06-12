function [fluid_IB, solid_IB] = ...
    getBoundaryCells(xgrid, ygrid, zgrid, fluid, solid, include_diagonals)

itot = length(xgrid);
jtot = length(ygrid);
ktot = length(zgrid);

% % check
nfluid = sum(fluid, 'all');
nsolid = sum(solid, 'all');
% assert(nsolid == length(xgrid)*length(ygrid)*length(zgrid) - nfluid)

%fluid_xyz = zeros(nfluid,3);
%solid_xyz = zeros(nsolid,3);

fluid_IB = false(length(xgrid), length(ygrid), length(zgrid));
solid_IB = false(length(xgrid), length(ygrid), length(zgrid));
bound = false(length(xgrid), length(ygrid), length(zgrid));

fluid_IB_xyz = [];
solid_IB_xyz = [];
bound_xyz = [];

nfluid_count = 0;
nsolid_count = 0;

for i=1:itot
    for j=1:jtot
        for k=1:ktot
            if fluid(i,j,k)
                nfluid_count = nfluid_count + 1;
                %fluid_xyz(nfluid_count,:) = [xgrid(i), ygrid(j), zgrid(k)];

                % Identify fluid IB points
                if i~=1
                    if solid(i-1,j,k)
                        fluid_IB(i,j,k) = true;
                    end
                end

                if i~=itot
                    if solid(i+1,j,k)
                        fluid_IB(i,j,k) = true;
                    end
                end
                if j~=1
                    if solid(i,j-1,k)
                        fluid_IB(i,j,k) = true;
                    end
                end
                if j~=jtot
                    if solid(i,j+1,k)
                        fluid_IB(i,j,k) = true;
                    end
                end

                if k~=1
                    if solid(i,j,k-1)
                        fluid_IB(i,j,k) = true;
                    end
                end

                if k~=ktot
                    if solid(i,j,k+1)
                        fluid_IB(i,j,k) = true;
                    end
                end

                if (include_diagonals)
                    if (i~=1 && j~=1)
                        if solid(i-1,j-1,k)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && j~=jtot)
                        if solid(i-1,j+1,k)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && j~=1)
                        if solid(i+1,j-1,k)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && j~=jtot)
                        if solid(i+1,j+1,k)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && k~=1)
                        if solid(i-1,j,k-1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && k~=ktot)
                        if solid(i-1,j,k+1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && k~=1)
                        if solid(i+1,j,k-1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && k~=ktot)
                        if solid(i+1,j,k+1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (j~=1 && k~=1)
                        if solid(i,j-1,k-1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (j~=1 && k~=ktot)
                        if solid(i,j-1,k+1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (j~=jtot && k~=1)
                        if solid(i,j+1,k-1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (j~=jtot && k~=ktot)
                        if solid(i,j+1,k+1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && j~=1 && k~=1)
                        if solid(i-1,j-1,k-1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && j~=1 && k~=1)
                        if solid(i+1,j-1,k-1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && j~=jtot && k~=1)
                        if solid(i-1,j+1,k-1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && j~=jtot && k~=1)
                        if solid(i+1,j+1,k-1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && j~=1 && k~=ktot)
                        if solid(i-1,j-1,k+1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && j~=1 && k~=ktot)
                        if solid(i+1,j-1,k+1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && j~=jtot && k~=ktot)
                        if solid(i-1,j+1,k+1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && j~=jtot && k~=ktot)
                        if solid(i+1,j+1,k+1)
                            fluid_IB(i,j,k) = true;
                        end
                    end
                end

%                 if fluid_IB(i,j,k)
%                     fluid_IB_xyz = [fluid_IB_xyz; [xgrid(i), ygrid(j), zgrid(k)]];
%                     bound(i,j,k) = true;
%                     bound_xyz = [bound_xyz; [xgrid(i), ygrid(j), zgrid(k), 1]];
%                 end


            % Identify solid IB points (ghost cells)
            % This could be put into fluid sorting above but easier to log
            % if done like this, and no less efficient despite more code.
            elseif solid(i,j,k)
                nsolid_count = nsolid_count + 1;
                solid_xyz(nsolid_count,:) = [xgrid(i), ygrid(j), zgrid(k)];
                if i~=1
                    if fluid(i-1,j,k)
                        solid_IB(i,j,k) = true;
                    end
                end
                if i~=itot
                    if fluid(i+1,j,k)
                        solid_IB(i,j,k) = true;
                    end
                end
                if j~=1
                    if fluid(i,j-1,k)
                        solid_IB(i,j,k) = true;
                    end
                end
                if j~=jtot
                    if fluid(i,j+1,k)
                        solid_IB(i,j,k) = true;
                    end
                end
                if k~=1
                    if fluid(i,j,k-1)
                        solid_IB(i,j,k) = true;
                    end
                end
                if k~=ktot
                    if fluid(i,j,k+1)
                        solid_IB(i,j,k) = true;
                    end
                end

                if (include_diagonals)
                    if (i~=1 && j~=1)
                        if fluid(i-1,j-1,k)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && j~=jtot)
                        if fluid(i-1,j+1,k)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && j~=1)
                        if fluid(i+1,j-1,k)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && j~=jtot)
                        if fluid(i+1,j+1,k)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && k~=1)
                        if fluid(i-1,j,k-1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && k~=ktot)
                        if fluid(i-1,j,k+1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && k~=1)
                        if fluid(i+1,j,k-1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && k~=ktot)
                        if fluid(i+1,j,k+1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (j~=1 && k~=1)
                        if fluid(i,j-1,k-1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (j~=1 && k~=ktot)
                        if fluid(i,j-1,k+1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (j~=jtot && k~=1)
                        if fluid(i,j+1,k-1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (j~=jtot && k~=ktot)
                        if fluid(i,j+1,k+1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && j~=1 && k~=1)
                        if fluid(i-1,j-1,k-1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && j~=1 && k~=1)
                        if fluid(i+1,j-1,k-1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && j~=jtot && k~=1)
                        if fluid(i-1,j+1,k-1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && j~=jtot && k~=1)
                        if fluid(i+1,j+1,k-1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && j~=1 && k~=ktot)
                        if fluid(i-1,j-1,k+1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && j~=1 && k~=ktot)
                        if fluid(i+1,j-1,k+1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=1 && j~=jtot && k~=ktot)
                        if fluid(i-1,j+1,k+1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                    if (i~=itot && j~=jtot && k~=ktot)
                        if fluid(i+1,j+1,k+1)
                            solid_IB(i,j,k) = true;
                        end
                    end
                end

%                 if solid_IB(i,j,k)
%                     solid_IB_xyz = [solid_IB_xyz; [xgrid(i), ygrid(j), zgrid(k)]];
%                     bound(i,j,k) = true;
%                     bound_xyz = [bound_xyz; [xgrid(i), ygrid(j), zgrid(k), 0]];
%                 end
            end
        end
    end
end

end