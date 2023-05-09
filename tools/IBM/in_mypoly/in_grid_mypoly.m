function solid = in_grid_mypoly(vertices,facets,incenters,faceNormals,xgrid,ygrid,zgrid,L_char,max_height)

solid = false(length(xgrid),length(ygrid),length(zgrid));
for ix = 1:length(xgrid)
    for iy = 1:length(ygrid)
        for iz = 1:length(zgrid)
            if (zgrid(iz)>max_height)
                solid(ix,iy,iz) = false;
            else
                solid(ix,iy,iz) = in_mypoly(vertices,facets,incenters,faceNormals,xgrid(ix),ygrid(iy),zgrid(iz),L_char);
            end

        end
    end
end