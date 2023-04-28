function inside = in_mypoly(vertices,facets,incenters,faceNormals,xgrid,ygrid,zgrid,L_char)

nfcets = length(facets(:,1));
inside = false(length(xgrid),length(ygrid),length(zgrid));

%Take x axis, or y axis or positive z axis
% even with carefull cohice of this ray direction there might be wrong
% prediction. To resolve this one must shoot at the least three rays along
% the three cartresian direction and decide based on the majority outcome.
Dir_ray = [0 0 1]';       %% not in the user choice. should be decided carefully

for ix = 1:length(xgrid)
    for iy = 1:length(ygrid)
        for iz = 1:length(zgrid)
            counter = 0;
            Origin = [xgrid(ix) ygrid(iy) zgrid(iz)]';   
            for i_facet = 1:nfcets
                if (incenters(i_facet,1) > Origin(1)-L_char && incenters(i_facet,1) < Origin(1)+L_char && incenters(i_facet,2) > Origin(2)-L_char && incenters(i_facet,2) < Origin(2)+L_char)  %% not in the user choice. should be decided carefully
                    x_tri = [vertices(facets(i_facet,1),1) vertices(facets(i_facet,2),1) vertices(facets(i_facet,3),1)]';
                    y_tri = [vertices(facets(i_facet,1),2) vertices(facets(i_facet,2),2) vertices(facets(i_facet,3),2)]';
                    z_tri = [vertices(facets(i_facet,1),3) vertices(facets(i_facet,2),3) vertices(facets(i_facet,3),3)]';
                    flag = point_triangle_intersect(Origin,Dir_ray,x_tri,y_tri,z_tri);
                    if (flag==true)
                        
                        
                        % Special case need to be included here origin points lying on the surface of spherical objects
                        isOnFacet = point_on_triangle(Origin,x_tri,y_tri,z_tri,incenters(i_facet,:)',faceNormals(i_facet,:)');
                        if (isOnFacet==true)
                            inside(ix,iy,iz) = true; %%inside
                            break;
                        else
                        
                            if (counter == 0)
                                counter=counter+1;
                                facet_matrix_old(:,:,counter) = [x_tri y_tri z_tri];
                            else
                                facet_matrix = [x_tri y_tri z_tri];
                                count2 = 0;
                                for kct = 1:3
                                    for ic=1:counter
                                        for jc = 1:3
                                            if (facet_matrix(kct,:) == facet_matrix_old(jc,:,ic))
                                                count2 = 1;
                                            end
                                        end
                                    end
                                end
                                if (count2 == 0)
                                    counter=counter+1;
                                    facet_matrix_old(:,:,counter) = [x_tri y_tri z_tri];
                                end
                            end
                        end
                    end
                end
            end
            
            if(inside(ix,iy,iz)~=true)
                if (mod(counter,2) == 0)
                    inside(ix,iy,iz) = false; %% outside
                else
                    inside(ix,iy,iz) = true; %%inside
                end
            end
        end
    end
end