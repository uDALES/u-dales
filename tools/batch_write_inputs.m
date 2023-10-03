%% Batch write inputs

%expnrs2 = [145, 147, 149, 151, 153, 155, 157, 159, 163, 165, 167, 169, 171, 173, 175, 177, 179, 183, 185, 187, 189, 191, 193, 195, 197, 199];
%expnrs2 = [211:2:219,231:2:239,251:2:259];

expnrs2 = [999];
expnrs1 = expnrs2-1;
ncpus =1;
for e = 1:length(expnrs1)
    expnr = num2str(expnrs1(e));
    expnr2 = num2str(expnrs2(e));
    %     fpath = [DA_EXPDIR '/' expnr '/'];
    %     fpath2 = [DA_EXPDIR '/' expnr2 '/'];
    %     cd(fpath)
    %     ttest = stlread(['geom.' expnr '.stl'])
    %     cd(fpath2)
    %     t2test = stlread(['geom.' expnr2 '.stl'])
    xtiles = 3;
    ytiles = 2;
    tic
    % DA_EXPDIR = getenv('DA_EXPDIR');
    % DA_TOOLSDIR = getenv('DA_TOOLSDIR');
    DA_EXPDIR = '/media/chris/Project3/uDALES2.0/experiments'
    DA_TOOLSDIR = '/media/chris/Project3/uDALES2.0/u-dales/tools'
    addpath(genpath([DA_TOOLSDIR '/']));
    addpath([DA_TOOLSDIR '/IBM/'])
    addpath([DA_TOOLSDIR '/SEB/'])
    exppath = [DA_EXPDIR '/'];
    fpath = [DA_EXPDIR '/' expnr '/'];
    cd(fpath)
    fpath2 = [DA_EXPDIR '/' expnr2 '/'];
    addpath(fpath2);
    
    r = preprocessing(expnr, exppath); % reads namoptions file and creates the object r
    preprocessing.set_defaults(r);
    preprocessing.generate_xygrid(r);
    preprocessing.generate_zgrid(r);
    preprocessing.generate_lscale(r)
    preprocessing.write_lscale(r)
    preprocessing.generate_prof(r);
    preprocessing.write_prof(r);
    
    
    % if r.nsv>0
    %     preprocessing.generate_scalar(r);
    %     preprocessing.write_scalar(r);
    %     disp(['Written scalar.inp.', r.expnr])
    % end
    
    if isfile(['factypes.inp.', expnr])
        r.factypes = dlmread(['factypes.inp.', r.expnr],'',3,0);
    else
        preprocessing.write_factypes(r)
        disp(['Written factypes.inp.', r.expnr])
    end
    
    cd(fpath2)
    r2 = preprocessing(expnr2, exppath);
    preprocessing.set_defaults(r2);
    preprocessing.generate_xygrid(r2);
    preprocessing.generate_zgrid(r2);
    preprocessing.generate_lscale(r2)
    preprocessing.write_lscale(r2)
    preprocessing.generate_prof(r2);
    preprocessing.write_prof(r2);
    if isfile(['factypes.inp.', expnr2])
        r2.factypes = dlmread(['factypes.inp.', r2.expnr],'',3,0);
    else
        preprocessing.write_factypes(r2)
        disp(['Written factypes.inp.', r2.expnr])
    end
    cd(fpath)
    
    if r.libm
        %% Read the .stl file and write necessary ibm files
        TR = stlread(r.stl_file);
        TR2 = stlread(r2.stl_file);
        F = TR.ConnectivityList;
        F2 =TR2.ConnectivityList;
        V = TR.Points;
        V2 = TR2.Points;
        area_facets = facetAreas(F, V); % Useful for checking if area_fluid_IB_c == sum(area_facets)
        area_facets2 = facetAreas(F2, V2);
        %% Set facet types
        nfcts = size(TR.ConnectivityList,1);
        nfcts2 = size(TR2.ConnectivityList,1);
        preprocessing.set_nfcts(r, nfcts);
        preprocessing.set_nfcts(r2, nfcts2);
        facet_types = ones(nfcts,1); % facet_types are to be user-defined - defaults to type 1 (concrete)
        facet_types2 = ones(nfcts2,1);
        %%
        lamstile = lamcal(expnr,area_facets);
        lamsfull = lamcal(expnr2,area_facets2);
    
        %%
        facet_types = type_setter(expnr,TR,r);
        facet_types2 = type_setter(expnr2,TR2,r2);
        %%
        preprocessing.write_facets(r, facet_types, TR.faceNormal);
        preprocessing.write_facets(r2, facet_types2,TR2.faceNormal);
    
    %     lamdba_calculation
    %     setting_types
    
        disp(['Written facets.inp.', r.expnr])
    
        calculate_facet_sections_uvw = r.iwallmom > 1;
        calculate_facet_sections_c = r.ltempeq || r.lmoist;
        if r.gen_geom
            % c-grid (scalars/pressure)
            xgrid_c = r.xf;
            ygrid_c = r.yf;
            zgrid_c = r.zf;
            [X_c,Y_c,Z_c] = ndgrid(xgrid_c,ygrid_c,zgrid_c);
    
            % u-grid
            xgrid_u = r.xh(1:end-1);
            ygrid_u = r.yf;
            zgrid_u = r.zf;
            [X_u,Y_u,Z_u] = ndgrid(xgrid_u,ygrid_u,zgrid_u);
    
            % v-grid
            xgrid_v = r.xf;
            ygrid_v = r.yh(1:end-1);
            zgrid_v = r.zf;
            [X_v,Y_v,Z_v] = ndgrid(xgrid_v,ygrid_v,zgrid_v);
    
            % w-grid
            xgrid_w = r.xf;
            ygrid_w = r.yf;
            zgrid_w = r.zh(1:end-1);
            [X_w,Y_w,Z_w] = ndgrid(xgrid_w,ygrid_w,zgrid_w);
    
            diag_neighbs = r.diag_neighbs;
            stl_ground = r.stl_ground;
            periodic_x = r.BCxm == 1;
            periodic_y = r.BCym == 1;
            xsize = r.xlen;
            ysize = r.ylen;
            zsize = r.zsize;
            
            lmypolyfortran = 1; lmypoly = 0;		% remove eventually
            lwindows = false;
            Dir_ray_u = [0 0 1];
            Dir_ray_v = [0 0 1];
            Dir_ray_w = [0 0 1];
            Dir_ray_c = [0 0 1];
            tol_mypoly = 5e-4;
            toc
            writeIBMFiles_Dipanjan; % Could turn into a function and move writing to this script
            toc
        else
            if isempty(r.geom_path)
                error('Need to specify the path to geometry files')
            end
            copy_command = ['cp ' r.geom_path 'solid_* ' r.geom_path 'fluid_boundary_* ' fpath];
            system(copy_command);
            copy_command = ['cp ' r.geom_path 'fluid_boundary_* ' fpath];
            system(copy_command);
            if calculate_facet_sections_uvw
                copy_command = ['cp ' r.geom_path 'facet_sections_u* ' fpath];
                system(copy_command);
                copy_command = ['cp ' r.geom_path 'facet_sections_v* ' fpath];
                system(copy_command);
                copy_command = ['cp ' r.geom_path 'facet_sections_w* ' fpath];
                system(copy_command);
            end
            if calculate_facet_sections_c
                copy_command = ['cp ' r.geom_path 'facet_sections_c* ' fpath];
                system(copy_command);
            end
        end
        %% Maniupulate data for the larger geom
    
        %% solid_U
        nfcts2 = nfcts*xtiles*ytiles;
        solid_ijk_u_original = solid_ijk_u;
        solid_ijk_u_final = [];
        for i = 0:xtiles-1
            for j = 0:ytiles-1
                solid_ijk_u_new = solid_ijk_u_original;
                solid_ijk_u_new(:,1) = solid_ijk_u_original(:,1)+i*(length(xgrid_u));
                solid_ijk_u_new(:,2) = solid_ijk_u_original(:,2)+j*(length(ygrid_u));
                solid_ijk_u_final = [solid_ijk_u_final; solid_ijk_u_new];
            end 
        end
        solid_ijk_u = solid_ijk_u_final;
        filename_u = [fpath2 'solid_u.txt'];
        fileID_u = fopen(filename_u,'W');
        fprintf(fileID_u, '# position (i,j,k)\n');
        fclose(fileID_u);
        dlmwrite(filename_u, solid_ijk_u, '-append','delimiter',' ');
        disp('Written solid_u.txt')
    
        %% Fluid boundary u
        fluid_IB_ijk_u_original = fluid_IB_ijk_u;
        fluid_IB_ijk_u_final = [];
        ubps_tot = [];
        for i = 0:xtiles-1
            for j = 0:ytiles-1   
                fluid_IB_ijk_u_new = fluid_IB_ijk_u_original;
                fluid_IB_ijk_u_new(:,1) = fluid_IB_ijk_u_original(:,1)+i*(length(xgrid_u));
                fluid_IB_ijk_u_new(:,2) = fluid_IB_ijk_u_original(:,2)+j*(length(ygrid_u));
                if i+j==0
                    ubps_added = 0;
                else 
                    ubps_added = length(fluid_IB_ijk_u_new(:,1));
                end    
                ubps_tot = [ubps_tot, ubps_added];
                fluid_IB_ijk_u_final = [fluid_IB_ijk_u_final; fluid_IB_ijk_u_new];
            end 
        end
        
        ubps_tot = cumsum(ubps_tot);
        fluid_boundary_u = fluid_IB_ijk_u_final;
        fluid_IB_ijk_u = fluid_boundary_u;
        filename_u = [fpath2 'fluid_boundary_u.txt'];
        fileID_u = fopen(filename_u,'W');
        fprintf(fileID_u, '# position (i,j,k)\n');
        fclose(fileID_u);
        dlmwrite(filename_u, fluid_boundary_u, '-append','delimiter',' ');
        disp('Written fluid_boundary_u.txt')
    
        %% Facet sections u
        facet_sections_u_original = facet_sections_u;
        facet_sections_u_final = [];
        count = 1;
        nfcts_u_original = max(facet_sections_u(:,1));
        nfbp_u_original = max(facet_sections_u(:,5));
        while count <= xtiles*ytiles
            facet_sections_u_new = facet_sections_u_original;
            facet_sections_u_new(:,1) = facet_sections_u_original(:,1)+(count-1)*nfcts_u_original;
            facet_sections_u_new(:,5) = facet_sections_u_original(:,5)+ubps_tot(count);
            facet_sections_u_final = [facet_sections_u_final; facet_sections_u_new];
            count = count + 1; 
        end 
        facet_sections_u = facet_sections_u_final;
        filename_u = [fpath2 'facet_sections_u.txt'];
        fileID_u = fopen(filename_u,'W');
        fprintf(fileID_u, '# facet, area, fluid boundary point, distance\n');
        fprintf(fileID_u, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_u(:,[1,2,5,6])');
        fclose(fileID_u);
        disp('Written facet_sections_u.txt')
        
        %% solid v
        solid_ijk_v_original = solid_ijk_v;
        solid_ijk_v_final = [];
        for i = 0:xtiles-1
            for j = 0:ytiles-1
                solid_ijk_v_new = solid_ijk_v_original;
                solid_ijk_v_new(:,1) = solid_ijk_v_original(:,1)+i*(length(xgrid_v));
                %solid_ijk_v_new(:,2) = solid_ijk_v_original(:,2)+j*(length(ygrid_v)-1);
                solid_ijk_v_new(:,2) = solid_ijk_v_original(:,2)+j*(length(ygrid_v));
                solid_ijk_v_final = [solid_ijk_v_final; solid_ijk_v_new];
            end 
        end
        solid_ijk_v = solid_ijk_v_final;
        filename_v = [fpath2 'solid_v.txt'];
        fileID_v = fopen(filename_v,'W');
        fprintf(fileID_v, '# position (i,j,k)\n');
        fclose(fileID_v);
        dlmwrite(filename_v, solid_ijk_v, '-append','delimiter',' ');
        disp('Written solid_v.txt')
    
        %% Fluid IB v
        fluid_IB_ijk_v_original = fluid_IB_ijk_v;
        fluid_IB_ijk_v_final = [];
        vbps_tot = [];
        for i = 0:xtiles-1
            for j = 0:ytiles-1   
                fluid_IB_ijk_v_new = fluid_IB_ijk_v_original;
                fluid_IB_ijk_v_new(:,1) = fluid_IB_ijk_v_original(:,1)+i*(length(xgrid_v));
                fluid_IB_ijk_v_new(:,2) = fluid_IB_ijk_v_original(:,2)+j*(length(ygrid_v));
                if i+j==0
                    vbps_added = 0;
                else 
                    vbps_added = length(fluid_IB_ijk_v_new(:,2));
                end    
                vbps_tot = [vbps_tot, vbps_added];
                fluid_IB_ijk_v_final = [fluid_IB_ijk_v_final; fluid_IB_ijk_v_new];
            end 
        end
        vbps_tot = cumsum(vbps_tot);
        fluid_boundary_v = fluid_IB_ijk_v_final;
        fluid_IB_ijk_v = fluid_boundary_v;
        filename_v = [fpath2 'fluid_boundary_v.txt'];
        fileID_v = fopen(filename_v,'W');
        fprintf(fileID_v, '# position (i,j,k)\n');
        fclose(fileID_v);
        dlmwrite(filename_v, fluid_boundary_v, '-append','delimiter',' ');
        disp('Written fluid_boundary_v.txt')
    
        %% Facet sections v
        facet_sections_v_original = facet_sections_v;
        facet_sections_v_final = [];
        nfcts_v_original = max(facet_sections_v(:,1));
        nfbp_v_original = max(facet_sections_v(:,5));
        count = 1;
        while count <= xtiles*ytiles
            facet_sections_v_new = facet_sections_v_original;
            facet_sections_v_new(:,1) = facet_sections_v_original(:,1)+(count-1)*nfcts_v_original;
            facet_sections_v_new(:,5) = facet_sections_v_original(:,5)+vbps_tot(count);
            facet_sections_v_final = [facet_sections_v_final; facet_sections_v_new];
            count = count + 1; 
        end 
        facet_sections_v = facet_sections_v_final;
        filename_v = [fpath2 'facet_sections_v.txt'];
        fileID_v = fopen(filename_v,'W');
        fprintf(fileID_v, '# facet, area, fluid boundary point, distance\n');
        fprintf(fileID_v, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_v(:,[1,2,5,6])');
        fclose(fileID_v);
        disp('Written facet_sections_v.txt')
    
        %% solid w
        solid_ijk_w_original = solid_ijk_w;
        solid_ijk_w_final = [];
        for i = 0:xtiles-1
            for j = 0:ytiles-1
                solid_ijk_w_new = solid_ijk_w_original;
                solid_ijk_w_new(:,1) = solid_ijk_w_original(:,1)+i*(length(xgrid_w));
                solid_ijk_w_new(:,2) = solid_ijk_w_original(:,2)+j*(length(ygrid_w));
                solid_ijk_w_final = [solid_ijk_w_final; solid_ijk_w_new];
            end 
        end
        solid_ijk_w = solid_ijk_w_final;
        filename_w = [fpath2 'solid_w.txt'];
        fileID_w = fopen(filename_w,'W');
        fprintf(fileID_w, '# position (i,j,k)\n');
        fclose(fileID_w);
        dlmwrite(filename_w, solid_ijk_w, '-append','delimiter',' ');
        disp('Written solid_w.txt')
    
        %% Fluid IB w
        fluid_IB_ijk_w_original = fluid_IB_ijk_w;
        fluid_IB_ijk_w_final = [];
        for i = 0:xtiles-1
            for j = 0:ytiles-1
                fluid_IB_ijk_w_new = fluid_IB_ijk_w_original;
                fluid_IB_ijk_w_new(:,1) = fluid_IB_ijk_w_original(:,1)+i*(length(xgrid_w));
                fluid_IB_ijk_w_new(:,2) = fluid_IB_ijk_w_original(:,2)+j*(length(ygrid_w));
                fluid_IB_ijk_w_final = [fluid_IB_ijk_w_final; fluid_IB_ijk_w_new];
            end 
        end
        fluid_boundary_w = fluid_IB_ijk_w_final;
        fluid_IB_ijk_w = fluid_boundary_w;
        filename_w = [fpath2 'fluid_boundary_w.txt'];
        fileID_w = fopen(filename_w,'W');
        fprintf(fileID_w, '# position (i,j,k)\n');
        fclose(fileID_w);
        dlmwrite(filename_w, fluid_boundary_w, '-append','delimiter',' ');
        disp('Written fluid_boundary_w.txt')
    
        %% Facet sections w
        facet_sections_w_original = facet_sections_w;
        facet_sections_w_final = [];
        count = 0;
        nfcts_w_original = max(facet_sections_w(:,1));
        nfbp_w_original = max(facet_sections_w(:,5));
        for i = 0:xtiles-1
            for j = 0:ytiles-1
                facet_sections_w_new = facet_sections_w_original;
                facet_sections_w_new(:,1) = facet_sections_w_original(:,1)+count*nfcts_w_original;
                facet_sections_w_new(:,5) = facet_sections_w_original(:,5)+count*nfbp_w_original;
                count = count + 1;
                facet_sections_w_final = [facet_sections_w_final; facet_sections_w_new];
            end 
        end
        
        facet_sections_w = facet_sections_w_final;
        filename_w = [fpath2 'facet_sections_w.txt'];
        fileID_w = fopen(filename_w,'W');
        fprintf(fileID_w, '# facet, area, fluid boundary point, distance\n');
        
        fprintf(fileID_w, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_w(:,[1,2,5,6])');
        fclose(fileID_w);
        disp('Written facet_sections_w.txt')
    
        %% Solid points c
        solid_ijk_c_original = solid_ijk_c;
        solid_ijk_c_final = [];
        for i = 0:xtiles-1
            for j = 0:ytiles-1
                solid_ijk_c_new = solid_ijk_c_original;
                solid_ijk_c_new(:,1) = solid_ijk_c_original(:,1)+i*(length(xgrid_c));
                solid_ijk_c_new(:,2) = solid_ijk_c_original(:,2)+j*(length(ygrid_c));
                solid_ijk_c_final = [solid_ijk_c_final; solid_ijk_c_new];
            end 
        end
        solid_ijk_c = solid_ijk_c_final;
        
        filename_c = [fpath2 'solid_c.txt'];
        fileID_c = fopen(filename_c,'W');
        fprintf(fileID_c, '# position (i,j,k)\n');
        fclose(fileID_c);
        dlmwrite(filename_c, solid_ijk_c, '-append','delimiter',' ');
        disp('Written solid_c.txt')
    
        %% Fluid Boundary c
        fluid_IB_ijk_c_original = fluid_IB_ijk_c;
        fluid_IB_ijk_c_final = [];
        for i = 0:xtiles-1
            for j = 0:ytiles-1
                fluid_IB_ijk_c_new = fluid_IB_ijk_c_original;
                fluid_IB_ijk_c_new(:,1) = fluid_IB_ijk_c_original(:,1)+i*(length(xgrid_c));
                fluid_IB_ijk_c_new(:,2) = fluid_IB_ijk_c_original(:,2)+j*(length(ygrid_c));
                fluid_IB_ijk_c_final = [fluid_IB_ijk_c_final; fluid_IB_ijk_c_new];
            end 
        end
        fluid_boundary_c = fluid_IB_ijk_c_final;
        fluid_IB_ijk_c = fluid_boundary_c;
        filename_c = [fpath2 'fluid_boundary_c.txt'];
        fileID_c = fopen(filename_c,'W');
        fprintf(fileID_c, '# position (i,j,k), distance to surface, reconstruction point location\n');
        fclose(fileID_c);
        dlmwrite(filename_c, fluid_boundary_c, '-append','delimiter',' ');
        disp('Written fluid_boundary_c.txt')
    
        %% Facet sections c
        facet_sections_c_original = facet_sections_c;
        facet_sections_c_final = [];
        count = 0;
        nfcts_c_original = max(facet_sections_c(:,1));
        nfbp_c_original = max(facet_sections_c(:,5));
        for i = 0:xtiles-1
            for j = 0:ytiles-1
                facet_sections_c_new = facet_sections_c_original;
                facet_sections_c_new(:,1) = facet_sections_c_original(:,1)+count*nfcts_c_original;
                facet_sections_c_new(:,5) = facet_sections_c_original(:,5)+count*nfbp_c_original;
                count = count + 1;
                facet_sections_c_final = [facet_sections_c_final; facet_sections_c_new];
            end 
        end
        
        facet_sections_c = facet_sections_c_final;
        filename_c = [fpath2 'facet_sections_c.txt'];
        fileID_c = fopen(filename_c,'W');
        fprintf(fileID_c, '# facet, area, fluid boundary point, distance\n');
        fprintf(fileID_c, '%-2d %-4.4f %-4d %-4.8f\n', facet_sections_c(:,[1,2,5,6])');
        fclose(fileID_c);
        disp('Written facet_sections_c.txt')
    
        %% Writing info file
        filename_info = [fpath2 'info.txt'];
        fileID_info = fopen(filename_info,'W');
        fprintf(fileID_info, ['nfcts = ', num2str(nfcts2), '\n']);
        fprintf(fileID_info, ['nsolpts_u = ', num2str(size(solid_ijk_u_final,1)), '\n']);
        fprintf(fileID_info, ['nsolpts_v = ', num2str(size(solid_ijk_v_final,1)), '\n']);
        fprintf(fileID_info, ['nsolpts_w = ', num2str(size(solid_ijk_w_final,1)), '\n']);
        fprintf(fileID_info, ['nsolpts_c = ', num2str(size(solid_ijk_c_final,1)), '\n']);
        %fprintf(fileID_info, ['nbndpts_u = ', num2str(size(fluid_IB_xyz_u,1)), '\n']);
        fprintf(fileID_info, ['nbndpts_u = ', num2str(size(fluid_IB_ijk_u_final,1)), '\n']);
        % fprintf(fileID_info, ['nbndpts_v = ', num2str(size(fluid_IB_xyz_v,1)), '\n']);
        % fprintf(fileID_info, ['nbndpts_w = ', num2str(size(fluid_IB_xyz_w,1)), '\n']);
        % fprintf(fileID_info, ['nbndpts_c = ', num2str(size(fluid_IB_xyz_c,1)), '\n']);
        fprintf(fileID_info, ['nbndpts_v = ', num2str(size(fluid_IB_ijk_v_final,1)), '\n']);
        fprintf(fileID_info, ['nbndpts_w = ', num2str(size(fluid_IB_ijk_w_final,1)), '\n']);
        fprintf(fileID_info, ['nbndpts_c = ', num2str(size(fluid_IB_ijk_c_final,1)), '\n']);
        fprintf(fileID_info, ['nfctsecs_u = ', num2str(size(facet_sections_u_final,1)), '\n']);
        fprintf(fileID_info, ['nfctsecs_v = ', num2str(size(facet_sections_v_final,1)), '\n']);
        fprintf(fileID_info, ['nfctsecs_w = ', num2str(size(facet_sections_w_final,1)), '\n']);
        fprintf(fileID_info, ['nfctsecs_c = ', num2str(size(facet_sections_c_final,1)), '\n']);
        fprintf(fileID_info, ['sinkbase = ', num2str(max(solid_ijk_c(:,3))), '\n']);
        fclose(fileID_info);
    
        %%
        if r.lEB
            toc
            preprocessing.write_facetarea(r, area_facets);
    
    
            %% Write STL in View3D input format
            fpath_facets_view3d = [fpath2 'facets.vs3'];
            STLtoView3D(r2.stl_file, fpath_facets_view3d);
    
            %% Calculate view factors
            % remember to build View3D in local system windows/linux
            % Add check to see if View3D exists in the tools directory.
            if lwindows
                view3d_exe = [DA_TOOLSDIR '/View3D/src/View3D.exe'];
            else
                view3d_exe = [DA_TOOLSDIR '/View3D/build/src/view3d'];
            end
            fpath_vf = [fpath2 'vf.txt'];
            vf = view3d(view3d_exe, fpath_facets_view3d, fpath_vf);
    %         vftile = view3d(view3d_exe, fpath_facets_view3d, fpath_vf);
    %         vf = diagonal_matrix_maker(vftile,xtiles,ytiles);
    %         svftile = max(1 - sum(vftile, 2), 0);
            svf = max(1 - sum(vf, 2), 0);
    %         preprocessing.write_svf(r2, svf);
    %         preprocessing.write_facetarea(r2, area_facets2);
    %         if ~r2.lvfsparse
    %             preprocessing.write_vf(r2, vf)
    %             disp(['Written vf.nc.inp.', r2.expnr])
    %         else
    %             vfsparse = sparse(double(vf));
    %             preprocessing.write_vfsparse(r2, vfsparse);
    %             disp(['Written vfsparse.inp.', r2.expnr])
    %         end
    
    
    
            %% Calculate direct solar radiation (Sdir)
            disp('Calculating direct solar radiation.')
            azimuth = r2.solarazimuth - r2.xazimuth;
            nsun = [sind(r2.solarzenith)*cosd(azimuth), -sind(r2.solarzenith)*sind(azimuth), cosd(r2.solarzenith)];
            show_plot_2d = false; % User-defined
            show_plot_3d = true;  % User-defined
            Sdir = directShortwave(F2, V2, nsun, r2.I, r2.psc_res, show_plot_2d, show_plot_3d);
    
            %% Calculate net shortwave radiation (Knet)
            disp('Calculating net shortwave radiation.')
            albedos = preprocessing.generate_albedos(r2, facet_types2);
            %Knet_tile = netShortwave(Sdir, r.Dsky, vftile, svftile, albedos);
            Knet = netShortwave(Sdir, r2.Dsky, vf, svf, albedos);
            preprocessing.write_netsw(r2, Knet);
            %tiles = xtiles*ytiles;
    %         Knet = zeros(length(Knet_tile)*tiles,1);
    %         for i = 0:tiles-1
    %             xl = 1+i*length(Knet_tile);
    %             xu = length(Knet_tile) + i*length(Knet_tile);
    %             Knet(xl:xu) = Knet_tile;
    %         end 
    %         preprocessing.write_netsw(r2, Knet);
    %         disp(['Written netsw.inp.', r2.expnr])
        end
    
        %% Write initial facet temperatures
        if (r.lEB || r.iwallmom == 2 || r.iwalltemp == 2)
            disp('Setting initial facet temperatures.')
            facT = r.facT;
            facT2 = r2.facT;
            nfaclyrs = r.nfaclyrs;
            nfaclyrs2 = r2.nfaclyrs;
            facT_file = r.facT_file;
            facT_file2 = r2.facT_file;
            lfacTlyrs = r.lfacTlyrs;
            if ~r.lfacTlyrs
                Tfacinit = ones(nfcts,1) .* r.facT;
                preprocessing.write_Tfacinit(r, Tfacinit)
                disp(['Written Tfacinit.inp.', r.expnr])
                % Could always read in facet temperature as layers, defaulting to linear?
            else
                Tfac = ncread(r.facT_file, 'T');
                Tfacinit_layers = Tfac(:, :, end);
                preprocessing.write_Tfacinit_layers(r, Tfacinit_layers)
                disp(['Written Tfacinit_layers.inp.', r.expnr])
            end
            cd(fpath2)
            if ~r2.lfacTlyrs
                Tfacinit2 = ones(nfcts2,1) .* r2.facT;
                preprocessing.write_Tfacinit(r2, Tfacinit2)
                disp(['Written Tfacinit.inp.', r2.expnr])
                % Could always read in facet temperature as layers, defaulting to linear?
            else
                Tfac2 = ncread(r2.facT_file, 'T');
                Tfacinit_layers2 = Tfac2(:, :, end);
                preprocessing.write_Tfacinit_layers(r2, Tfacinit_layers2)
                disp(['Written Tfacinit_layers.inp.', r2.expnr])
            end
    
        end
    
        %% Writing last inputs for expnr2 
        cd(fpath2)
        preprocessing.write_svf(r2, svf);
        preprocessing.write_facetarea(r2, area_facets2);
        if ~r2.lvfsparse
            preprocessing.write_vf(r2, vf)
            disp(['Written vf.nc.inp.', r2.expnr])
        else
            vfsparse = sparse(double(vf));
            preprocessing.write_vfsparse(r2, vfsparse);
            disp(['Written vfsparse.inp.', r2.expnr])
        end
        preprocessing.write_netsw(r2, Knet);
        disp(['Written netsw.inp.', r2.expnr])
    end
    %% setting effective albedo
    %cd(fpath2)
    %efctvalb = 1-sum(Knet)/sum((Sdir+r.Dsky*svf));
    efctvalb = 1-sum(area_facets2.*Knet)/(sum(area_facets2.*(Sdir))+sum(r.Dsky*area_facets2.*svf))
    %preprocessing.write_efalb(r,efctvalb);
    dlmwrite([fpath2 'Sdir.' expnr2], Sdir);
    preprocessing.write_efalb(r2,efctvalb);
    % cd(fpath2)
    % efctvalb2 = 1-sum(Knet)/sum(Sdir+r.Dsky*svf);
    % preprocessing.write_efalb(r2,efctvalb2);
    %%
    fac_type_table = r2.factypes;
    ems = [];
    for i = 1:r2.nfcts
        typ = facet_types2(i);
        typind = find(fac_type_table(:,1)==typ);
        em = fac_type_table(typind,6);
        ems = [ems,em];
    end 
    dlmwrite(['emissivity.' expnr2], ems');
    
    toc
end 

function newmatrix = diagonal_matrix_maker(m, nxtiles, nytiles)
    % Determine the size of the matrix m
    [rowm, colm] = size(m);
    % Specify the desired size of matrix B (m*m)
    tiles = nxtiles*nytiles;
    totx = rowm*tiles;
    toty = colm*tiles;
    B = zeros(totx,toty);
    for i = 0:tiles-1
            il = 1+i*rowm;
            iu = rowm+i*rowm;
            jl = 1+i*colm;
            ju = colm+i*colm;
            B(il:iu,jl:ju) = m;
    end
    % Display the resulting matrix B
    newmatrix = B;
end 

function lams = lamcal(expnr, area_facets)

    stem = '/media/chris/Project3/uDALES2.0/experiments/';
    fpathstl = [stem expnr '/geom_new.' expnr '.stl'];
    %fpathstl = [stem expnr '/geom.' expnr '.stl'];
    fpathg = [stem  expnr '/green_outline.' expnr];
    fpathfacarea = [stem expnr '/facetarea.inp.' expnr];
    fpathgnfcts = [stem expnr '/green_facets.' expnr];
    TR = stlread(fpathstl);
    %area_facets = dlmread(fpathfacarea,'', 1,0);
    con = TR.ConnectivityList;
    ps = TR.Points;
    norms = TR.faceNormal;
    Area = max(ps(:,1))*max(ps(:,2));
    west_area = 0;
    east_area = 0;
    green_area = 0;
    top_area = 0;
    green_outline = dlmread(fpathg);
    green_facets = [];
    for i = 1:length(con(:,1))
        verts =con(i,:);
        p1 = ps(verts(1),:);
        p2 = ps(verts(2),:);
        p3 = ps(verts(3),:);
        zsum = p1(3) + p2(3) + p3(3);
        if zsum == 0
            for g = 1:length(green_outline(:,1))
                gsp = green_outline(g,:);
                gxl = gsp(1);
                gxu = gsp(2);
                gyl = gsp(3);
                gyu = gsp(4);
                if all([and(p1(1)>=gxl,p1(1)<=gxu),...
                    and(p1(2)>=gyl,p1(2)<=gyu),...
                    and(p2(1)>=gxl,p2(1)<=gxu),...
                    and(p2(2)>=gyl,p2(2)<=gyu),...
                    and(p3(1)>=gxl,p3(1)<=gxu),...
                    and(p3(2)>=gyl,p3(2)<=gyu)])
                    area = area_facets(i);
                    green_area = green_area + area;
                    green_facets = [green_facets; i];
                end
            end
        else 
            fnorm = norms(i,:);
            if fnorm == [-1,0,0]
                area = area_facets(i);
                west_area = west_area + area;
    
            elseif fnorm == [1,0,0]
                area = area_facets(i);
                east_area = east_area + area;
            elseif fnorm == [0,0,1]
                area = area_facets(i);
                top_area = top_area + area;
            end     
    
        end
    end
    lp = top_area/Area;
    lv = green_area/Area;
    lf = west_area/Area;
    dlmwrite(fpathgnfcts, green_facets);
    disp('Done calculating lambdas.');
    lams = [lp,lv,lf];
end

function types = type_setter(expnr,TR,r)
    %% Setting green facets
    fpath = ['/media/chris/Project3/uDALES2.0/experiments/' expnr '/']; 
    green_facets = dlmread([fpath '/green_facets.' expnr]);
    nfcts = r.nfcts;
    facet_types = ones(nfcts,1);
    facet_types(green_facets) = 12;
    
    filename_facets = [fpath 'facets.inp.' expnr];
    fileID_facets = fopen(filename_facets,'W');
    fprintf(fileID_facets, '# type, normal\n');
    fclose(fileID_facets);
    dlmwrite(filename_facets, [facet_types, TR.faceNormal], '-append','delimiter',' ','precision','%1.0f');
    
    disp('Done setting types.');
    types = facet_types;
end