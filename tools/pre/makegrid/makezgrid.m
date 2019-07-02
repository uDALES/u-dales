if any(ismember(['exp','tanh','2tanh'], stretch))
    disp(['z-grid will be stretched: ' stretch])
else
    disp('z-grid will be linear')
end


%if source == 2
    
    
    if dz*nk==dh  %linear
        if any(ismember(['exp','tanh','2tanh'], stretch))
            disp('no stretching needed')
        end
        zf=dz/2:dz:dh-dz/2;
        zh=0:dz:dh;
        
        dzf = zeros(nk,1);

        for k=1:nk
            dzf(k) = zh(k+1)-zh(k);  
        end
        
    else %streched
        if stretch=="tanh"
            tanh_mesh
        elseif stretch=="exp"
            exp_mesh
        elseif stretch=="2tanh"
            doubletanh_mesh
        end
    end
%end


if lwritefiles
    fname = [outputdir '/zgrid.inp.' num2str(expnr)];
    fileID = fopen(fname,'w');
    fprintf(fileID,'# %4s \n','z-grid');
    fprintf(fileID,'# %4s \n','');
    fprintf(fileID,'%8.4f\n', zf);
    fclose(fileID);
end