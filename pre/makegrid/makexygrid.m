if source == 2
% actual resulting resolution, overwriting the desired one
dx=round(dx/dxinp*ni)/(dxinp*ni);
dy=round(dy/dyinp*nj)/(dyinp*nj);
end
   
% create grid
xf=dx/2:dx:(ni*dx-dx/2);
yf=dy/2:dy:(nj*dy-dy/2);
xh=0:dx:(ni*dx+dx/2);
yh=0:dy:(nj*dy+dy/2);
    
if lwritefiles
    fname = [outputdir '/xgrid.inp.' num2str(expnr)]; 
    fileID = fopen(fname,'w');
    fprintf(fileID,'# %4s \n','x-grid');
    fprintf(fileID,'# %4s \n','');
    fprintf(fileID,'%8.4f\n', xf);
    fclose(fileID);
end
