%clear all
%close all

%% add a bounding wall around the domain used in radiation calculations only
%% THIS DOES NOT YET DEAL WITH PERIODIC GEOMETRY (I.E. CANYONS)



%% read files
%blocks
nheader=2;
try %in case file is empty -> no blocks
    
    %can use blocks instead of sliced blocks (bbri)  [I guess]
    %B = dlmread([tempdir '/blocks.inp.' num2str(expnr)],'',nheader,0);  %#il   iu   jl    ju   kl   ku  ttop  twest teast tnor tsou
    %B = dlmread([tempdir '/blocks.inp.' num2str(expnr)],'',nheader,0);  %#il   iu   jl    ju   kl   ku  ttop  twest teast tnor tsou
    B = obj.blocks;
catch
    B =[];
end


%facets
nheader=1;
try
    %F = dlmread([outputdir '/facets.inp.' num2str(expnr)],'',nheader,0);  %#   or     wl    blk    bld
    F = obj.facets;
catch
    F=[];
end

[nblocks, nbi]=size(B);
[nfcts, nfacprop]=size(F);

height=floor(median(B(:,6)));

if nblocks>0
    if any(B(:,1)==1)
        disp('Building(s) at lower x domain edge')
    end
    if any(B(:,2)==nx)
        disp('Building(s) at upper x domain edge')
    end
    if any(B(:,3)==1)
        disp('Building(s) at lower y domain edge')
    end
    if any(B(:,4)==ny)
        disp('Building(s) at upper y domain edge')
    end
end

nxwalls=ceil(ny/maxsize);
remx=rem(ny,maxsize);
nywalls=ceil(nx/maxsize);
remy=rem(nx,maxsize);
nzw=ceil((height+1)/maxsize);
remz=rem((height+1),maxsize);

boundingwalls=zeros(2*nzw*(nxwalls+nywalls),7);
boundingwallfacets=zeros(2*nzw*(nxwalls+nywalls),4);

%%

%in the west east north south
%facing east west south north
if remx>0
    for j=1:nzw
        if ((j==nzw) && (remz>0))
            lh=height-remz+1;
            uh=height;
        else
            lh=(j-1)*maxsize;
            uh=j*maxsize-1;
        end
        for i=1:(nxwalls-1)
            boundingwalls((i-1)*nzw+j,:)=[1 1 (i-1)*maxsize+1 i*maxsize lh uh -101];
            boundingwalls((i-1)*nzw+j+nxwalls*nzw,:)=[nx nx (i-1)*maxsize+1 i*maxsize lh uh -101];
        end
        boundingwalls((nxwalls-1)*nzw+j,:)=[1 1 ny-remx+1 ny lh uh -101];
        boundingwalls((nxwalls-1)*nzw+j+nxwalls*nzw,:)=[nx nx ny-remx+1 ny lh uh -101];
    end
else
    for j=1:nzw
        if ((j==nzw) && (remz>0))
            lh=height-remz+1;
            uh=height;
        else
            lh=(j-1)*nzw;
            uh=j*nzw-1;
        end
        for i=1:nxwalls
            boundingwalls((i-1)*nzw+j,:)=[1 1 (i-1)*maxsize+1 i*maxsize lh uh -101];
            boundingwalls((i-1)*nzw+j+nxwalls*nzw,:)=[nx nx (i-1)*maxsize+1 i*maxsize lh uh -101];
        end
    end
end

if remy>0
    for j=1:nzw
        if ((j==nzw) && (remz>0))
            lh=height-remz+1;
            uh=height;
        else
            lh=(j-1)*maxsize;
            uh=j*maxsize-1;
        end
        for i=1:(nywalls-1)
            boundingwalls(2*nzw*nxwalls+(i-1)*nzw+j,:)=[(i-1)*maxsize+1 i*maxsize ny ny lh uh -101];
            boundingwalls(2*nzw*nxwalls+(i-1)*nzw+j+nywalls*nzw,:)=[(i-1)*maxsize+1 i*maxsize 1 1 lh uh -101];
        end
        boundingwalls(2*nzw*nxwalls+(nywalls-1)*nzw+j,:)=[nx-remy+1 nx ny ny lh uh -101];
        boundingwalls(2*nzw*nxwalls+(nywalls-1)*nzw+j+nywalls*nzw,:)=[nx-remy+1 nx 1 1 lh uh -101];
    end
else
    for j=1:nzw
        if ((j==nzw) && (remz>0))
            lh=height-remz+1;
            uh=height;
        else
            lh=(j-1)*nzw;
            uh=j*nzw-1;
        end
        for i=1:nywalls
            boundingwalls(2*nzw*nxwalls+(i-1)*nzw+j,:)=[(i-1)*maxsize+1 i*maxsize ny ny lh uh -101];
            boundingwalls(2*nzw*nxwalls+(i-1)*nzw+j+nywalls*nzw,:)=[(i-1)*maxsize+1 i*maxsize 1 1 lh uh -101];
        end
    end
end


for i=1:(nxwalls*nzw)
    boundingwallfacets(i,:)=[3 -101 i -101]; %west, facing east
    boundingwallfacets((nxwalls*nzw)+i,:)=[2 -101 i+nxwalls*nzw -101]; %east facing west
end
for i=1:(nywalls*nzw)
    boundingwallfacets(2*(nxwalls*nzw)+i,:)=[5 -101 2*nxwalls*nzw+i -101]; %north, facing south
    boundingwallfacets(2*(nxwalls*nzw)+(nywalls*nzw)+i,:)=[4 -101 2*nxwalls*nzw+i+nywalls*nzw -101]; %south, facing north
end

%if lwritefiles
%    fname = [tempdir '/boundingwalls.txt']; %it's not an input to the les
%    fileID = fopen(fname,'w');
%    fprintf(fileID,'# %4s\n','bounding wall facets');
%    fprintf(fileID,'# %4s\n','indeces & walltype');
%    fprintf(fileID,'# %4s %6s %6s %6s %6s %6s %6s\n','il', 'iu', 'jl', 'ju', 'kl', 'ku', 't');
%    fprintf(fileID,'%6d %6d %6d %6d %6d %6d %3d\n', boundingwalls(:, 1:7)');
%    fclose(fileID);

obj.facets(end+1:end+size(boundingwallfacets,1),1:4) = boundingwallfacets;

%    dlmwrite([outputdir '/facets.inp.' num2str(expnr)],boundingwallfacets,'delimiter',' ','precision','%7d','-append')
%end

%% append fctl
% fctl format: orientation, walltype, blockid, buildingid, isinternal
%              il, iu, jl, ju, kl, ku
%if lplot
top = 1; west = 2; east = 3; north = 4; south=5; bot = 6;

il = 1; iu = 2;
jl = 3; ju = 4;
kl = 5; ku = 6;
nwallfcts=size(boundingwalls,1);
if nfcts==0 %create fctl
    fctl=zeros(nwallfcts,11);
    for j=1:nwallfcts
        fctl(j, :) = [2 , 1, -101, -101, 0, boundingwalls(j, iu)+1, boundingwalls(j, iu)+1, boundingwalls(j, jl), boundingwalls(j, ju)+1, boundingwalls(j, kl)+1, boundingwalls(j, ku)+1+1];
    end
else %append
    for j=1:nwallfcts
        switch(boundingwallfacets(j,1))
            case west
                fctl(end+1, :) = [2 , 1, -101, -101, 0, boundingwalls(j, iu)+1, boundingwalls(j, iu)+1, boundingwalls(j, jl), boundingwalls(j, ju)+1, boundingwalls(j, kl)+1, boundingwalls(j, ku)+1+1];
            case east
                fctl(end+1, :) = [3 , 1, -101, -101, 0, boundingwalls(j, il), boundingwalls(j, il), boundingwalls(j, jl), boundingwalls(j, ju)+1, boundingwalls(j, kl)+1, boundingwalls(j, ku)+1+1];
            case north
                fctl(end+1, :) = [4 , 1, -101, -101, 0, boundingwalls(j, il), boundingwalls(j, iu)+1, boundingwalls(j, jl),  boundingwalls(j, jl), boundingwalls(j, kl)+1, boundingwalls(j, ku)+1+1];
            case south
                fctl(end+1, :) = [5 , 1, -101, -101, 0, boundingwalls(j, il), boundingwalls(j, iu)+1, boundingwalls(j, ju)+1, boundingwalls(j, ju)+1,  boundingwalls(j, kl)+1, boundingwalls(j, ku)+1+1];
        end
    end
end
%end