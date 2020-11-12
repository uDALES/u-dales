%clear all
%close all

%% read blocks
nheader = 2;
try %in case file is empty -> no blocks
blk = dlmread([tempdir '/blocks.inp.' num2str(expnr)],'',nheader,0);
catch
blk=[];
end

    
nblks=size(blk, 1);
if nblks>0
blk(:, [5,6]) = blk(:,[5,6])+1; % indices for k start at zero?
end
%some dummy grid properties since these are currently not loaded


xc=dlmread([outputdir '/xgrid.inp.' num2str(expnr)],'',2,0);
nx=length(xc);
zc=dlmread([outputdir '/zgrid.inp.' num2str(expnr)],'',2,0);
nz=length(zc);

xb=xh;
zb=zh;
yb=yh;
dx=ones(ni,1)*dx;
dy=ones(nj,1)*dy;
dz=dzf;
nx=ni;
ny=nj;
nz=nk;




%% create Mask-matrix
% this new mask is in x,y coordinates, not y,x coordinates as before
M=zeros(nx,ny);
IM=zeros(nx,ny);
for i=1:size(blk,1)
    xl=blk(i,1);
    xu=blk(i,2);
    yl=blk(i,3);
    yu=blk(i,4);
    M(xl:xu,yl:yu)=1;
    IM(xl:xu,yl:yu)=i;
end

% figure
% imagesc(M)
% set(gca,'YDir','normal')

%% define all facets

% fctl format: orientation, walltype, blockid, buildingid, isinternal
%              il, iu, jl, ju, kl, ku

top = 1; west = 2; east = 3; north = 4; south=5; bot = 6;

il = 1; iu = 2;
jl = 3; ju = 4;
kl = 5; ku = 6;  

nfcts=nblks*6;
fctl=int32(zeros(nfcts,10));
for j=1:nblks
    for k=top:bot
        i=(j-1)*6+k;
        fctl(i,1)=k; % orientation
        fctl(i,2)=blk(j, 6 + min(k,5)); % wall id
        fctl(i,3)=j; % blockid
        
        switch(k)
            case top
                fctl(i, 5 + il) = blk(j, il);
                fctl(i, 5 + jl) = blk(j, jl);
                fctl(i, 5 + kl) = blk(j, ku)+1;
                fctl(i, 5 + iu) = blk(j, iu)+1;
                fctl(i, 5 + ju) = blk(j, ju)+1;
                fctl(i, 5 + ku) = blk(j, ku)+1;
            case west
                fctl(i, 5 + il) = blk(j, il);
                fctl(i, 5 + jl) = blk(j, jl);
                fctl(i, 5 + kl) = blk(j, kl);
                fctl(i, 5 + iu) = blk(j, il);
                fctl(i, 5 + ju) = blk(j, ju)+1;
                fctl(i, 5 + ku) = blk(j, ku)+1;
            case east
                fctl(i, 5 + il) = blk(j, iu)+1;
                fctl(i, 5 + jl) = blk(j, jl);
                fctl(i, 5 + kl) = blk(j, kl);
                fctl(i, 5 + iu) = blk(j, iu)+1;
                fctl(i, 5 + ju) = blk(j, ju)+1;
                fctl(i, 5 + ku) = blk(j, ku)+1;
            case north
                fctl(i, 5 + il) = blk(j, il);
                fctl(i, 5 + jl) = blk(j, ju)+1;
                fctl(i, 5 + kl) = blk(j, kl);
                fctl(i, 5 + iu) = blk(j, iu)+1;
                fctl(i, 5 + ju) = blk(j, ju)+1;
                fctl(i, 5 + ku) = blk(j, ku)+1;
            case south
                fctl(i, 5 + il) = blk(j, il);
                fctl(i, 5 + jl) = blk(j, jl);
                fctl(i, 5 + kl) = blk(j, kl);
                fctl(i, 5 + iu) = blk(j, iu)+1;
                fctl(i, 5 + ju) = blk(j, jl);
                fctl(i, 5 + ku) = blk(j, ku)+1;
            case bot
                fctl(i, 5 + il) = blk(j, il);
                fctl(i, 5 + jl) = blk(j, jl);
                fctl(i, 5 + kl) = blk(j, kl);
                fctl(i, 5 + iu) = blk(j, iu)+1;
                fctl(i, 5 + ju) = blk(j, ju)+1;
                fctl(i, 5 + ku) = blk(j, kl);
        end
    end
end


%% alternative after merging internal blocks
disp('determine internal facets')
intern=zeros(nfcts,2);
j=1;
for i=1:nfcts
    switch(fctl(i,1))
        %only one test necessary, since the whole edge is internal, or not
        %remember that fctl stores facet coordinates, not block coordinates
        case 2
            if fctl(i,5+il)-1>=1 %not at the domain edge
            if M(fctl(i,5+il)-1,fctl(i,5+jl))
                fctl(i,5)=1;
                intern(j,1)=fctl(i,3);
                intern(j,2)=IM(fctl(i,5+il)-1,fctl(i,5+jl));
                j=j+1;
            end
            end
        case 3
            if fctl(i,5+iu)<=nx %not at the domain edge
            if M(fctl(i,5+iu),fctl(i,5+jl))
                fctl(i,5)=1;
                intern(j,1)=fctl(i,3);
                intern(j,2)=IM(fctl(i,5+iu),fctl(i,5+jl));
                j=j+1;
            end
            end
        case 4
            if fctl(i,5+ju)<=ny %not at the domain edge
            if M(fctl(i,5+il),fctl(i,5+ju))
                fctl(i,5)=1;
                intern(j,1)=fctl(i,3);
                intern(j,2)=IM(fctl(i,5+il),fctl(i,5+ju));
                j=j+1;
            end
            end
        case 5
            if fctl(i,5+jl)-1>=1 %not at the domain edge
            if M(fctl(i,5+il),fctl(i,5+jl)-1)
                fctl(i,5)=1;
                intern(j,1)=fctl(i,3);
                intern(j,2)=IM(fctl(i,5+il),fctl(i,5+jl)-1);
                j=j+1;
            end
            end
    end
end
nintern=sum(fctl(:,5));
intern = intern(1:nintern,:);

%%
tic
disp('remove downward and internal facets')
bblk0 = fctl(:, 3);
bblk1 = bblk0;

for n=1:nintern
    % find all facets belonging to block associated with the facet in
    % the second column of intern, and replace these with the block id
    % from the facet in the first column (contagion)
    list=find(bblk0==intern(n,2));
    list2=find(bblk1==bblk1(list(1)));
    bblk1(list2)=intern(n,1);
end
toc

% assign building id to fctl
bblku = unique(bblk1);
nbld = length(bblku);
for n = 1:nbld
    fctl(bblk1 == bblku(n), 4) = n;
end
toc

%% remove downward facets and internal facets
% fctl format: orientation, walltype, blockid, buildingid, isinternal
%              il, iu, jl, ju, kl, ku
sel = find(fctl(:,1)~=bot & fctl(:,5)~=1);
fctl = fctl(sel, :);
nfcts=size(fctl,1);

%% remove facets at domain edge (not actually done)

sel = find(~((fctl(:,6)==1 & fctl(:,7)==1) | (fctl(:,6)==ni+1 & fctl(:,7)==ni+1) | (fctl(:,8)==1 & fctl(:,9)==1) | (fctl(:,8)==nj+1 & fctl(:,9)==nj+1)));

if lradiation && (length(sel)~=nfcts)
myicon = imread('flamingos.jpg');
h=msgbox("Flamingos!!",'Flamingos','custom',myicon);
    
myicon = imread('llama.jpg');
h=msgbox("llama!!",'llama','custom',myicon);

myicon = imread('sherlock.jpg');
h=msgbox("sherlock is a good boy!!",'dog','custom',myicon);
    
cdata = get(0,'DefaultImageCData');
cdata2=uint8(zeros(size(cdata,1),size(cdata,2),3));
blub=(cdata - floor(cdata));
cdata2(:,:,1)= uint8(blub/max(blub(:))*255);
cdata2(:,:,2)= uint8(blub/max(blub(:))*255);
cdata2(:,:,3)= uint8(blub/max(blub(:))*255);
myicon=cdata2;
h=msgbox("Rest of the code will likely crash",'Fluff','custom',myicon);    
    
myicon = imread('peppers.png');
h=msgbox("Don't have blocks on the edge, when using radiation and/or the energy balance!!",'Ratatouille','custom',myicon);



pause()
end

%this 2 lines would actually remove the facets on domain edge
% fctl = fctl(sel, :);
% nfcts=size(fctl);




%% plot orientation, blocks, buildings and wall type
if lhqplot
    figure
    cmap = colormap('parula');
    for i=1:nfcts
        switch fctl(i, 1)
            case {top, bot}
                x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+il))];
                y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+ju)), yb(fctl(i, 5+ju))];
                z = [zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl))];
            case {west, east}
                x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+il)), xb(fctl(i, 5+il)), xb(fctl(i, 5+il))];
                y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+ju)), yb(fctl(i, 5+ju))];
                z = [zb(fctl(i, 5+kl)), zb(fctl(i, 5+ku)), zb(fctl(i, 5+ku)), zb(fctl(i, 5+kl))];
            case {north, south}
                x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+il))];
                y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl))];
                z = [zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl)), zb(fctl(i, 5+ku)), zb(fctl(i, 5+ku))];
        end
        
        subplot(1,3,1)
        ci = min(floor(double(fctl(i, 1))/6*length(cmap))+1, length(cmap));
        patch(x,y,z, cmap(ci, :),'FaceLighting','none');
        hold on
        title('orientation')
        
        subplot(1,3,2)
        ci = min(floor(double(fctl(i, 3))/nblks*length(cmap))+1, length(cmap));
        patch(x,y,z, cmap(ci, :),'FaceLighting','none');
        hold on
        title('blocks')
        
        subplot(1,3,3)
        ci = min(floor(double(fctl(i, 4))/nbld*length(cmap))+1, length(cmap));
        patch(x,y,z, cmap(ci, :),'FaceLighting','none');
        hold on
        title('buildings')
    end
    
    for n =1:3
        subplot(1,3,n)
        view(3)
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis equal
        xlim([0 xb(end)])
        ylim([0 yb(end)])
        zlim([0 zb(end)])
    end
    
    
    %% plot building map
    figure
    subplot(1,2,1)
    for i = find(fctl(:,1) == top)'
        x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+il))];
        y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+ju)), yb(fctl(i, 5+ju))];
        
        ci = min(floor(double(fctl(i, 4))/nbld*length(cmap))+1, length(cmap));
        patch(x,y, cmap(ci, :),'FaceLighting','none');
        hold on
        text(mean(x), mean(y), num2str(fctl(i, 4), '%8d'), ...
            'horizontalalignment', 'center')
    end
    xlabel('x')
    ylabel('y')
    axis equal
    title('building id')
    
    subplot(1,2,2);
    for i=1:nfcts
        switch fctl(i, 1)
            case {top, bot}
                x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+il))];
                y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+ju)), yb(fctl(i, 5+ju))];
                z = [zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl))];
            case {west, east}
                x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+il)), xb(fctl(i, 5+il)), xb(fctl(i, 5+il))];
                y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+ju)), yb(fctl(i, 5+ju))];
                z = [zb(fctl(i, 5+kl)), zb(fctl(i, 5+ku)), zb(fctl(i, 5+ku)), zb(fctl(i, 5+kl))];
            case {north, south}
                x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+il))];
                y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl))];
                z = [zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl)), zb(fctl(i, 5+ku)), zb(fctl(i, 5+ku))];
        end
        
        ci = min(floor(double(fctl(i, 2))/double(max(fctl(:, 2)))*length(cmap))+1, length(cmap));
        patch(x,y,z, cmap(ci, :),'FaceLighting','none');
        d = [0, 0, 0]; a = 0.25;
        switch(fctl(i, 1))
            case top
                d(3) = a * dz(1);
            case west
                d(1) = -a*dx(1);
            case east
                d(1) = a*dx(1);
            case south
                d(2) = -a*dy(1);
            case north
                d(2) = a*dy(1);
        end
        text(mean(x)+d(1), mean(y)+d(2), mean(z)+d(3), num2str(fctl(i, 2)), ...
            'horizontalalignment', 'center')
        hold on
        title('wall type')
    end
    view(3)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    xlim([0 xb(end)])
    ylim([0 yb(end)])
    zlim([0 zb(end)])
    % colorbar
    % return
    
    
    %% plot facets
    
    figure
    for i=1:nfcts
        switch fctl(i, 1)
            case {top, bot}
                x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+il))];
                y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+ju)), yb(fctl(i, 5+ju))];
                z = [zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl))];
            case {west, east}
                x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+il)), xb(fctl(i, 5+il)), xb(fctl(i, 5+il))];
                y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+ju)), yb(fctl(i, 5+ju))];
                z = [zb(fctl(i, 5+kl)), zb(fctl(i, 5+ku)), zb(fctl(i, 5+ku)), zb(fctl(i, 5+kl))];
            case {north, south}
                x = [xb(fctl(i, 5+il)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+iu)), xb(fctl(i, 5+il))];
                y = [yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl)), yb(fctl(i, 5+jl))];
                z = [zb(fctl(i, 5+kl)), zb(fctl(i, 5+kl)), zb(fctl(i, 5+ku)), zb(fctl(i, 5+ku))];
        end
        
        ci = min(floor(double(fctl(i, 2))/double(max(fctl(:, 2)))*length(cmap))+1, length(cmap));
        patch(x,y,z, cmap(ci, :),'FaceLighting','none');
        d = [0, 0, 0]; a = 0.25;
        switch(fctl(i, 1))
            case top
                d(3) = a * dz(1);
            case west
                d(1) = -a*dx(1);
            case east
                d(1) = a*dx(1);
            case south
                d(2) = -a*dy(1);
            case north
                d(2) = a*dy(1);
        end
        text(mean(x)+d(1), mean(y)+d(2), mean(z)+d(3), num2str(i), ...
            'horizontalalignment', 'center')
        hold on
        title('facet nr')
    end
    view(3)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal
    xlim([0 xb(end)])
    ylim([0 yb(end)])
    zlim([0 zb(end)])
    %%
end
%% write list of all facets


% type = type(sel);
nfcts = size(fctl,1);
nblockfcts=nfcts;

if lwritefiles
    
    fname = [outputdir '/facets.inp.' num2str(expnr)];
    
    fileID = fopen(fname,'w');
    fprintf(fileID,'# %4s %6s %6s %6s\n','or', 'wl', 'blk', 'bld');
    fprintf(fileID,'%6d %6d %6d %6d\n', fctl(:, 1:4)');
    fclose(fileID);
end