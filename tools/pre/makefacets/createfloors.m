%% Create floor
%appends facets.inp.xxx  (alternatively write floorfacets.inp.xxx)
%fill space between blocks with floors (streets etc.)
%create file floors.inp.xxx identical to facets.inp.101
%floors only have x and y coordinates, with building ID = 0

%CAREFUL, DON'T RUN MORE THAN ONCE. IT APPENDS TO FACETS.INP.XXX
%                                      -------



%% read blocks
%read blocks

nheader=2;
try %in case file is empty -> no blocks
    
%can use blocks instead of sliced blocks (bbri)  [I guess]
%B = dlmread([tempdir '/blocks.inp.' num2str(expnr)],'',nheader,0);  %#il   iu   jl    ju   kl   ku  ttop  twest teast tnor tsou
%B = dlmread([tempdir '/blocks.inp.' num2str(expnr)],'',nheader,0);  %#il   iu   jl    ju   kl   ku  ttop  twest teast tnor tsou
B = obj.blocks;
catch
B =[];
end

[nblocks, nbi]=size(B);  %number of blocks, number of block parameters


%% create Mask-matrix
%
M=ones(nx,ny);
BI=zeros(nx,ny); %block index mask
corm=zeros(nx,ny); %mask with all wall-floor corners
cornm=zeros(nx,ny); %mask with all the wall-wall-floor corners, 8 = NW, 10 = SW, 12 = NE, 15=SE

for i=1:nblocks
    xl=B(i,1);
    xu=B(i,2);
    yl=B(i,3);
    yu=B(i,4);
    M(xl:xu,yl:yu)=0;
    BI(xl:xu,yl:yu)=i;
end
NM=1-M;

if ltestplot
    % figure
    % imagesc(xb,yb,M')
    % axis equal tight
    % set(gca,'YDir','normal')
    % title('building mask')
    figure
    imagesc(xb+1,yb-1,BI')
    axis equal tight
    set(gca,'YDir','normal')
    title('building index mask')
end
%% make them around blocks first, with 1 blocksize wide
%numbers indicate floors, "--" and "|" walls  (only 1 blocksize wide, just two
%number used to make diagonal clear)
%
%   ----------------------
% |
% |    -------------------
% |   |3311111111111111111
% |   |3311111111111111111
% |   |22
% |   |22
% |   |22
%after detsub
%   ----------------------
% |
% |    -------------------
% |   |2111111111111111111
% |   |2211111111111111111
% |   |22
% |   |22
% |   |22

maxblocks=sum(M(:)); %there cant be more blocks then number of grid cells
floors=NaN(maxblocks,4); %allocate a maximum size for floors, reduce size later. xy coordinates of corners, counterclockwise, starting nortwest

c=0;
M2=M;
iM=zeros(size(M)); %save indeces
for i=1:nblocks
    xl=B(i,1);
    xu=B(i,2);
    yl=B(i,3);
    yu=B(i,4);
    
    %west
    if xl-1>=1 %not at domain edge
        if BI(xl-1,yl)==0 % left neighbour is a floor
            c=c+1;
            M2(xl-1,yl:yu)=0; %set to 2 for later check if it is a corner
            iM(xl-1,yl:yu)=c;
            
            floors(c,:)=[xl-1 xl-1 yl yu];
            if (yl-1>=1) && (xu+1<=nx) && (yu+1<=ny) && (xl-1>=1)      
                if BI(xl-1,yl-1)>0 %corner with a north wall
                    cornm(xl-1,yl)=8;
                elseif BI(xl-1,yu+1)>0 %corner with a south wall
                    cornm(xl-1,yu)=10;
                end
            end
        end
    end
    %east
    if xu+1<=nx %not at domain edge
        if BI(xu+1,yl)==0
            c=c+1;
            M2(xu+1,yl:yu)=0;
            iM(xu+1,yl:yu)=c;
            floors(c,:)=[xu+1 xu+1 yl yu];
            if (yl-1>=1) && (xu+1<=nx) && (yu+1<=ny) && (xl-1>=1)   
                if BI(xu+1,yl-1)>0 %corner with a north wall
                    cornm(xu+1,yl)=12;
                elseif BI(xu+1,yu+1)>0  %corner with a south wall
                    cornm(xu+1,yu)=15;
                end
            end
        end
    end
    %north
    if yu+1<=ny %not aat domain edge
        if BI(xu,yu+1)==0
            c=c+1;
            M2(xl:xu,yu+1)=0;
            iM(xl:xu,yu+1)=c;
            floors(c,:)=[xl xu yu+1 yu+1];
            if (yl-1>=1) && (xu+1<=nx) && (yu+1<=ny) && (xl-1>=1)     
                if BI(xl-1,yu+1)>0 %corner with an east wall
                    cornm(xl,yu+1)=12;
                elseif BI(xu+1,yu+1)>0 %corner with a west wall
                    cornm(xu,yu+1)=8;
                end
            end
        end
    end
    %south
    if yl-1>=1 %not at domain edge
        if BI(xu,yl-1)==0
            c=c+1;
            M2(xl:xu,yl-1)=0;
            iM(xl:xu,yl-1)=c;
            floors(c,:)=[xl xu yl-1 yl-1];
            if (yl-1>=1) && (xu+1<=nx) && (yu+1<=ny) && (xl-1>=1)
                if BI(xl-1,yl-1)>0 %corner with an east wall
                    cornm(xl,yl-1)=15;
                elseif BI(xu+1,yl-1)>0 %corner with a west wall
                    cornm(xu,yl-1)=10;
                end
            end
        end
    end
end

corm(iM>0)=1;

if ltestplot
    % figure
    % imagesc(xb,yb,M2')
    % axis equal tight
    % set(gca,'YDir','normal')
    % title('floor mask after adding floors around buildings')
    
    figure
    imagesc(xb,yb,iM')
    axis equal tight
    set(gca,'YDir','normal')
    title('floor indeces after adding floors around buildings')
    
    figure
    blub=corm'+NM'.*0.5;
    imagesc(xb,yb,blub)
    axis equal tight
    set(gca,'YDir','normal')
    title('mask of all corners floor-wall')
    
    figure
    imagesc(xb,yb,cornm')
    axis equal tight
    set(gca,'YDir','normal')
    title('mask of all corners floor-wall-wall')
end

%% remove identical facets in corners (if it's a 1x1 facet in both cases)
%truncate matrix
lnan=find(isnan(floors(:,1)));
if ~isempty(lnan)
    floors(lnan(1):lnan(end),:)=[];
end
floors = unique(floors,'rows','stable');

indexarea=(floors(:,2)-floors(:,1)+1).*(floors(:,4)-floors(:,3)+1);

nfloors=size(floors,1);
count=1;
while count<=nfloors
    i=count  ;
    if sum(floors(:,1)<=floors(i,1) & floors(:,2)>=floors(i,2) & floors(:,3)<=floors(i,3) & floors(:,4)>=floors(i,4))>1
        floors(i,:)=[]; %this floor is contained within another and can be removed
        nfloors=nfloors-1;
    else
        count=count+1;
    end
end

c=size(floors,1);
%% Make floors
% make them in 1D first (fixed x, along y)

while any(M2(:)>0)
    for i=1:nx
        ls=find(M2(i,:)==1);
        if ~isempty(ls)
            first=ls(1);
            if length(ls)>1
                last=ls(find(diff(ls)~=1,1));
                if isempty(last)
                    last=min(ny,first+maxsize-1);
                else
                    last=min(last,first+maxsize-1);
                end
            else
                last=first;
            end
            c=c+1;
            floors(c,:)=[i i first last];
            M2(i,first:last)=0;
            iM(i,first:last)=c;
        end
    end
end

% if ltestplot
% figure
% imagesc(xb,yb,M2')
% axis equal tight
% set(gca,'YDir','normal')
% title('floormask after filling')
% figure
% imagesc(xb,yb,iM')
% axis equal tight
% set(gca,'YDir','normal')
% title('floor indeces')
% end

%truncate matrix
lnan=find(isnan(floors(:,1)));
if ~isempty(lnan)
    floors(lnan(1):lnan(end),:)=[];
end

% figure
% subplot(1,2,1)
% imagesc(M')
% set(gca,'YDir','normal')
% axis equal tight
% title('before')
% subplot(1,2,2)
% imagesc(M2')
% set(gca,'YDir','normal')
% axis equal tight
% title('after')

%%
% combine same size in other dimension
nslice=size(floors,1);
%floors2=NaN(size(floors));
floors2=floors;
c=0;

dsize=1;
sizeold=nslice;
while dsize>0
    i=1;
    while 1
        a=floors2(i,2); %this floors xu
        if corm(a,floors2(i,3)) %his is a floor that belongs to a corner with a wall, don't merge with others
            i=i+1;
            continue
        end
        bv=find(floors2(:,1)==(a+1)); %all floors with xl == this floors xu+1
        b2=bv(find( (floors2(bv,3)==floors2(i,3)) & (floors2(bv,4)==floors2(i,4)) )); %all of these floors witch also have the same y dimensions (should only be one)
        if ~isempty(b2) && ~corm(floors2(b2,1),floors2(b2,3))
            floors2(i,2)=floors2(b2,2);
            floors2(b2,:)=[];
            %  floors2(b2,2)=[];floors2(b2,4)=[];floors2(b2,1)=[];floors2(b2,3)=[];
        end
        i=i+1;
        if i>=size(floors2,1)
            break
        end
    end
    dsize=sizeold-size(floors2,1);
    sizeold=size(floors2,1);
end


ls=999;

while ~isempty(ls)
    nfloors=size(floors2,1);
    ls=find(floors2(:,2)-floors2(:,1)>maxsize);
    floors2=[floors2; NaN(length(ls),4)];
    for i=1:length(ls)
        ind=ls(i);
        floors2(nfloors+i,:)=floors2(ind,:);
        floors2(ind,2)=floors2(ind,1)+maxsize-1;
        floors2(nfloors+i,1)=floors2(ind,2)+1;
    end
end

floors3=zeros(size(floors2,1),5);
floors3(:,1:4)=floors2; %indeces
floors3(:,5)=-1; %type


%% merge floors in y, where possible and as long as smaller than maxsize, don't merge triple corners
change=true;
while change
    change=false;
    for j=1:size(floors3,1)
        il=floors3(j,1);
        iu=floors3(j,2);
        jl=floors3(j,3);
        ju=floors3(j,4);
        if sum(sum(cornm(il:iu,jl:ju)))==0 %no triple corner somewhere on this floor facet, try to merge along y
            
            flu=find(floors3(:,1)==il & floors3(:,2)==iu & floors3(:,3)==ju+1); %floor with same x dimension on ju+1
            fll=find(floors3(:,1)==il & floors3(:,2)==iu & floors3(:,4)==jl-1); %floor with same x dimension on jl-1
            
            if ~isempty(flu)
                ilu=floors3(flu,1);
                iuu=floors3(flu,2);
                jlu=floors3(flu,3);
                juu=floors3(flu,4);
                if sum(sum(cornm(ilu:iuu,jlu:juu)))==0 && (floors3(flu,4)-floors3(j,3)+1<maxsize)
                    floors3(j,4)=floors3(flu,4);
                    floors3(flu,:)=[];
                    change=true;
                end  
            elseif ~isempty(fll)
                ill=floors3(fll,1);
                iul=floors3(fll,2);
                jll=floors3(fll,3);
                jul=floors3(fll,4);
                if sum(sum(cornm(ill:iul,jll:jul)))==0 && (floors3(j,4)-floors3(fll,3)+1<maxsize)
                    floors3(j,3)=floors3(fll,3);
                    floors3(fll,:)=[];
                    change=true;
                end
            end
            
        end
        if change
        break    
        end
    end
end
nfloors=size(floors3,1);

if ltestplot
    %rebuild indexmask and plot
    xc=xb+0.5; xc(end)=[];
    yc=yb+0.5; yc(end)=[];
    
    blub=zeros(nx,ny);
    for i=1:size(floors3,1)
        blub(floors3(i,1):floors3(i,2),floors3(i,3):floors3(i,4))=i;
    end
    
    figure
    imagesc(xc,yc,blub')
    %axis equal tight
    title('floor indeces with borders outlined')
    set(gca,'YDir','normal')
    hold on
    for i=1:size(floors3,1)
        rectangle('Position',[xc(floors3(i,1))-0.5 yc(floors3(i,3))-0.5  xc(floors3(i,2))-xc(floors3(i,1))+1 yc(floors3(i,4))-yc(floors3(i,3))+1])
    end
end

%% append fctl
% fctl format: orientation, walltype, blockid, buildingid, isinternal
%              il, jl, kl, iu, ju, ku
%if ltestplot
for j=1:size(floors3,1)
    fctl(end+1, :) = [1 , 1, -99, -99, 0, floors3(j,1), floors3(j,2)+1, floors3(j,3), floors3(j,4)+1, 1, 1];
end
%end


%% write
nfcts=size(fctl,1);

disp([num2str(nfcts) ' facets, of which: ' num2str(nblockfcts) ' from buildings, ' num2str(nwallfcts) ' from walls, ' num2str(nfloors) ' from floors.'])

%if lwritefiles
%    
%    fname = [tempdir '/facetnumbers.txt']; %it's not an input to the les
%    fileID = fopen(fname,'w');
%    fprintf(fileID,'# %4s %4s %4s %4s %4s \n','nblocks','nfloors','nblockfcts', ' nwallfcts', ' nfloorfcts');
%    fprintf(fileID,'%6i %6i %6i %6i %6i\n', [nblocks nfloors nblockfcts nwallfcts nfloors]);
%    fclose(fileID);
    
%    fname = [tempdir '/floors.txt'];
%    %fname = 'floors.txt'; %it's not an input to the les
%    fileID = fopen(fname,'w');
%    fprintf(fileID,'# %4s\n','floor facets');
%    fprintf(fileID,'# %4s\n','indeces & walltype');
%    fprintf(fileID,'# %4s %6s %6s %6s %6s\n','il', 'iu', 'jl', 'ju', 't');
%    fprintf(fileID,'%6d %6d %6d %6d %3d\n', floors3(:, 1:5)');
%    fclose(fileID);
       
    floors2=zeros(size(floors3,1),4);
    floors2(:,1)=1;   %orientation
    floors2(:,2)=-1;  %type
    floors2(:,4)=-1;   %corresponding building (=none)
    
    for i=1:size(floors3,1)
       % floors2(i,3)=i+nblocks;
        floors2(i,3)=i;
    end
    
%    dlmwrite([outputdir '/facets.inp.' num2str(expnr)],floors2,'delimiter',' ','precision','%7d','-append')
 
obj.facets(end+1:end+size(floors2,1),1:4) = floors2;

%end

%% write blocks & floorblocks
%if lwritefiles
    
    if nblocks>0
    Btw=B(:,1:6);
    %increase the z index of all blocks by one
    Btw(:,5:6)=Btw(:,5:6)+1;
    else
    Btw=zeros(nfloors,6);
    end
    

    
    %set z index of floors to 0
    Btw((nblocks+1):(nblocks+nfloors),:)=zeros(nfloors,6);
    Btw((nblocks+1):(nblocks+nfloors),1:4)=floors3(:,1:4);
    Btw(:,7:11)=zeros(nblocks+nfloors,5);
    
%     % add floors below buildings
%     % do we need/want this? DALES will loop over more blocks but not really
%     % do anything, on the other hand the statistics might look better?
%     % If we use this, then the number of blocks in modibm should
%     be different from the number of blocks in masking matrices to avoid
%     looping. Also these blocks don't have corresponding facets and thus
%     access element 0 of any facet array in DALES
%     Btw(end+1:end+nblocks,:)=Btw(1:nblocks,:)
%     Btw(end-nblocks+1:end,5:end)=zeros(nblocks,7)


    
    %#il   iu   jl    ju   kl   ku  fatop  fawest faeast fanor fasou
    %add the corresponding facets
    for i=1:nblockfcts %for blocks
        Btw(fctl(i,3),fctl(i,1)+6)=i;
    end
    j=nblocks+1;
    for i=(nblockfcts+nwallfcts+1):nfcts %for floors
        Btw(j,7)=i;
        j=j+1;
    end
    
   % cd(tempdir)
    %copyfile bbri.inp bbri2.inp
    %fname = [tempdir '/blocks.inp.' num2str(expnr)];
    %fileID = fopen(fname,'W');
    %fprintf(fileID, '# block location\n');
    %fprintf(fileID, '#il    iu    jl     ju    kl    ku   corresp. facets\n');
    %fclose(fileID);
    %dlmwrite(fname,Btw,'-append','delimiter','\t','precision','%4i')
    %[SUCCESS,MESSAGE,MESSAGEID] = copyfile([tempdir '/blocks.inp.' num2str(expnr)], [outputdir '/blocks.inp.' num2str(expnr)]);
%     [SUCCESS,MESSAGE,MESSAGEID] = movefile('bbri.inp', ['blocks.inp.' num2str(expnr)]) ;
%     [SUCCESS,MESSAGE,MESSAGEID] = movefile('bbri2.inp', 'bbri.inp') ;
%     [SUCCESS,MESSAGE,MESSAGEID] = copyfile(['blocks.inp.' num2str(expnr)], [outputdir '/blocks.inp.' num2str(expnr)]) ;
    %cd(parentdir)
%end

obj.blocks(end+1:end+size(Btw,1),1:11) = Btw;

%% plot all facets

%%
scale=2;
scalef=1.5;
if lhqplot
        il = 1; iu = 2;
    jl = 3; ju = 4;
    kl = 5; ku = 6;
    %   F = dlmread(['facets.inp.' num2str(expnr)],'',nheader,0);  %#   or     wl    blk    bld
    %  [nfcts, nfacprop]=size(F);
    h= figure;
    set(gcf,'units','centimeters','position',[0 0 14.5*scale 14.5*scale]);
    set(h,'PaperPosition',[0 0 14.5*scale 14.5*scale]);
    set(h,'PaperUnits','centimeters');
    set(h,'renderer','painters');
    %set(h,'renderer','opengl');
    cmap = colormap('parula');
    title("facets")
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
%         if (i==99)
%             x=[9, 10, 11 , 9];
%             y=[8,  8, 9 , 9];
%         elseif (i==103)
%             x=[10, 11, 11, 10];
%             y=[7, 7, 9, 8];
%         end
        if(fctl(i,3)<0)
            add=5;
        else
            add=0;
        end
        ci = min(floor(double((fctl(i,1)+add))/11*length(cmap)), length(cmap));
        p=patch(x,y,z, cmap(ci, :),'FaceLighting','none');
        d = [0, 0, 0]; a = 0.2; b=-0.2;
        switch(fctl(i, 1))
            case top
                %d(3) = a * dz(1);
                d=[b*dx(1), b*dy(1), a*dz(1)];
            case west
                %d(1) = -a*dx(1);
                d=[-0.27*dx(1) b*dy(1) -b/2*dz(1)];
            case east
                %d(1) = a*dx(1);
                d=[0.42*dx(1) b*dy(1) -b/2*dz(1)];
            case south
                %d(2) = -a*dy(1);
                d=[b/2*dx(1) -a*dy(1) -b/2*dz(1)];
            case north
                d(2) = a*dy(1);
        end
        t=text(mean(x)+d(1), mean(y)+d(2), mean(z)+d(3), num2str(i), 'horizontalalignment', 'center','Interpreter','latex','Fontsize',12);
        hold on
    end
    view(3)
    
    axis equal
    xlim([0 xb(end)])
    ylim([0 yb(end)])
    %zlim([0 72])   % 
    zlim([0 zb(end)])
    
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12*scalef)
    %h1=gca;
    % h1.Position=[0.08 0.1100 0.78 0.8150];
    
    %caxis([0 71])
    xlabel('x [-]','Interpreter','latex','FontSize',12*scalef)
    ylabel('y [-]','Interpreter','latex','FontSize',12*scalef)
    zlabel('z [-]','Interpreter','latex','FontSize',12*scalef)
    
    %colormap(flipud(bone)) %colormap(flipud(gray))
    
    %hcb=colorbar('Position',[0.92 0.15 0.03 0.69]);
    %hcb.Label.Interpreter='latex';
    %hcb.TickLabelInterpreter='latex';
    %title(hcb,'height [m]','Interpreter','latex','FontSize',12)
    
    
    hold off
    
    %print -depsc2 facetindeces.eps
    %print -dpng facetindeces.png
    
    set(gcf, 'Color', 'w');
  %  export_fig facetindeces.eps
    
end
