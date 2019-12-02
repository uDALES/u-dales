%just make y slices of width 1
%add y slices of same dimension
%split slices up according to xy rules
%determine and merge internal blocks along y
%split slicers accroding to z rules

%%testblock
% topomask=[0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
%         0 1 1 0 0 0 0 0 0 0 0 0 0 0;...
%         0 1 1 1 0 0 0 0 0 0 0 0 0 0;...
%         0 0 1 1 0 0 0 0 0 0 0 1 1 0;...
%         0 0 1 1 1 1 0 0 1 1 1 1 1 0;...
%         0 1 1 1 1 1 1 1 1 1 1 1 0 0;...
%         0 1 1 1 1 1 1 1 1 0 1 1 1 0;...
%         0 0 0 0 1 1 1 0 0 0 0 1 1 0;...
%         0 0 0 0 1 1 0 0 0 0 0 0 0 0;...
%         0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
%         0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% max=topomask;
% ni=14; nj=11;


%% just make y slices of width 1
maxnrblocks=sum(topomask(:)); %allocate arrays with maximum size they can possibly have, reduce size later
xmin=zeros(maxnrblocks,1); %store lower x bound of blocks
xmax=zeros(maxnrblocks,1); %store upper x bound of blocks
ymin=zeros(maxnrblocks,1); %store lower y bound of blocks
ymax=zeros(maxnrblocks,1); %store upper y bound of blocks
zmin=zeros(maxnrblocks,1); %store lower z bound of blocks
zmax=zeros(maxnrblocks,1); %store upper z bound of blocks

blchecked=zeros(size(topomask)); %mask of blocks which have already been checked
indexmask=zeros(size(topomask)); %mask of the indeces of the (new) blocks

count=0;
for i=1:ni %loop over all x
    xminl=i;
    xmaxl=i;
    j=1;
    while j<=nj %loop along j
        if topomask(j,i)==0 %not a building
            j=j+1; %check next j
            continue
        else
            
            yminl=j;
            
            heightgradient=diff(topo(j:end,i));
            
            if isempty(heightgradient) % if at the end of the y-direction (je)
            heightchangey=j;
            elseif all(heightgradient==0) % same block/ floor until the end of the domain
            heightchangey=nj;    
            else
            heightchangey=find(heightgradient~=0,1)+j-1;  %last cell with same height as j (i.e. there is a height change betwen this and the next cell
            end
            
%            if topo(heightchangey,1)>0 && topo(heightchangey+1,1)>0 && heightgradient(heightchangey)<0 % same building different block height! % not interested as currently setting LOWER nonzero block change

            for jj=yminl+1:heightchangey
                
                % if there is a smaller block on the left or right and
                % there is also a change in height between that block and
                % it's neighbour.
                
%                 if ( i==1 ) 
%                     
%                     if ( ( any(topo(yminl-1:yminl,i+1)<topo(yminl-1:yminl,i) & topo(yminl-1:yminl,i+1)>0) ) && ( topo(yminl,i+1)~= topo(yminl-1,i+1) ) )
%                 
%                         heightchangey = jj-1;
%                     
%                     end
%                     
%                 elseif ( i==ni )
%                     
%                     if ( any(topo(yminl-1:yminl,i-1)<topo(yminl-1:yminl,i) & topo(yminl-1:yminl,i-1)>0) && ( topo(yminl,i-1)~= topo(yminl-1,i-1) ) )
%                     
%                         heightchangey = jj-1;
%                     
%                     end
%                     
%                 else
%                     if ( any(topo(yminl-1:yminl,i-1)<topo(yminl-1:yminl,i) & topo(yminl-1:yminl,i-1)>0) && ( topo(yminl,i-1)~= topo(yminl-1,i-1) ) ) || ( ( any(topo(yminl-1:yminl,i+1)<topo(yminl-1:yminl,i) & topo(yminl-1:yminl,i+1)>0) ) && ( topo(yminl,i+1)~= topo(yminl-1,i+1) ) )
%                     
%                         heightchangey = jj-1; % overwrite heightchangey as we need to truncate block earlier!
%                     
%                     end
%                     
%                 end

                if ( i==1 ) 
                    
                    if ( ( any(topo(jj-1:jj,i+1)>0) ) && ( topo(jj,i+1)~= topo(jj-1,i+1) ) )
                
                        1
                        
                        heightchangey = jj-1;
                        
                        break
                    
                    end
                    
                elseif ( i==ni )
                    
                    if ( any(topo(jj-1:jj,i-1)>0) && ( topo(jj,i-1)~= topo(jj-1,i-1) ) )
                    
                        2
                        
                        heightchangey = jj-1;
                        
                        break
                    
                    end
                    
                else
                    
                    if ( any(topo(jj-1:jj,i-1)>0) && ( topo(jj,i-1)~= topo(jj-1,i-1) ) ) || ( ( any(topo(jj-1:jj,i+1)>0) ) && ( topo(jj,i+1)~= topo(jj-1,i+1) ) )
                        
                        heightchangey = jj-1; % overwrite heightchangey as we need to truncate block earlier!
                        
                        break
                    
                    end
                    
                end

                
            end
                
            ymaxl=heightchangey;  % end of the current block (either floor or block with different height comes next)
            
            ztemp = zeros(6,1);
            
            % tg3315 changed from commented below as can have yminl==1 and
            % ymaxl==nj
            if yminl==1
                ztemp(2,1) = NaN;
            else
                ztemp(2,1) = topo(yminl-1,xminl);
            end
            
            if xminl==1
                ztemp(3,1) = NaN;
            else
                ztemp(3,1) = topo(yminl,xminl-1);
            end
            
            if xmaxl==ni
                ztemp(4,1) = NaN;
            else
                ztemp(4,1) = topo(yminl,xminl+1);
            end
            
            if ymaxl==nj
                ztemp(5,1) = NaN;
            else
                ztemp(5,1) = topo(yminl+1,xminl);
            end
           
%             if xminl==1 && yminl==1
%                 ztemp(2,1) = NaN; ztemp(3,1) = NaN; ztemp(4,1) = topo(yminl,xmaxl+1); ztemp(5,1) = topo(ymaxl+1,xmaxl);
%             elseif xmaxl==ni && yminl==1
%                 ztemp(2,1) = NaN; ztemp(3,1) = topo(yminl,xminl-1); ztemp(4,1) = NaN; ztemp(5,1) = topo(ymaxl+1,xmaxl);
%             elseif xminl==1 && ymaxl==nj
%                 ztemp(2,1) = topo(yminl-1,xminl); ztemp(3,1) = NaN; ztemp(4,1) = topo(yminl,xmaxl+1); ztemp(5,1) = NaN;
%             elseif xmaxl==ni && ymaxl==nj
%                 ztemp(2,1) = topo(yminl-1,xminl); ztemp(3,1) = topo(yminl,xminl-1); ztemp(4,1) = NaN; ztemp(5,1) = NaN;
%             elseif xminl==1
%                 ztemp(2,1) = topo(yminl-1,xminl); ztemp(3,1) = NaN; ztemp(4,1) = topo(yminl,xmaxl+1); ztemp(5,1) = topo(ymaxl+1,xmaxl);
%             elseif xaxl==ni
%                 ztemp(2,1) = topo(yminl-1,xminl); ztemp(3,1) = topo(yminl,xminl-1); ztemp(4,1) = NaN; ztemp(5,1) = topo(ymaxl+1,xmaxl);
%             elseif yminl==1
%                 ztemp(2,1) = NaN; ztemp(3,1) = topo(yminl,xminl-1); ztemp(4,1) = topo(yminl,xmaxl+1); ztemp(5,1) = topo(ymaxl+1,xmaxl);
%             elseif ymaxl==nj
%                 ztemp(2,1) = topo(yminl-1,xminl); ztemp(3,1) = topo(yminl,xminl-1); ztemp(4,1) = topo(yminl,xmaxl+1); ztemp(5,1) = NaN;
%             else
%                 ztemp(2,1) = topo(yminl-1,xminl); ztemp(3,1) = topo(yminl,xminl-1); ztemp(4,1) = topo(yminl,xmaxl+1); ztemp(5,1) = topo(ymaxl+1,xmaxl);
%             end
            
            zcuttemp = sort(ztemp(ztemp<topo(yminl, xminl) & ztemp>0));
            zcut = [0; zcuttemp; topo(yminl,xminl)];
            
            for kc=1:size(zcut,1)-1
                
                count=count+1;
                
                blchecked(yminl:ymaxl,xminl:xmaxl)=blchecked(yminl:ymaxl,xminl:xmaxl)+1;
                indexmask(yminl:ymaxl,xminl:xmaxl)=count;
                xmin(count)=xminl;
                xmax(count)=xmaxl;
                ymin(count)=yminl;
                ymax(count)=ymaxl;
                if zcut(kc,1)==0
                    zmin(count)=0;
                else
                    zmin(count)=zcut(kc,1)/dz+1;
                end
                zmax(count)=zcut(kc+1,1)/dz;
                
            end

            j=heightchangey+1;  % move to index after current block
            
        end
    end
end
if ltestplot
    figure
    imagesc(xf,yf,indexmask)
    set(gca,'YDir','normal')
    xlim([xh(1) xh(end)])
    ylim([yh(1) yh(end)])
title('after cutting into y slices')
end

%shorten arrays
ymin((count+1):end)=[];
ymax((count+1):end)=[];
xmin((count+1):end)=[];
xmax((count+1):end)=[];
zmin((count+1):end)=[];
zmax((count+1):end)=[];
disp(['Nr of blocks after cutting into y slices: ' num2str(count)])



%% add y slices of same dimension, aggregate along x
%keep old xmax etc.. to keep better track whats happening
xmax2=xmax; ymax2=ymax; zmax2=zmax;
xmin2=xmin; ymin2=ymin; zmin2=zmin;
dsize=1;
sizeold=count;
while dsize>0  %do as long as there are unmerged blocks % tg3315 need to also check here that you do not merge if they have a block above...!
    i=1;   
    while 1  %do along x
        a=xmax2(i);  %upper x index of this block
        bv=find(xmin2==(a+1));  %all blocks with a lower x bound 1 bigger than this blocks upper x bound
        b2=bv(ymin2(bv)==ymin2(i)); %all of those blocks with also the same lower y bound
        b3=b2(zmin2(b2)==zmin2(i));
        if ~isempty(b3);
            if all(ymax2(b3)==ymax2(i)) && all(zmax2(b3)==zmax2(i)) && topo(ymin2(b3),xmin2(b3))==topo(ymin2(i),xmin2(i)) %&& all(zmin2(b2)==zmin2(i)) %if they also have the same upper y bound and the same height
                % additional check to make sure we do not merge blocks in a way that causes internal-external facets at this point
                if (topo(max(1,ymin2(b3)-1),xmax2(b3))==topo(max(1,ymin2(i)-1),xmax2(i))) && (topo(min(ymax2(b3)+1,nj),xmax2(b3))==topo(min(ymax2(i)+1,nj),xmax2(i)))
                    xmax2(i)=xmax2(b3); %merge 
                    xmax2(b3)=[];ymax2(b3)=[];zmax2(b3)=[];xmin2(b3)=[];ymin2(b3)=[];zmin2(b3)=[]; %remove the just merged block from list of blocks
                end
            end
        end
        i=i+1; %check furhter along x
        if i>=length(xmin2)
            break
        end
    end
    dsize=sizeold-length(xmax2);
    sizeold=length(xmax2);
    
end
count2=sizeold;
disp(['Nr of blocks after merging y slices of same size along x: ' num2str(count2)])

%make fields again
datamean1=zeros(size(topomask));
datamean2=zeros(size(topomask));
datamean3=zeros(size(topomask));
indexmask2=zeros(size(topomask));

%make new matrices
for i=1:count
    datamean1(ymin(i):ymax(i),xmin(i):xmax(i))=zmax(i);
end
for i=1:count2
    indexmask2(ymin2(i):ymax2(i),xmin2(i):xmax2(i))=i;
    datamean2(ymin2(i):ymax2(i),xmin2(i):xmax2(i))=zmax2(i);
end

if ltestplot
    figure
    imagesc(xf,yf,indexmask2)
    set(gca,'YDir','normal')
    xlim([xh(1) xh(end)])
    ylim([yh(1) yh(end)])
title('after merging y slices of same size along x')
end

%% tg and bs 30.05.19 - try to split blocks vertically too

% xmax25=xmax2; ymax25=ymax2; zmax25=zmax2;
% xmin25=xmin2; ymin25=ymin2; zmin25=zmin2;
% 
% for i=1:ni %loop over all x


%% At this point can form block file that works for old cold

% Edit for NY in 1930s topology
% dummy=ones(5,count2);
% 
% for k=1:length(zmax2)
%     if zmax2(k)>12
%         zmax2(k) = 12+round(12*rand(1));
%     end
% end
% 
% fileID = fopen([outputdir '/blocksOLD.inp.' num2str(expnr)],'w');
% fprintf(fileID,'# Block data\n');
% fprintf(fileID,'#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t dtop\t dwest\t deast\t dnor\t dsou\n');
% fprintf(fileID,'%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n',[xmin2';xmax2';ymin2';ymax2';zmin2';zmax2';dummy]);
% fclose(fileID);

%% split slices up according to rules in x and y, do z later
xmin3=zeros(maxnrblocks,1);
xmax3=zeros(maxnrblocks,1);
ymin3=zeros(maxnrblocks,1);
ymax3=zeros(maxnrblocks,1);
zmin3=zeros(maxnrblocks,1);
zmax3=zeros(maxnrblocks,1);

xmin3(1:count2)=xmin2;
xmax3(1:count2)=xmax2;
ymin3(1:count2)=ymin2;
ymax3(1:count2)=ymax2;
zmin3(1:count2)=zmin2;
zmax3(1:count2)=zmax2;

if lradiation

indexmask3=indexmask2;
change=true;
count3=count2;

maxx = 0;
cmax = 0;

while change  %do until there is no changes anymore
    
    change=false;
    
    for c=1:count3 %check if all blocks have same dimension as neighbour
        %only have to check x-1 and x+1, since there can't be a
        %neightbour on y-1 or y+1 (due to the inital slicing along y)
        %index of blocks to the left
        
        if ~(xmin3(c)==1) %on left domain boundary, don't need to check for neighbouring blocks
            illb=indexmask3(ymin3(c):ymax3(c),xmin3(c)-1);
            iillb=find(illb>0,1);
            if ~isempty(illb)
                ilb=illb(iillb);
            else
                ilb=0;
            end
        else
            ilb=0;
        end
        
        if ~(xmax3(c)==ni) %on right domain boundary, don't check for neighbouring blocsk
            irrb=indexmask3(ymin3(c):ymax3(c),xmax3(c)+1);
            iirrb=find(irrb>0,1);
            if ~isempty(irrb)
                irb=irrb(iirrb);
            else
                irb=0;
            end
        else
            irb=0;
        end
        %   irt=indexmask(ymin(c),xmax(c)+1);
        
        if (irb==0) %no block to the right
            yminr=ymin3(c); %use y coordinates of this block
            ymaxr=ymax3(c);
        else
            yminr=ymin3(irb);  %use y coordinates of block to the right
            ymaxr=ymax3(irb);
        end
        
        if (ilb==0) %no block to the left
            ymaxl=ymax3(c); %use y coordinates of this block
            yminl=ymin3(c);
        else
            ymaxl=ymax3(ilb); %use y coordinates of block to the left
            yminl=ymin3(ilb);
        end
        
        if xmin(3)>maxx
        maxx = xmin3(c);
        maxx/ni
        end

        if c/count3>cmax
            cmax=c/count3
        end
        
        %
        %
            
        if (ymaxl<ymax3(c))% & zmax3(ilb)==zmax3(c)  %this block extends further (up) than left neighbour, this block has to be split
            %x  x  x  x      %x  x  x  x
            %x    [c] x      %x    [n] x
            %x [l][c] x      %x [l][c] x
            %x [l][c] x ===> %x [l][c] x
            %x [l][c] x      %x [l][c] x
            %x  x  x  x      %x  x  x  x
            change=true;            %keep checking
            count3=count3+1;        %add new block
            ymin3(count3)=ymaxl+1;  %new block starts 1 above xmax of left neighbour
            ymax3(count3)=ymax3(c); %new block ends at xmax of this block
            xmin3(count3)=xmin3(c); %new block has same x coordinates as this block
            xmax3(count3)=xmax3(c);
            zmin3(count3)=zmin3(c); %keep the same heights
            zmax3(count3)=zmax3(c);
            indexmask3(ymin3(count3):ymax3(count3),xmin3(count3):xmax3(count3))=count3;
            ymax3(c)=ymaxl;         %shorten this block
            break                   %continue forloop from the start
        
        elseif (yminl>ymin3(c))% & zmin3(ilb)==zmin3(c)  %this block extends further (down) than left neighbour, this block has to be split
            %x  x  x  x      %x  x  x  x
            %x [l][c] x      %x [l][c] x
            %x [l][c] x      %x [l][c] x
            %x [l][c] x ===> %x [l][c] x
            %x    [c] x      %x    [n] x
            %x  x  x  x      %x  x  x  x
            change=true; %keep checking
            count3=count3+1; %add new block
            ymin3(count3)=ymin3(c); %new block starts at same ymin as this
            ymax3(count3)=yminl-1;  %new block ends 1 below ymin of left neighbour
            xmin3(count3)=xmin3(c); %new block has same x coordinates as this block
            xmax3(count3)=xmax3(c);
            zmin3(count3)=zmin3(c); %keep the same heights
            zmax3(count3)=zmax3(c);
            indexmask3(ymin3(count3):ymax3(count3),xmin3(count3):xmax3(count3))=count3;  %overwrite the index in the indexmask with the new one
            ymin3(c)=yminl;         %shorten this block
            break                   %continue forloop from the start
        
        elseif (yminr>ymin3(c))% & zmin3(irb)==zmin3(c) %this block extends further (down) than reight neighbour, this block has to be split
            %x  x  x  x      %x  x  x  x
            %x [c][r] x      %x [c][r] x
            %x [c][r] x      %x [c][r] x
            %x [c][r] x ===> %x [c][r] x
            %x [c]    x      %x [n]    x
            %x  x  x  x      %x  x  x  x
            change=true;            %keep checking
            count3=count3+1;        %add new block
            ymin3(count3)=ymin3(c); %new block starts at same ymin as this
            ymax3(count3)=yminr-1;  %new block ends 1 below ymin of reight neighbour
            xmin3(count3)=xmin3(c); %new block has same x coordinates as this block
            xmax3(count3)=xmax3(c);
            zmin3(count3)=zmin3(c); %keep the same heights
            zmax3(count3)=zmax3(c);
            indexmask3(ymin3(count3):ymax3(count3),xmin3(count3):xmax3(count3))=count3; %overwrite the index in the indexmask with the new one
            ymin3(c)=yminr;         %shorten this block
            break                   %continue for-loop from the start
            
        elseif (ymaxr<ymax3(c))% & zmax3(irb)==zmax3(c) %this block extends further (up) than reight neighbour, this block has to be split
            %x  x  x  x      %x  x  x  x
            %x [c]    x      %x [n]    x
            %x [c][r] x      %x [c][r] x
            %x [c][r] x ===> %x [c][r] x
            %x [c][r] x      %x [c][r] x
            %x  x  x  x      %x  x  x  x
            change=true;             %keep checking
            count3=count3+1;         %add new block
            ymin3(count3)=ymaxr+1 ;  %new block starts 1 above xmax of reight neighbour
            ymax3(count3)=ymax3(c);  %new block ends at xmax of this block
            xmin3(count3)=xmin3(c);  %new block has same x coordinates as this block
            xmax3(count3)=xmax3(c);
            zmin3(count3)=zmin3(c); %keep the same heights
            zmax3(count3)=zmax3(c);
            indexmask3(ymin3(count3):ymax3(count3),xmin3(count3):xmax3(count3))=count3;
            ymax3(c)=ymaxr; %shorten this block
            break %continue forloop from the start
            
        end
        
    end
    
    
end
disp(['Nr of blocks after applying rules for radiation in x and y : ' num2str(count3)])
ymin3((count3+1):end)=[];
ymax3((count3+1):end)=[];
xmin3((count3+1):end)=[];
xmax3((count3+1):end)=[];
zmin3((count3+1):end)=[];
zmax3((count3+1):end)=[];


if ltestplot
figure
imagesc(indexmask3)
hold on
for i=1:count3-1
    rectangle('Position',[xmin3(i)-0.5 ymin3(i)-0.5 xmax3(i)-xmin3(i)+1 ymax3(i)-ymin3(i)+1])
end
axis equal tight
hold off
title('after applying rules for radiation in x and y')

end


else
    
   count3 = count2;
   xmin3 = xmin2;
   xmax3 = xmax2;
   ymin3 = ymin2;
   ymax3 = ymax2;
   zmin3 = zmin2;
   zmax3 = zmax2;
    
end

%% aggregate building internal blocks  (only rectangles)
%  1 2 3 4     1 2 3 4
%  5 6 7 8     5i 6 6 7
%  9 a b c  -> 8 6 6 9
%  d e f g     a b c d

%use that all neighbours have the same size in 1 dimension % tg3315 mot
%true now...

datamean4=zeros(size(topomask));
indexmask4=zeros(size(topomask));
internalmask=zeros(size(topomask));


for i=1:count3 
    indexmask4(ymin3(i):ymax3(i),xmin3(i):xmax3(i))=i;
    datamean4(ymin3(i):ymax3(i),xmin3(i):xmax3(i))=zmax3(i);
end
if ltestplot
figure
imagesc(datamean4)
end

for k=1:count3
    i=indexmask4(ymin3(k),xmin3(k)); %index of this block (k = i ?might be unnecessary)
    %check if it is internal
    xind=xmin3(i):xmax3(i);
    yind=ymin3(i):ymax3(i);
    
    temp1=zeros(length(yind),2);
    temp2=zeros(2,length(xind));
    if xmin3(i)==1 %at edge
    temp1(:,1)=999999;    %dummy value
    else
    temp1(:,1)=topo(yind,xmin3(i))-topo(yind,xmin3(i)-1); % tg3315 changed so do not get internal blocks with adjacent different sizes %indexmask4(yind, xmin3(i)-1);    
    end
    if xmax3(i)==ni %at edge
    temp1(:,2)=999999;     
    else
    temp1(:,2)=topo(yind,xmax3(i))-topo(yind,xmax3(i)+1); %indexmask4(yind, xmax3(i)+1);       
    end
    if ymin3(i)==1 %at edge
    temp2(1,:) = 999999;    %dummy value   
    else 
    temp2(1,:) = topo(ymin3(i),xind)-topo(ymin3(i)-1,xind); %indexmask4(ymin3(i)-1,xind)    ;
    end
    if ymax3(i)==nj %at edge
    temp2(2,:) = 999999;    %dummy value  
    else 
    temp2(2,:) = topo(ymax3(i),xind)-topo(ymax3(i)+1,xind); %indexmask4(ymax3(i)+1,xind)  ; 
    end
    
    temp=[temp1(:)' temp2(:)'];
    if all(temp==0) % tg3315 switched this indexing around... %temp>0) %it's internal
        internalmask(yind,xind)=1;
    end
end
externalmask=topomask-internalmask;

if ltestplot
figure
imagesc(indexmask4.*internalmask)
hold on
for i=1:count3-1
    rectangle('Position',[xmin3(i)-0.5 ymin3(i)-0.5 xmax3(i)-xmin3(i)+1 ymax3(i)-ymin3(i)+1])
end
axis equal tight
hold off
title('internal blocks after applying rules for radiation in x and y')
end

% tg3315 may need to develop this further for geometries with VERY uneven
% roof tops... Currently will not allow internal blocks below but this does
% satisy the radiation laws in ther vertical direction...

% merge internal blocks
count5=count3;
xmin5=xmin3;
xmax5=xmax3;
ymin5=ymin3;
ymax5=ymax3;
zmax5=zmax3;
zmin5=zmin3;
indexmask5=indexmask4;
change=true;
%in y first
while change
    change=false;
    for k=1:count5
        i=indexmask5(ymin5(k),xmin5(k));
        if internalmask(ymin5(i),xmin5(i))
            if internalmask(ymin5(i)-1,xmin5(i))  %can only have 1 neighbour in y direction at this stage
                i2=indexmask5(ymin5(i)-1,xmin5(i));
                
                ymin5(i)=ymin5(i2);
                
                indexmask5(ymin5(i):ymax5(i),xmin5(i):xmax5(i))=i;
                indexmask5(indexmask5>i2)=indexmask5(indexmask5>i2)-1;
                
                count5=count5-1;
                xmin5(i2)=[];
                xmax5(i2)=[];
                ymin5(i2)=[];
                ymax5(i2)=[];
                zmax5(i2)=[];
                zmin5(i2)=[];
                
                change=true;
                %count5
                break               
            end
        end
        
    end
end
disp(['Nr of blocks after merging internal blocks along y: ' num2str(length(xmin5))])

if ltestplot
figure
imagesc(indexmask5.*internalmask)
hold on
for i=1:count5
    rectangle('Position',[xmin5(i)-0.5 ymin5(i)-0.5 xmax5(i)-xmin5(i)+1 ymax5(i)-ymin5(i)+1])
end
axis equal tight
hold off
title('internal blocks after merging internal blocks')
end

%make blocks

datamean5=zeros(size(topomask));
dataind=zeros(size(topomask));
for i=1:count5
    datamean5(ymin5(i):ymax5(i),xmin5(i):xmax5(i))=zmax5(i);
     dataind(ymin5(i):ymax5(i),xmin5(i):xmax5(i))=i;
end



if ltestplot
figure
subplot(1,3,1)
imagesc(indexmask5)
hold on
for i=1:count5
    rectangle('Position',[xmin5(i)-0.5 ymin5(i)-0.5 xmax5(i)-xmin5(i)+1 ymax5(i)-ymin5(i)+1])
end
axis equal tight
hold off
title('index: blocks after merging merging internal blocks')
set(gca,'YDir','normal')



subplot(1,3,2)
imagesc(datamean5)
axis equal tight
hold off
title('height: blocks after merging merging internal blocks')
set(gca,'YDir','normal')


subplot(1,3,3)
imagesc(dataind)
axis equal tight
hold off
title('index: blocks after merging merging internal blocks')
set(gca,'YDir','normal')
end


 ymin6=ymin5;
 ymax6=ymax5;
 xmin6=xmin5;
 xmax6=xmax5;
 zmin6=zmin5;
 zmax6=zmax5;


%%
if lhqplot
cd(outputdir)
h=figure;
set(gcf,'units','centimeters','position',[0 0 14.5 14.5]);
set(h,'PaperPosition',[0 0 14.5 14.5]);
set(h,'PaperUnits','centimeters');
imagesc(xf,yf,flipud(datamean2))
set(gca,'YDir','normal','TickLabelInterpreter','latex')
h1=gca;
h1.Position=[0.08 0.1100 0.78 0.8150];
xlim([xh(1) xh(end)])
ylim([yh(1) yh(end)])
caxis([0 1])
xlabel('x [m]','Interpreter','latex','FontSize',12)
ylabel('y [m]','Interpreter','latex','FontSize',12)

colormap(flipud(bone)) %colormap(flipud(gray))
hcb=colorbar('Position',[0.92 0.15 0.03 0.69]);
hcb.Label.Interpreter='latex';
hcb.TickLabelInterpreter='latex';
title(hcb,'height [m]','Interpreter','latex','FontSize',12)

hold on
for i=1:count5-1
    rectangle('Position',[2*(xmin6(i))-4 960-2*(ymin6(i))+4-2*(ymax6(i)-ymin6(i)+1) 2*(xmax6(i)-xmin6(i)+1) 2*(ymax6(i)-ymin6(i)+1)])
end
axis equal tight

hold off

print -depsc2 topo.eps
print -dpng topo.png
cd(parentdir)
end

%% flip the whole y-dimension!!

ymin6f=nj-ymax6+1;
ymax6f=nj-ymin6+1;
ymin2f=nj-ymax2+1;
ymax2f=nj-ymin2+1;

%%
% wtmean5=ones(5,length(xmin6));
% 
% fileID = fopen([tempdir '/blocksmeannointernalflip.inp'],'w');
% fprintf(fileID,'# Block data\n');
% fprintf(fileID,'# block indices,                     wall type\n');
% fprintf(fileID,'#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t ttop\t twest\t teast\t tnor\t tsou\n');
% fprintf(fileID,'%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n',[xmin6';xmax6';ymin6f';ymax6f';zmin6';zmax6';wtmean5]);
% fclose(fileID);
%%
% wtmean=ones(5,count2);
% 
% fileID = fopen([tempdir '/bbri.inp'],'w');
% fprintf(fileID,'# Block data\n');
% fprintf(fileID,'# block indices,                     wall type\n');
% fprintf(fileID,'#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t ttop\t twest\t teast\t tnor\t tsou\n');
% fprintf(fileID,'%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n',[xmin2';xmax2';ymin2f';ymax2f';zmin2';zmax2';wtmean]);
% fclose(fileID);



%write block file to extend with roads and bounding walls etc

dummy=ones(count5,5);

%fileID = fopen([tempdir '/blocks.inp.' num2str(expnr)],'w');
%fprintf(fileID,'# Block data\n');
%fprintf(fileID,'#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t dtop\t dwest\t deast\t dnor\t dsou\n');
%fprintf(fileID,'%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n',[xmin5';xmax5';ymin6';ymax6';zmin5';zmax5';dummy]);
%fclose(fileID);

obj.blocks = [xmin5 xmax5 ymin6 ymax6 zmin5 zmax5 dummy];

%write block file to use for test of ray intersection (without roads that will be added to blocks.inp)
%already raise by 1 (because of roads, is done to blocks.inp.xxx later)
%give each building a unique number, might differ from results in block2fac
%but it doesn't matter

%if source == 1 
%if (size(B,1)<size(xmin5,1))
%whichsource=1;
%nbl=size(B,1);
%blockindexmask=topoind;
%else
blockindexmask=dataind;   
nbl=size(xmin5,1);
%whichsource=2;
%end
%else
%blockindexmask=dataind;
%nbl=size(xmin5,1);
%whichsource=2;
%end
buildingindexmask=bwlabel(blockindexmask);
buildingindexlist=zeros(1,nbl);

if ltestplot
figure
imagesc(buildingindexmask)
axis equal tight
hold off
title('building indeces')
set(gca,'YDir','normal')
end

%if whichsource == 1
%    xuu=B(:,2);
%    yuu=B(:,4);
%    zuu=B(:,6);
%elseif whichsource == 2
    xuu=xmax5;
    yuu=ymax6;
    zuu=zmax5;
%end
        
for i=1:nbl
buildingindexlist(i)=buildingindexmask(yuu(i),xuu(i));    
end

%fileID = fopen([tempdir '/bbri.inp'],'w');
%fprintf(fileID,'# Block data for ray intersection\n');
%fprintf(fileID,'#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t buildingindex\n');
%%if whichsource == 2
% fprintf(fileID,'%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n',[xmin5';xmax5';ymin6';ymax6';zmin5'+1;zmax5'+1;buildingindexlist]);
%%elseif whichsource == 1
%%fprintf(fileID,'%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n',[B(:,1)';B(:,2)';B(:,3)';B(:,4)';B(:,5)';B(:,6)';buildingindexlist]);    
%%end
%fclose(fileID);

bbri = [xmin5  xmax5  ymin6 ymax6 zmin5+1 zmax5+1 buildingindexlist'];
