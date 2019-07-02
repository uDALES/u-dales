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
            count=count+1;
            yminl=j;
            
            heightgradient=diff(topo(j:end,i));
            
            if isempty(heightgradient) 
            heightchangey=j;
            elseif all(heightgradient==0)
            heightchangey=nj;    
            else
            heightchangey=find(heightgradient~=0,1)+j-1;  %last cell with same height as j (i.e. there is a height change betwen this and the next cell
            end

            ymaxl=heightchangey;
            
            j=heightchangey+1;
            
            blchecked(yminl:ymaxl,xminl:xmaxl)=blchecked(yminl:ymaxl,xminl:xmaxl)+1;
            indexmask(yminl:ymaxl,xminl:xmaxl)=count;
            xmin(count)=xminl;
            xmax(count)=xmaxl;
            ymin(count)=yminl;
            ymax(count)=ymaxl;
            zmin(count)=0;
            zmax(count)=topo(yminl,xminl)/dz;
            
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
while dsize>0  %do as long as there are unmerged blocks
    i=1;   
    while 1  %do along x
        a=xmax2(i);  %upper x index of this block
        bv=find(xmin2==(a+1));  %all blocks with a lower x bound 1 bigger than this blocks upper x bound
        b2=bv(ymin2(bv)==ymin2(i)); %all of those blocks with also the same lower y bound
        if ~isempty(b2);
            if (ymax2(b2)==ymax2(i)) && (zmax2(b2)==zmax2(i)) %if they also have the same upper y bound and the same height

                xmax2(i)=xmax2(b2); %merge 
                xmax2(b2)=[];ymax2(b2)=[];zmax2(b2)=[];xmin2(b2)=[];ymin2(b2)=[];zmin2(b2)=[]; %remove the just merged block from list of blocks
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
%%

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

indexmask3=indexmask2;
change=true;
count3=count2;
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
        
        %
        %
        
        if (yminr>ymin3(c)) %this block extends further (down) than reight neighbour, this block has to be split
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
            
        elseif (yminl>ymin3(c))  %this block extends further (down) than left neighbour, this block has to be split
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
            
        elseif (ymaxl<ymax3(c))  %this block extends further (up) than left neighbour, this block has to be split
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
            
        elseif(ymaxr<ymax3(c)) %this block extends further (up) than reight neighbour, this block has to be split
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

%% aggregate building internal blocks  (only rectangles)
%  1 2 3 4     1 2 3 4
%  5 6 7 8     5 6 6 7
%  9 a b c  -> 8 6 6 9
%  d e f g     a b c d

%use that all neighbours have the same size in 1 dimension

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
    temp1(:,1)=indexmask4(yind, xmin3(i)-1);    
    end
    if xmax3(i)==ni %at edge
    temp1(:,2)=999999;     
    else
    temp1(:,2)=indexmask4(yind, xmax3(i)+1);       
    end
    if ymin3(i)==1 %at edge
    temp2(1,:) = 999999;    %dummy value   
    else 
    temp2(1,:) = indexmask4(ymin3(i)-1,xind)    ;
    end
    if ymax3(i)==nj %at edge
    temp2(2,:) = 999999;    %dummy value  
    else 
    temp2(2,:) = indexmask4(ymax3(i)+1,xind)  ; 
    end
    
    temp=[temp1(:)' temp2(:)'];
    if all(temp>0) %it's internal
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
caxis([0 maxh])
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


wtmean=ones(5,count5);

fileID = fopen([tempdir '/bbri.inp'],'w');
fprintf(fileID,'# Block data\n');
fprintf(fileID,'# block indices,                     wall type\n');
fprintf(fileID,'#  il\t   iu\t   jl\t   ju\t   kl\t   ku\t ttop\t twest\t teast\t tnor\t tsou\n');
fprintf(fileID,'%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\t%5d\n',[xmin5';xmax5';ymin6';ymax6';zmin5';zmax5';wtmean]);
fclose(fileID);
