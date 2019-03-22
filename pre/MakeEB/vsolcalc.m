%clear all
%close all
%feature('numcores')
%delete(gcp('nocreate'))
% if isempty(isempty(gcp))
%    pool=parpool
% end


%% read files

%blocks
nheader=3;
try %in case file is empty -> no blocks
BB = dlmread([outputdir '/blocks.inp.' num2str(expnr)],'',nheader,0);  %#il   iu   jl    ju   kl   ku  ttop  twest teast tnor tsou
catch
BB =[];
end

%blocks for intersection
nheader=3;
try %in case file is empty -> no blocks
B = dlmread([tempdir '/bbri.inp'],'',nheader,0);  %#il   iu   jl    ju   kl   ku  ttop  twest teast tnor tsou
catch
B =[];
end
%floors
nheader=3;
%G = dlmread(['floors.inp.' num2str(expnr)],'',nheader,0);  %#il   iu   jl    ju  type
G = dlmread([outputdir '/floors.txt'],'',nheader,0);  %#il   iu   jl    ju  type
%bounding walls
nheader=3;
W = dlmread([outputdir '/boundingwalls.txt'],'',nheader,0);  %#il   iu   jl    ju  type
%facets
nheader=1;
F = dlmread([outputdir '/facets.inp.' num2str(expnr)],'',nheader,0);  %#   or     wl    blk    bld


[nblocks, nbi]=size(B);
[nfcts, nfacprop]=size(F);
[nfl, nflprop]=size(G);
[nbw, nbwprop]=size(W);

%% assigne azimuthal angle based on orientation
wallaz=zeros(nfcts,1); %wall azimuthal angles

for i=1:nfcts
    if F(i,1)==1 %horizontal surface
        wallaz(i)=solaz; %azimuthal angle irrelevant for horizontal surface. Just set it to "solaz" so it passes self-shading test.
    elseif F(i,1)==2 % west
        wallaz(i)=90;
    elseif F(i,1)==3 % east
        wallaz(i)=270;
    elseif F(i,1)==4 % north
        wallaz(i)=0;
    elseif F(i,1)==5
        wallaz(i)=180;  %=5, south
    end
end
if lwritefiles
dlmwrite([outputdir '/wallaz.txt'],wallaz,'delimiter',' ','precision','%4f')
end

%% vector to sun
x = sin( Z/360*2*pi ) * cos( (solaz+90)/360*2*pi );  %add 90 since it's from north, not from east (= our x coordinate)
y = sin( Z/360*2*pi ) * sin( (solaz+90)/360*2*pi );
z = cos( Z/360*2*pi );
v1=[x,y,z];

%% create blocks to test for intersection
%
bl=zeros(nblocks+nbw,6);
for k=1:nblocks
    xl=xb(B(k,1));
    xu=xb(B(k,2)+1);
    yl=yb(B(k,3));
    yu=yb(B(k,4)+1);
    zl=zb(B(k,5)+1);
    zu=zb(B(k,6)+2);
    bl(k,:)=[xl xu yl yu zl zu];
end
for k=1:nbw
    fi3=F(k+(nfcts-nfl-nbw),1);
    il=W(k,1);
    iu=W(k,2);
    jl=W(k,3);
    ju=W(k,4);
    kl=W(k,5)+1;
    ku=W(k,6)+1;
    if (fi3==2)
        xl=xb(end);
        xu=xb(end)+0.1;
        yl=yb(jl);
        yu=yb(ju+1);
        zl=zb(kl);
        zu=zb(ku+1);
    elseif (fi3==3)
        xl=xb(1)-0.1;
        xu=xb(1);
        yl=yb(jl);
        yu=yb(ju+1);
        zl=zb(kl);
        zu=zb(ku+1);
    elseif (fi3==4)
        xl=xb(il);
        xu=xb(iu+1);
        yl=yb(1)-0.1;
        yu=yb(1);
        zl=zb(kl);
        zu=zb(ku+1);
    else %if (fi==5)
        xl=xb(il);
        xu=xb(iu+1);
        yl=yb(end);
        yu=yb(end)+0.1;
        zl=zb(kl);
        zu=zb(ku+1);
    end
    bl(k+nblocks,:)=[xl xu yl yu zl zu];
end
gl=zeros(nfl,6); %[xl xu yl yu zl zu] %space coordinates of floors
for j=1:nfl
    xl=xb(G(j,1));
    xu=xb(G(j,2)+1);
    yl=yb(G(j,3));
    yu=yb(G(j,4)+1);
    zl=0-0.1;
    zu=0;
    gl(j,:)=[xl xu yl yu zl zu];
end

toc
disp('done creating blocks to test for intersection')

disp('Calculate Ray-Block-Intersection')
%% Calculate Ray-Block-Intersection
%% do it in physical space (not index)
wsl=ones(nfcts,1); %potentially sunlit surface (not self shaded)
asl=zeros(nfcts,1); %sunlit area of facet
walltheta=zeros(nfcts,1); %theta according to URBAN4


for i=1:nfcts
    if abs(wallaz(i)-solaz)>=90  %angle between sun and facet-normal larger 90
        wsl(i)=0; %self shading
        continue
    end


    % get facet center and corners
    [ndim,~,co]=detsub(i,F,BB,G,W,cornm,xb,yb,zb,delta); %determine corners
    %count how much area cannot see sun (i.e. view is blocked)
    [ as ] = prblckd(i,-999,co,ndim,true,v1,-999,-999,F,centerweight,cornerweight,nblocks,nbw,bl); 
 
    asl(i)=1-as; %fraction of facet in sunlight
    
end



%%
toc
disp('done calculating Ray-Block-Intersection')
disp('assign wall theta')

for i=1:nfcts
    if wsl(i)==0 %shaded
        walltheta(i)=90;
    elseif F(i,1)==1 %horizontal
        walltheta(i)=0;
    else
        walltheta(i)=abs(solaz-wallaz(i));
    end
end
toc
disp('done assigning wall theta')


%write them
%dlmwrite('wss.txt',wss,'delimiter',' ','precision','%2d')
if lwritefiles
dlmwrite([outputdir '/walltheta.txt'],walltheta,'delimiter',' ','precision','%2d')


fname = [outputdir '/vsol.txt'];
fileID = fopen(fname,'w');
fprintf(fileID,'# %4s\n','sun lit fraction of facet');
fprintf(fileID,'%6d\n', asl);
fclose(fileID);
end
%% shortwave
%

disp('calculate incoming shortwave ')

Sdir=zeros(nfcts,1); %direct solar radiation on every facet

for i=1:nfcts  %direct from sky
    if F(i,1)==1 %horizontal
        phi=0;
    else
        phi=90;
    end
    Sdir(i)=I*cos((Z-phi)/360*2*pi)*cos(walltheta(i)/360*2*pi)*asl(i);
end
if lwritefiles
fname = [outputdir '/Sdir.txt'];
fileID = fopen(fname,'w');
fprintf(fileID,'# %4s\n','Direct shortwave radiation reaching facets [W/m2]');
fprintf(fileID,'%6d\n', Sdir);
fclose(fileID);
end
toc
disp('done calculating incoming shortwave ')

%%

if ltestplot == 500
  

scale=2;
scalef=1.5;
h=figure;
set(gcf,'units','centimeters','position',[2 2 2+14.5*scale 2+14.5*scale]);
set(h,'PaperPosition',[2 2 2+14.5*scale 2+14.5*scale]);
set(h,'PaperUnits','centimeters');
set(h,'renderer','painters');

    for i=1:nfcts
        bi=F(i,3); %block index
        fi=F(i,1); %facet index
        
        if (F(i,4)==-1) %it is a floor, not a building
            il=G(bi,1);
            iu=G(bi,2);
            jl=G(bi,3);
            ju=G(bi,4);
            xl=xc(il)-dx(il)/2;
            xu=xc(iu)+dx(iu)/2;
            yl=yc(jl)-dy(jl)/2;
            yu=yc(ju)+dy(ju)/2;
            zl=0;
            zu=0;
        elseif (F(i,4)==-99)  %it is a bounding wall
            il=W(bi,1);
            iu=W(bi,2);
            jl=W(bi,3);
            ju=W(bi,4);
            kl=W(bi,5)+1;
            ku=W(bi,6)+1;
            
            if (fi==2)
                xl=xc(iu)+dx(iu)/2;
                xu=xc(iu)+dx(iu);
                yl=yc(jl)-dy(jl)/2;
                yu=yc(ju)+dy(ju)/2;
                zl=zc(kl)-dz(kl)/2;
                zu=zc(ku)+dz(ku)/2;
            elseif (fi==3)
                xl=xc(il)-dx(il);
                xu=xc(il)-dx(il)/2;
                yl=yc(jl)-dy(jl)/2;
                yu=yc(ju)+dy(ju)/2;
                zl=zc(kl)-dz(kl)/2;
                zu=zc(ku)+dz(ku)/2;
            elseif (fi==4)
                xl=xc(il)-dx(il)/2;
                xu=xc(iu)+dx(iu)/2;
                yl=yc(jl)-dy(jl);
                yu=yc(jl)-dy(jl)/2;
                zl=zc(kl)-dz(kl)/2;
                zu=zc(ku)+dz(ku)/2;
            else %if (fi==5)
                xl=xc(il)-dx(il)/2;
                xu=xc(iu)+dx(iu)/2;
                yl=yc(ju)+dy(ju)/2;
                yu=yc(ju)+dy(ju);
                zl=zc(kl)-dz(kl)/2;
                zu=zc(ku)+dz(ku)/2;
            end
        else
            il=B(bi,1);
            iu=B(bi,2);
            jl=B(bi,3);
            ju=B(bi,4);
            kl=B(bi,5)+1;
            ku=B(bi,6)+1;
            
            xl=xc(il)-dx(il)/2;
            xu=xc(iu)+dx(iu)/2;
            yl=yc(jl)-dy(jl)/2;
            yu=yc(ju)+dy(ju)/2;
            zl=zc(kl)-dz(kl)/2;
            zu=zc(ku)+dz(ku)/2;
        end
        
        switch F(i, 1)
            case {1} %top
                x = [xl, xu, xu, xl];
                y = [yl, yl, yu, yu];
                z = [zu, zu, zu, zu];
            case {2} %west
                x = [xl, xl, xl, xl];
                y = [yl, yl, yu, yu];
                z = [zl, zu, zu, zl];
            case {3} %east
                x = [xu, xu, xu, xu];
                y = [yl, yl, yu, yu];
                z = [zl, zu, zu, zl];
            case {4} %north
                x = [xl, xl, xu, xu];
                y = [yu, yu, yu, yu];
                z = [zl, zu, zu, zl];
            case {5} %south
                x = [xl, xl, xu, xu];
                y = [yl, yl, yl, yl];
                z = [zl, zu, zu, zl];
        end
        if (i==99)
        x=[8, 9, 10 , 8];
        y=[8,  8, 9 , 9];
        elseif (i==103)
        x=[9, 10, 10,9];
        y=[7, 7, 9, 8];
        end
        %shaded vs non-shaded
        cs=asl(i)*[0.85 0.85 0.85] + [0.15 0.15 0.15];
        
        patch(x,y,z, cs,'FaceLighting','none');
        hold on
    end

view(15,33)
axis equal 
xlim([0 xb(end)])
ylim([0 yb(end)])
%IC25
%zlim([0 72])  
zlim([0 100])
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12*scalef)
xlabel('x [m]','Interpreter','latex','FontSize',12*scalef)
ylabel('y [m]','Interpreter','latex','FontSize',12*scalef)
zlabel('z [m]','Interpreter','latex','FontSize',12*scalef)
set(gcf, 'Color', 'w');
% export_fig shading.eps
end

%%
%cs
if lhqplot
scale=2;
scalef=1.5;
h=figure;
set(gcf,'units','centimeters','position',[2 2 2+14.5*scale 2+14.5*scale/2]);
set(h,'PaperPosition',[2 2 2+14.5*scale 2+14.5*scale/2]);
set(h,'PaperUnits','centimeters');
set(h,'renderer','painters');

plotshadowline

for i=1:nfcts
        bi=F(i,3); %block index
        fi=F(i,1); %facet index
        
        if (F(i,4)<0 && F(i,4)>-100) %it is a floor, not a building
            il=G(bi,1);
            iu=G(bi,2);
            jl=G(bi,3);
            ju=G(bi,4);
            xl=xc(il)-dx(il)/2;
            xu=xc(iu)+dx(iu)/2;
            yl=yc(jl)-dy(jl)/2;
            yu=yc(ju)+dy(ju)/2;
            zl=0;
            zu=0;
        elseif (F(i,4)<=-100)  %it is a bounding wall
            il=W(bi,1);
            iu=W(bi,2);
            jl=W(bi,3);
            ju=W(bi,4);
            kl=W(bi,5)+1;
            ku=W(bi,6)+1;
            
            if (fi==2)
                xl=xc(iu)+dx(iu)/2;
                xu=xc(iu)+dx(iu);
                yl=yc(jl)-dy(jl)/2;
                yu=yc(ju)+dy(ju)/2;
                zl=zc(kl)-dz(kl)/2;
                zu=zc(ku)+dz(ku)/2;
            elseif (fi==3)
                xl=xc(il)-dx(il);
                xu=xc(il)-dx(il)/2;
                yl=yc(jl)-dy(jl)/2;
                yu=yc(ju)+dy(ju)/2;
                zl=zc(kl)-dz(kl)/2;
                zu=zc(ku)+dz(ku)/2;
            elseif (fi==4)
                xl=xc(il)-dx(il)/2;
                xu=xc(iu)+dx(iu)/2;
                yl=yc(jl)-dy(jl);
                yu=yc(jl)-dy(jl)/2;
                zl=zc(kl)-dz(kl)/2;
                zu=zc(ku)+dz(ku)/2;
            else %if (fi==5)
                xl=xc(il)-dx(il)/2;
                xu=xc(iu)+dx(iu)/2;
                yl=yc(ju)+dy(ju)/2;
                yu=yc(ju)+dy(ju);
                zl=zc(kl)-dz(kl)/2;
                zu=zc(ku)+dz(ku)/2;
            end
        else
            il=B(bi,1);
            iu=B(bi,2);
            jl=B(bi,3);
            ju=B(bi,4);
            kl=B(bi,5)+1;
            ku=B(bi,6)+1;
            
            xl=xc(il)-dx(il)/2;
            xu=xc(iu)+dx(iu)/2;
            yl=yc(jl)-dy(jl)/2;
            yu=yc(ju)+dy(ju)/2;
            zl=zc(kl)-dz(kl)/2;
            zu=zc(ku)+dz(ku)/2;
        end
        
        switch F(i, 1)
            case {1} %top
                x = [xl, xu, xu, xl];
                y = [yl, yl, yu, yu];
                z = [zu, zu, zu, zu];
            case {2} %west
                x = [xl, xl, xl, xl];
                y = [yl, yl, yu, yu];
                z = [zl, zu, zu, zl];
            case {3} %east
                x = [xu, xu, xu, xu];
                y = [yl, yl, yu, yu];
                z = [zl, zu, zu, zl];
            case {4} %north
                x = [xl, xl, xu, xu];
                y = [yu, yu, yu, yu];
                z = [zl, zu, zu, zl];
            case {5} %south
                x = [xl, xl, xu, xu];
                y = [yl, yl, yl, yl];
                z = [zl, zu, zu, zl];
        end
        if (i==99)
            x=[9, 10, 11 , 9];
            y=[8,  8, 9 , 9];
        elseif (i==103)
            x=[10, 11, 11, 10];
            y=[7, 7, 9, 8];
        end
        %shaded vs non-shaded
        greymax=0.15; %maximum percent of blackness, [1 1 1] is white;
        cs=asl(i)*[1-greymax 1-greymax 1-greymax] + [greymax greymax greymax];
        patch(x,y,z, cs,'FaceLighting','none');
        hold on
        %quiver3(face(92,1,2,1),face(92,1,2,2),face(92,1,2,3), v1(1),v1(2),v1(3), 8);
        %quiver3(face(92,1,3,1),face(92,1,3,2),face(92,1,3,3), v1(1),v1(2),v1(3), 8);
        % quiver3(face(92,2,2,1),face(92,2,2,2),face(92,2,2,3), v1(1),v1(2),v1(3), 8);
        %  quiver3(face(92,2,3,1),face(92,2,3,2),face(92,2,3,3), v1(1),v1(2),v1(3), 8);
        
        
     %   title('sunlit fraction')
end
HS=scatter3(shadowscatter(:,1),shadowscatter(:,2),shadowscatter(:,3),'filled','dk');

    ffv=134;
    quiver3((xb(fctl(ffv, 6))+xb(fctl(ffv, 7)))/2,(yb(fctl(ffv, 8))+yb(fctl(ffv, 9)))/2,zb(fctl(ffv, 10))+0.05,v1(1),v1(2),v1(3),3,'r','LineWidth',2,'MaxHeadSize',1)
    quiver3(xb(fctl(ffv, 6)),yb(fctl(ffv, 8)),zb(fctl(ffv, 10))+0.05,v1(1),v1(2),v1(3),1.5,'r','LineWidth',2,'MaxHeadSize',1)
    quiver3(xb(fctl(ffv, 7)),yb(fctl(ffv, 8)),zb(fctl(ffv, 10))+0.05,v1(1),v1(2),v1(3),0.5,'r','LineWidth',2,'MaxHeadSize',1)
    quiver3(xb(fctl(ffv, 6)),yb(fctl(ffv, 9)),zb(fctl(ffv, 10))+0.05,v1(1),v1(2),v1(3),3,'r','LineWidth',2,'MaxHeadSize',1)
    quiver3(xb(fctl(ffv, 7)),yb(fctl(ffv, 9)),zb(fctl(ffv, 10))+0.05,v1(1),v1(2),v1(3),3,'r','LineWidth',2,'MaxHeadSize',1)

greycm=zeros(9,3);
cticks=zeros(9,1);
ctickl=zeros(9,1);
for l=0:8
cs=l*0.125*[1-greymax 1-greymax 1-greymax] + [greymax greymax greymax];
cticks(l+1)=0.15+(l)*0.85/9+0.85/18;
ctickl(l+1)=l*0.125;
greycm(l+1,:)=cs;
end
cm=colormap(greycm);
c=colorbar;
c.Title.Interpreter='latex';
set(c, 'TickLabelInterpreter', 'latex')
set(c,'FontSize',12*scalef)
caxis([0.15 1])
c.Ticks=cticks
c.TickDirection='out'
c.TickLabels=num2cell(ctickl)
c.Title.String='$f_{\mathrm{e}}$ $[-]$'
%c.Position(1)=c.Position(1)+0.00
c.Position(2)=c.Position(2)+0.02
c.Position(4)=c.Position(4)-0.2

view(3)
axis equal 
xlim([0 xb(end)])
ylim([0 yb(end)])
% zlim([0 72])
zlim([0 3])  
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12*scalef)
%h1=gca;
% h1.Position=[0.08 0.1100 0.78 0.8150];
%caxis([0 71])
xlabel('x [-]','Interpreter','latex','FontSize',12*scalef)
ylabel('y [-]','Interpreter','latex','FontSize',12*scalef)
zlabel('z [-]','Interpreter','latex','FontSize',12*scalef) 
set(gcf, 'Color', 'w');
export_fig facetshadowalt.eps   


%%
scale=2;
scalef=1.5;
h=figure;
%cm=colormap(MPL_gist_heat(31:end,:));
%cm=colormap(MPL_viridis(30:end,:));
cm=colormap(flipud(perscolormaps('MPL_BuPu')));
cm=colormap(cm(1:end-2,:));

set(gcf,'units','centimeters','position',[2 2 2+14.5*scale 2+14.5*scale/2]);
set(h,'PaperPosition',[2 2 2+14.5*scale 2+14.5*scale/2]);
set(h,'PaperUnits','centimeters');
set(h,'renderer','painters');    
    for i=1:nfcts
        bi=F(i,3); %block index
        fi=F(i,1); %facet index
        
        if (F(i,4)<0 && F(i,4)>-100) %it is a floor, not a building
            il=G(bi,1);
            iu=G(bi,2);
            jl=G(bi,3);
            ju=G(bi,4);
            xl=xc(il)-dx(il)/2;
            xu=xc(iu)+dx(iu)/2;
            yl=yc(jl)-dy(jl)/2;
            yu=yc(ju)+dy(ju)/2;
            zl=0;
            zu=0;
        elseif (F(i,4)<=-100)  %it is a bounding wall
            il=W(bi,1);
            iu=W(bi,2);
            jl=W(bi,3);
            ju=W(bi,4);
            kl=W(bi,5)+1;
            ku=W(bi,6)+1;
            
            if (fi==2)
                xl=xc(iu)+dx(iu)/2;
                xu=xc(iu)+dx(iu);
                yl=yc(jl)-dy(jl)/2;
                yu=yc(ju)+dy(ju)/2;
                zl=zc(kl)-dz(kl)/2;
                zu=zc(ku)+dz(ku)/2;
            elseif (fi==3)
                xl=xc(il)-dx(il);
                xu=xc(il)-dx(il)/2;
                yl=yc(jl)-dy(jl)/2;
                yu=yc(ju)+dy(ju)/2;
                zl=zc(kl)-dz(kl)/2;
                zu=zc(ku)+dz(ku)/2;
            elseif (fi==4)
                xl=xc(il)-dx(il)/2;
                xu=xc(iu)+dx(iu)/2;
                yl=yc(jl)-dy(jl);
                yu=yc(jl)-dy(jl)/2;
                zl=zc(kl)-dz(kl)/2;
                zu=zc(ku)+dz(ku)/2;
            else %if (fi==5)
                xl=xc(il)-dx(il)/2;
                xu=xc(iu)+dx(iu)/2;
                yl=yc(ju)+dy(ju)/2;
                yu=yc(ju)+dy(ju);
                zl=zc(kl)-dz(kl)/2;
                zu=zc(ku)+dz(ku)/2;
            end
        else
            il=B(bi,1);
            iu=B(bi,2);
            jl=B(bi,3);
            ju=B(bi,4);
            kl=B(bi,5)+1;
            ku=B(bi,6)+1;
            
            xl=xc(il)-dx(il)/2;
            xu=xc(iu)+dx(iu)/2;
            yl=yc(jl)-dy(jl)/2;
            yu=yc(ju)+dy(ju)/2;
            zl=zc(kl)-dz(kl)/2;
            zu=zc(ku)+dz(ku)/2;
        end
        
        switch F(i, 1)
            case {1} %top
                x = [xl, xu, xu, xl];
                y = [yl, yl, yu, yu];
                z = [zu, zu, zu, zu];
            case {2} %west
                x = [xl, xl, xl, xl];
                y = [yl, yl, yu, yu];
                z = [zl, zu, zu, zl];
            case {3} %east
                x = [xu, xu, xu, xu];
                y = [yl, yl, yu, yu];
                z = [zl, zu, zu, zl];
            case {4} %north
                x = [xl, xl, xu, xu];
                y = [yu, yu, yu, yu];
                z = [zl, zu, zu, zl];
            case {5} %south
                x = [xl, xl, xu, xu];
                y = [yl, yl, yl, yl];
                z = [zl, zu, zu, zl];
        end
        if (i==99)
            x=[9, 10, 11 , 9];
            y=[8,  8, 9 , 9];
        elseif (i==103)
            x=[10, 11, 11, 10];
            y=[7, 7, 9, 8];
        end
        %shaded vs non-shaded
        %cs=Sdir(i)/I*[0.85 0.85 0.85] + [0.15 0.15 0.15];
        cs=cm(max(1,min(size(cm,1),round(Sdir(i)/I*size(cm,1)))),:);
        
        patch(x,y,z, cs,'FaceLighting','none');
        hold on
     %        quiver3(face(92,1,2,1),face(92,1,2,2),face(92,1,2,3), v1(1),v1(2),v1(3), 1);
     %        quiver3(face(92,1,3,1),face(92,1,3,2),face(92,1,3,3), v1(1),v1(2),v1(3), 1);
     %        quiver3(face(92,2,2,1),face(92,2,2,2),face(92,2,2,3), v1(1),v1(2),v1(3), 1);
     %        quiver3(face(92,2,3,1),face(92,2,3,2),face(92,2,3,3), v1(1),v1(2),v1(3), 1);
     %   title('Sdir')
    end
c=colorbar;
c.Title.Interpreter='latex';
set(c, 'TickLabelInterpreter', 'latex')
set(c,'FontSize',12*scalef)
c.Title.String='$S$ $[\mathrm{W}\mathrm{m}^{-2}]$'
%c.Position(1)=c.Position(1)+0.00
c.Position(2)=c.Position(2)+0.02
c.Position(4)=c.Position(4)-0.2


caxis([0 I])
view(3)
axis equal 
xlim([0 xb(end)])
ylim([0 yb(end)])
% zlim([0 72])
zlim([0 3])  
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',12*scalef)
%h1=gca;
% h1.Position=[0.08 0.1100 0.78 0.8150];
%caxis([0 71])
xlabel('x [-]','Interpreter','latex','FontSize',12*scalef)
ylabel('y [-]','Interpreter','latex','FontSize',12*scalef)
zlabel('z [-]','Interpreter','latex','FontSize',12*scalef) 
set(gcf, 'Color', 'w');
%export_fig facetSdir.eps    
    
end



