%% creates blocks.inp from greyscale-image 

%% initialising

A = imread(['lidar_data/' sourcename]);  %read topo image
[njorig, niorig, ~]=size(A);


% since its a greyscale image all 3 rgb channels have the same value
% substract value from 255 and scale by max value to get topography
topot=(255-double(A(:,:,1)))/255*maxh;

if ltestplot
figure; imagesc(topot) 
end

% interpolate to a coarser grid if necessary
if dx~=dxinp || dy~=dyinp %can't use imresize, have to do it manually since x and y scaling might be different
xxxo=dxinp/2:dxinp:(niorig*dxinp); % original pixel centres
yyyo=dyinp/2:dyinp:(njorig*dyinp);
xxx=dx/2:dx:(niorig*dxinp); %desired pixel centres
yyy=dy/2:dy:(njorig*dxinp);
[X,Y]=meshgrid(xxx,yyy);
topot=interp2(yyyo,xxxo,topot',Y,X,'nearest');
clear xxxo yyyo X Y
end

if ltestplot
figure; imagesc(topot);
end


% upper and lower coordinates to select
ip=round(centeri/dx)+ni/2-1-pad;
im=round(centeri/dx)-ni/2+pad;
jp=round(centerj/dy)+nj/2-1-pad;
jm=round(centerj/dy)-nj/2+pad;

topo=zeros(nj,ni);
topo(pad+(1:(nj-2*pad)),pad+(1:(ni-2*pad)))=topot(jm:jp,im:ip);

if ltestplot
figure; imagesc(xf,yf,topo) 
end

% round height to grid
topo=round(topo/dz)*dz;
topoinit=topo;


%% manual mainpulation of the topograpy come here, if necessary
%test for varying building height
% topo(170:180,173:209)=topo(170:180,173:209)+10;
% topo(180:184,291:331)=topo(180:184,291:331)+12.5;
% topo(143:148,173:214)=topo(143:148,173:214)+7.5;


%% processing
%% fill big holes
topoh=imfill(topo,'holes');
to=topoh-topo;
toi=imbinarize(to);

toi=bwareaopen(toi,smallarea);
topo=topoh-toi.*to;
clear to toi topoh topot

if ltestplot
figure; imagesc(xf,yf,topo); title('after filling holes')
end

% create building mask
topomask=imbinarize(topo);

%% remove small objects
topomask=bwareaopen(topomask,smallarea);

topo=topo.*topomask;

if ltestplot
    figure
    imagesc(xf,yf,topo)
    title('after removing small objects')
    xlim([xh(1) xh(end)])
    ylim([yh(1) yh(end)])
end

%% fill 1 cell gaps  %potentially problematic if the gap is between blocks of different height
[topomask, topo] = fillgaps(topomask,topo,nj,ni,pad);
if ltestplot
    figure
    imagesc(xf,yf,topo)
    title('after filling gaps (size 1)')
    xlim([xh(1) xh(end)])
    ylim([yh(1) yh(end)])
end
if pad==0  %buildings to the edge, make sure there is no weird gaps at domain edge
[topomask, topo] = smoothborders(topomask,topo,nj,ni,pad);
end
if ltestplot
    figure
    imagesc(xf,yf,topo)
    title('after smoothing borders')
    xlim([xh(1) xh(end)])
    ylim([yh(1) yh(end)])
end

%% remove cells with only 1 neighbour
data=topo;
datamask=topomask;

c2=1;
while c2~=0  %repeat until there is no cell with only one neighbour left
    c2=0;
    for i=2:ni-1
        for j=2:nj-1
            count=0;
            if datamask(j,i)>0
                c=datamask(j,i-1)+datamask(j,i+1)+datamask(j+1,i)+datamask(j-1,i);
                if c>1
                    continue
                else
                    c2=c2+1;
                    datamask(j,i)=0;
                    data(j,i)=0;
                end
            end
        end
    end
end

topo=data;
topomask=datamask;

if ltestplot
    figure
    title('after removing cells with only 1 neighbour')
    imagesc(xf,yf,topo)
    xlim([xh(1) xh(end)])
    ylim([yh(1) yh(end)])
end
%%
if lhqplot
    cd(outputdir)
  
    h=figure;
    hp1=subplot(1,2,1);
    imagesc(xf,yf,flipud(topoinit))
    set(gca,'YDir','normal','TickLabelInterpreter','latex')
    xlim([xh(1) xh(end)])
    ylim([yh(1) yh(end)])
    caxis([0 45])
    axis equal tight
    xlabel('x [m]','Interpreter','latex','FontSize',12)
    ylabel('y [m]','Interpreter','latex','FontSize',12)
    hp1.XTick=[0 200 400 600 800 960];
    hp1.YTick=[0 200 400 480];
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    %set(gca,XTick,[0 200 400 600 800 960])
    %set(gca,YTick,[0 200 400 480])
    %set(gca,'YTickLabel',[]);
    %set(gca,'XTickLabel',[]);
    
    hp2=subplot(1,2,2);
    imagesc(xf,yf,flipud(topo))
    set(gca,'YDir','normal','TickLabelInterpreter','latex')
    xlim([xh(1) xh(end)])
    ylim([yh(1) yh(end)])
    caxis([0 45])
    axis equal tight
    xlabel('x [m]','Interpreter','latex','FontSize',12)
    set(gca,'YTickLabel',[]);
    hp2.XTick=[0 200 400 600 800 960];
    hp1.Position=[0.09 0.1100 0.35 0.8150];
    hp2.Position=[0.511 0.1100 0.35 0.8150];
    hcb=colorbar('Position',[0.92 0.425 0.03 0.15]);
    colormap(flipud(bone)) %colormap(flipud(gray))
    hcb.Label.Interpreter='latex';
    hcb.TickLabelInterpreter='latex';
    title(hcb,'height [m]','Interpreter','latex','FontSize',12)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    set(gcf, 'Color', 'w');
    
    %ylabel('y [m]','Interpreter','latex','FontSize',12)
    %
    export_fig topocompbone -eps -png
    %print -depsc2 ICtopocompbone.eps
    %print -dpng ICtopocompbone.png
    %print -dpdf ICtopocompbone.pdf
    
    %%
    h=figure;
    set(gcf,'units','centimeters','position',[0 0 14.5 14.5]);
    set(h,'PaperPosition',[0 0 14.5 14.5]);
    set(h,'PaperUnits','centimeters');

    imagesc(xf,yf,flipud(topoinit))
    set(gca,'YDir','normal','TickLabelInterpreter','latex')
    xlim([xh(1) xh(end)])
    ylim([yh(1) yh(end)])
    caxis([0 45])
    axis equal tight
    xlabel('','Interpreter','latex','FontSize',10)
    ylabel('y [m]','Interpreter','latex','FontSize',10)
    h.Children.YTick=[0 200 400 480];
    set(gca,'XTickLabel',[]);
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    colormap(flipud(bone)) %colormap(flipud(gray))
    hcb=colorbar('northoutside')
    hcb.Label.Interpreter='latex';
    hcb.TickLabelInterpreter='latex';
    title(hcb,'height [m]','Interpreter','latex','FontSize',10)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',12)
    set(gcf, 'Color', 'w');
     %colormap(flipud(gray))

  export_fig inittopobone -eps -png

    
    h=figure;
    set(gcf,'units','centimeters','position',[0 0 14.5 14.5]);
    set(h,'PaperPosition',[0 0 14.5 14.5]);
    set(h,'PaperUnits','centimeters');
    imagesc(xf,yf,flipud(topo))
    set(gca,'YDir','normal','TickLabelInterpreter','latex')
    xlim([xh(1) xh(end)])
    ylim([yh(1) yh(end)])
    caxis([0 45])
    axis equal tight
    h.Children.XTick=[0 200 400 600 800 960];
    xlabel('x [m]','Interpreter','latex','FontSize',10)
    ylabel('y [m]','Interpreter','latex','FontSize',10)
    colormap(flipud(bone))
    h.Children.YTick=[0 200 400 480];
    %h.Position=[0.09 0.1100 0.35 0.8150];
    %h.Position=[0.511 0.1100 0.35 0.8150];
    %hcb=colorbar('Position',[0.92 0.425 0.03 0.15]);
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',10)
    set(gcf, 'Color', 'w');
    
    %ylabel('y [m]','Interpreter','latex','FontSize',12)
    %
   export_fig topobone -eps -png

 cd(parentdir)   
end


%should probably write block.inp. here