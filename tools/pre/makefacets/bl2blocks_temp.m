dx = obj.xsize/obj.imax;
dy = obj.ysize/obj.jtot;
dz = obj.zsize/obj.kmax;
dh = obj.zsize;
ni = obj.imax; 
nj = obj.jtot; 
nk = obj.kmax;
zh = obj.zh;
xh=obj.xh;
yh=obj.yh;
dzf=obj.zh(2:end)-obj.zh(1:end-1);

ltestplot = 0;
lhqplot = 0;
lradiation = obj.lEB;
expnr = obj.expnr;

%maximum size of floors and bounding walls (cells in each dimension)
if lradiation
    maxsize = 10; %ADD TO NAMOPTIONS EVENTUALLY
else
    maxsize = inf;
end

blocks = obj.bl;

topomask = zeros(nj,ni);
topo = zeros(nj,ni);

% if no blocks add lowest level
if isnan(blocks)
else
    for n = 1:size(blocks,1)
        topo(blocks(n,3):blocks(n,4),blocks(n,1):blocks(n,2)) = zh(blocks(n,6)+1);
        topomask(blocks(n,3):blocks(n,4),blocks(n,1):blocks(n,2)) = 1;
    end
end

