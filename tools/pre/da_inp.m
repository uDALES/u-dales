%% Code to produce *.inp.* files for DALES

close all
clear all

%DA_TOPDIR =system('$DA_EXPDIR') %cd('$DA_EXPDIR') %Set directory to create folder with

%% Setup

expnr = '010';
%expnr = num2str(iexpnr);

CPUS = 2;       % # cpus
% # cpus, set from execution

DA_EXPDIR = getenv('DA_EXPDIR');
DA_TOPDIR = getenv('DA_TOPDIR');
DA_PREDIR = getenv('DA_PREDIR');

addpath([DA_PREDIR '/']);
exppath=[DA_EXPDIR '/'];
r=da_pp(expnr,exppath); % Object to hold the simulation parameters

cd([DA_EXPDIR '/' expnr])

%% Domain

da_pp.addvar(r, 'imax',64)                    % # cells in x-direction
da_pp.addvar(r, 'xsize',64)                   % Domain size in x-direction

da_pp.addvar(r, 'jtot',64)                    % # cells in j-direction
da_pp.addvar(r, 'ysize',64)                   % Domain size in y-direction

da_pp.addvar(r, 'kmax',32)                    % # cells in k-direction
da_pp.addvar(r, 'zsize',32)                   % Domain size in z-direction

%% Blocks

da_pp.addvar(r, 'lcastro',0)             % switch for staggered cube pattern
da_pp.addvar(r, 'lcube',0)               % switch for linear cube pattern
da_pp.addvar(r, 'lblocks',0)             % switch for infinite blocks
if (not(r.lcastro) && not(r.lcube) && not(r.lblocks)) % switch
    da_pp.addvar(r, 'lflat',1)
    disp('No standard block config. setup so flat domain assumed.')
else
    da_pp.addvar(r, 'lflat',0)
end

da_pp.addvar(r, 'blockheight',16)        % block height in gridspaces
da_pp.addvar(r, 'blockwidth',16)         % block width in gridspaces
da_pp.addvar(r, 'canyonwidth',16)        % canyonwidth in gridspaces

%% Profiles

% forcing
da_pp.addvar(r, 'lmassflowr',0)    % switch for constant mass flow rate
da_pp.addvar(r, 'lprofforc',0)     % switch for 1D geostrophic forcing
da_pp.addvar(r, 'lcoriol',0)       % switch for coriolis forcing
if (not(r.lmassflowr) && not(r.lprofforc) && not(r.lcoriol)) % switch
    da_pp.addvar(r, 'ldp',1)
    disp('No forcing switch config. setup so initial velocities and pressure gradients applied.')
else
    da_pp.addvar(r, 'ldp',0)
end

da_pp.addvar(r, 'u0',0)             % initial u-velocity - also applied as geostrophic term where applicable
da_pp.addvar(r, 'v0',0)             % initial v-velocity - also applied as geostrophic term where applicable

da_pp.addvar(r, 'dpdx',0)           % dpdx
da_pp.addvar(r, 'dpdy',0)           % dpdy

% temperature
da_pp.addvar(r, 'thl0',288)         % temperature at lowest level
da_pp.addvar(r, 'qt0',0)
da_pp.addvar(r, 'sv10',0)
da_pp.addvar(r, 'sv20',0)
da_pp.addvar(r, 'lapse',0)          % lapse rate Ks-1

% other
da_pp.addvar(r, 'w_s',0)            % subsidence
da_pp.addvar(r, 'R',0)              % radiative forcing

% WALLS
da_pp.addvar(r, 'z0horiz',0.01)
da_pp.addvar(r, 'z0hhoriz',0.000067)
da_pp.addvar(r, 'Thoriz',288)
da_pp.addvar(r, 'Twest',288)
da_pp.addvar(r, 'Teast',288)
da_pp.addvar(r, 'Tnorth',288)
da_pp.addvar(r, 'Tsouth',288)

%% Extras

da_pp.addvar(r, 'ltrees',0)               % switch for trees capability
da_pp.addvar(r, 'lpurif',0)               % switch for purifiers
da_pp.addvar(r, 'lzstretch' ,0)           % switch for stretching z grid

%% Pollutants/ chemistry

da_pp.addvar(r, 'lchem' ,0)               % switch for chemistry
da_pp.addvar(r, 'NOb' ,0)
da_pp.addvar(r, 'NO2b' ,0)
da_pp.addvar(r, 'O3b' ,0)

%% Trees - units in grid size scale from actual size before inputting
if r.ltrees
    
    da_pp.addvar(r, 'tree_dz',0)
    da_pp.addvar(r, 'tree_dx',0)
    da_pp.addvar(r, 'tree_h',0)
    da_pp.addvar(r, 'tree_w',0)
    da_pp.addvar(r, 'tree_b',0)
    
    da_pp.addvar(r, 'nt1',0)
    da_pp.addvar(r, 'md',0)
    da_pp.addvar(r, 'ww',0)
    da_pp.addvar(r, 'lw',0)
    da_pp.addvar(r, 'nt2',0)
    
    if r.lblocks
        tree_dz = r.tree_dz;     % tree starting point from bottom
        tree_dx = r.tree_dx;      % distance from block
        tree_h = r.tree_h;      % tree height
        tree_w = r.tree_w;       % tree width
        tree_b = r.jtot;
        
        nt1 = r.nt1;      %one treesversion
        md = r.md;       %one tree: middle
        ww = r.ww;       %one tree: windward
        lw = r.lw;       %one tree: leeward
        nt2 = r.nt2;      %two treesversion
        
    elseif r.lflat
        tree_h = r.tree_h;
        tree_w = r.tree_w;
        tree_b = r.tree_b;
    end
    
end

%% Purifiers
if r.lpurif
    
    if r.lflat
        
    elseif r.lblocks
        purif_dz = 1;      % purifier starting point from bottom
        purif_dx = 3;      % distance from block
        purif_h = 3;       % purifier height
        purif_w = 0;       % purifier width
        purif_dy = 1;      % depth of purif (in y)
        purif_sp = 31;     % spacing between purifiers
        purif_i = 1;       % case for purifier (1 = +ve x, 2 = -ve x, 3 = +ve y etc.)
        
        npurif = je/(purif_dy+purif_sp);
        
        %error test
        if ceil(npurif) ~= floor(npurif)
            lp = 0:je/2;
            indp = rem(je/2,lp)==0; %// logical index that tells if remainder is zero or not
            errp = ([lp(indp); (je/2)./lp(indp)]);
            disp(['Block system does not fit grid'])
            disp(['sum widths to: ' num2str(errp(1,:))])
            disp(['Current width: ' num2str(purif_dy + purif_sp)])
            error('Incorrect block system')
        else
            disp(['Successful block network'])
            disp(['number of blocks ' num2str(npurif*2+1)])
        end
        
    else
        error('Must be r.lblocks to utilise purifiers')
    end
    
end

%% Domain processing

ib = 1; ie = r.imax; h = r.xsize;
jb = 1; je = r.jtot; hy = r.ysize;
kb = 1; ke = r.kmax; hz = r.zsize;

%% Domain Check j grid vs. CPUS

if ceil(je/CPUS) ~= floor(je/CPUS)
    
    disp(['Possible je: ' num2str([2 3 4 5 6 7 8]*CPUS)])
    error('No. CPUs does not fit j grid size')
    
else
    disp('CPUS and jtot successful')
end

%% Specify variables to be added to files

% forcing
if r.lmassflowr
    pqx = 0.0;
    pqy = 0.0;
    vg(1:ke) = 0;
    ug(1:ke) = 0;
    u(1:ke) = r.u0; %linspace(0,r.u0,ke);
    v(1:ke) = r.v0; %linspace(0,r.v0,ke);
    disp(['forcing : massflowr'])
elseif r.lprofforc
    pqx = 0.0;
    pqy = 0.0;
    vg(1:ke) = r.v0;
    ug(1:ke) = r.u0;
    u(1:ke) = r.u0;
    v(1:ke) = r.v0;
    disp(['forcing : profforc'])
elseif r.lcoriol
    pqx = 0.0;
    pqy = 0.0;
    vg(1:ke) = r.v0;
    ug(1:ke) = r.u0;
    u(1:ke) = r.u0;
    v(1:ke) = r.v0;
    disp(['forcing : r.lcoriol'])
elseif r.ldp
    pqx = r.dpdx;
    pqy = r.dpdy;
    vg(1:ke) = 0;
    ug(1:ke) = 0;
    u(1:ke) = r.u0;
    v(1:ke) = r.v0;
    disp(['forcing : ldp'])
end

if (r.lmassflowr+r.lprofforc+r.lcoriol+r.ldp)>1
    if r.lcoriol
        disp('coriolis')
    end
    if r.lprofforc
        disp('profforc')
    end
    if r.lmassflowr
        disp('massflowr')
    end
    error('More than one forcing specified (listed above)')
end

% other
qt0 = r.qt0;
sv10 = r.sv10;
sv20 = r.sv20;

% temperature^M
thl(1) = r.thl0;
disp(['r.thl0 = ' num2str(r.thl0)])

% r.lapse rate
for k = 1:ke-1
    thl(k+1) = thl(k)+ r.lapse*hz/ke;
end

thl_top = thl(ke)+(thl(ke)-thl(ke-1))/2;
display(['thl_top =' num2str(thl(ke)+(thl(ke)-thl(ke-1)))])

% others
wfls(1:ke) = r.w_s;
dthlrad(1:ke) = r.R;


%% Plots

figure
subplot(141)
plot(thl,1:ke)
title('temp')

subplot(142)
plot(dthlrad,1:ke)
title('radiative forcing')

subplot(143)
plot(wfls,1:ke)
title('subsidence')

subplot(144)
plot(u,1:ke)
hold on
plot(v,1:ke, 'r--')
title('velocity')
legend('u', 'v')


%% Grid stretching

% Non-equidistant grid variables
offset  = 0;
growthx = 1;    % initial value
x_top   = h;
offset  = 0;
growthz = 1;    % initial value
z_top   = hz;

%% xgrid calc
nit      =   100000;     % numbe

dxf      =   zeros(ie-ib+1,1);
xh       =   zeros(ie-ib+2,1);
%dxf(ib) =   0.025/26;  % first cell length
dxf(ib)  =   h/ie;  % first cell length  (from 'main_vertical_0toh_48.m)
xh(ib)   =   0;
xh(ib+1) =   dxf(ib); % first grid point

for n=1:nit
    
    % determine x-grid
    for i=ib+1:ie
        dxf(i)=dxf(i-1)*growthx;
        xh(i+1)=xh(i)+dxf(i);
    end
    
    if (xh(end) < x_top-0.00000000000000005 || xh(end) > x_top+0.00000000000000005)
        
        %growthx=growthx*(x_top/(xh(end)))^0.005;
        
    else
        display('converged!')
        display('xf(end)+0.5*dxf(end)= ')
        break
    end
end

xf = zeros(ie-ib+1,1);

for i=ib:ie
    xf(i) = 0.5*(xh(i+1)+xh(i));
end

% add the offset (top of the block)
xf  = xf  + offset; %, 15;
xh  = xh  + offset;

% figure(1)
% plot(xf,'-*')
% hold on
% plot(xh,'-og')xf

%% ygrid calc - scaled from xgrid... (not used except for visualising)

dyf      =   zeros(je-jb+1,1);
yh = zeros(je-jb+2,1);
dyf(jb)  =   hy/je;
yh(jb)   =   0;
yh(jb+1) =   dyf(jb);

for j=jb+1:je
    dyf(j)=dyf(j-1);
    yh(j+1)=yh(j)+dyf(j);
end

yf = zeros(je-jb+1,1);

for j=jb:je
    yf(j) = 0.5*(yh(j+1)+yh(j));
end

%% zgrid calc

if r.lzstretch

%     zgrid=r.dzlin/2:r.dzlin:r.nlin;
%     nzn=length(zgrid);
%     
%     z=zeros(ke-nzn,1);
%     z(1)=zgrid(end);
%     dzn=r.dzlin;
%     
%     gf=r.gf; % as r.gf is read only
%     
%     while true
%         dzn=r.dzlin;
%         for i=2:(ke-nzn)
%             dzn=dzn*gf;
%             z(i)=z(i-1)+dzn;
%         end
%         if z(end)+1.5*gf*dzn>hz
%             break
%         else
%             gf=gf+0.000001;
%         end
%         
%         gf
%         
%     end
%     zf=zeros(ke,1);
%     zf(1:nzn)=zgrid;
%     zf(nzn:ke-1)=z;
%     zf(end)=zf(end-1)+dzn*gf;
%     
%     zh(kb) = 0;
%     for k=kb+1:ke
%         zh(k) = 0.5*(zf(k)+zf(k-1));
%     end
%     
%     figure
%     plot(gradient(zf),kb:ke)
%     xlabel('dz [m]')
%     ylabel('zi [-]')

maxh=r.nlin;
dz = r.dzlin;
nk=ke;
dh=r.zsize;
stretchconst=2;
stretch = 'tanh';
source = 2;

makezgridOG
    
else
    
    nit=100000;     % numbe
    
    dzf   =  zeros(ke-kb+1,1);
    zh   =  zeros(ke-kb+2,1);
    %dxf(ib) =     0.025/26;  % first cell length
    dzf(kb) =    hz/ke;  % first cell length  (from 'main_vertical_0toh_48.m)
    zh(kb) = 0;
    zh(kb+1) = dzf(kb); % first grid point
    
    for n=1:nit
        
        % determine x-grid
        for i=kb+1:ke
            dzf(i)=dzf(i-1)*growthz;
            zh(i+1)=zh(i)+dzf(i);
        end
        
        if (zh(end) < z_top-0.00000000000000005 || zh(end) > z_top+0.00000000000000005)
            
            %growthx=growthx*(x_top/(xh(end)))^0.005;
            
        else
            display('converged!')
            display('zf(end)+0.5*dzf(end)= ')
            break
        end
    end
    
    zf = zeros(ke-kb+1,1);
    
    for i=kb:ke
        zf(i) = 0.5*(zh(i+1)+zh(i));
    end
    
    % add the offset (top of the block)
    zf  = zf  + offset;
    zh  = zh  + offset;
    
end

%% Create xgrid.inp.025
xgrid = fopen(['xgrid.inp.' expnr], 'w');
fprintf(xgrid, '%12s\n', '#     x-grid');
fprintf(xgrid, '%12s\n', '#           ');
fprintf(xgrid, '%-20.15f\n', xf);
fclose(xgrid);
disp(['... written xgrid.inp.' expnr])

%% Creat zgrid.inp.025
zgrid = fopen(['zgrid.inp.' expnr], 'w');
fprintf(zgrid, '%12s\n', '#     z-grid');
fprintf(zgrid, '%12s\n', '#           ');
fprintf(zgrid, '%-20.15f\n', zf);
fclose(zgrid);
disp(['... written zgrid.inp.' expnr])

%% Create lscale.inp.025
ls = zeros(length(zf),10);
ls(:,1) = zf;
ls(:,2) = ug;
ls(:,3) = vg;
ls(:,4) = pqx;
ls(:,5) = pqy;
ls(:,6) = wfls;
ls(:,10) = dthlrad;
lscale = fopen(['lscale.inp.' expnr], 'w');
fprintf(lscale, '%-12s\n', '# SDBL flow');
fprintf(lscale, '%-60s\n', '# z uq vq pqx pqy wfls dqtdxls dqtdyls dqtdtls dthlrad');
fprintf(lscale, '%-20.15f %-12.6f %-12.6f %-12.6f %-12.6f %-15.9f %-12.6f %-12.6f %-12.6f %-17.12f\n', ls');
fclose(lscale);
disp(['... written lscale.inp.' expnr])

%% Determine prof.inp
%qt0=0;
pr = zeros(length(zf),6);
pr(:,1) = zf;
pr(:,2) = thl;
pr(:,3) = qt0;
pr(:,4) = u;
pr(:,5) = v;
%pr(:,6) = tke;
prof = fopen(['prof.inp.' expnr], 'w');
fprintf(prof, '%-12s\n', '# SDBL flow');
fprintf(prof, '%-60s\n', '# z thl qt u v tke');
fprintf(prof, '%-20.15f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n', pr');
fclose(prof);
disp(['... written prof.inp.' expnr])

pwd

%% Determine scalar.inp

sc = zeros(length(zf),5);

if r.lchem
    
    % load('/projects/dales-urban/bin/profs')
    
    %    sc(:,4) = linspace(25e-6,25e-6,length(zf));  % linspace(86e-6,86e-6,length(zf)); %Initial O3
    %    sc(1:r.blockheight,4) = linspace(38.4e-6,38.4e-6,r.blockheight); %Initial NO
    %    sc(:,3) = linspace(2e-6,2e-6,length(zf)); %linspace(11.3e-6,11.3e-6,length(zf)); %Initial NO2
    %    sc(1:r.blockheight,3) = linspace(11.5e-6,11.5e-6,r.blockheight); %Initial NO
    %    sc(:,2) = linspace(20e-6,20e-6,length(zf));  %linspace(3.75e-6,3.75e-6,length(zf)); %Initial NO
    %    sc(1:r.blockheight,2) = linspace(75e-6,75e-6,r.blockheight); %Initial NO
    
    % was setting particular profile for a validation case
    %sc(kb:min(ke,length(NOxyprof)),2) = NOxyprof(kb:min(ke,length(NOxyprof))); % did all this for a specific case where I wanted to use a Case A profile but the domain that I was putting it on was bigger than the chem validation... 08/12/17
    %sc(kb:min(ke,length(NOxyprof)),3) = NO2xyprof(kb:min(ke,length(NOxyprof)));
    %sc(kb:min(ke,length(NOxyprof)),4) = O3xyprof(kb:min(ke,length(NOxyprof)));
    
    %if ke>length(NOxyprof)
    %    sc(length(NOxyprof)+1:ke,2) = sc(length(NOxyprof),2);
    %    sc(length(NOxyprof)+1:ke,3) = sc(length(NOxyprof),3);
    %    sc(length(NOxyprof)+1:ke,4) = sc(length(NOxyprof),4);
    %end
    
    sc(1:ke,2) = r.NOb;
   % sc(11:30,2) = linspace(18.54e-6,0,20);
   % sc(31:ke,2) = 0;
    
    sc(1:ke,3) = r.NO2b;
   % sc(11:30,3) = linspace(35.88e-6,0,20);
   % sc(31:ke,3) = 0;
    
    sc(1:ke,4) = r.O3b;
    
    sc(1:ke,5) = r.NOb + r.NO2b;
   % sc(11:30,4) = linspace(50.96e-6,0,20);
   % sc(31:ke,4) = 0;
    
    figure
    plot(sc(:,2))
    hold on
    plot(sc(:,3))
    hold on
    plot(sc(:,4))
    
    
end

sc(:,1) = zf;
sc(:,2) = sv10;
sc(:,3) = sv20;
scalar = fopen(['scalar.inp.' expnr], 'w');
fprintf(scalar, '%-12s\n', '# SDBL flow');
fprintf(scalar, '%-60s\n', '# z sca1 sca2 sca3 sca4');
fprintf(scalar, '%-20.15f %-14.10f %-14.10f %-14.10f %-14.10f\n', sc');
fclose(scalar);
disp(['... written scalar.inp.' expnr])

%% Determine blocks.inp
aspectratio = zh(r.blockheight+1)/(xh(r.canyonwidth+1)-xh(1));
nrows = ie/(r.blockwidth+r.canyonwidth);

%error test

if r.lflat
    
elseif ceil(nrows) ~= floor(nrows)
    l = 0:ie/2;
    ind = rem(ie/2,l)==0; %// logical index that tells if remainder is zero or not
    err = ([l(ind); (ie/2)./l(ind)]);
    disp(['Block system does not fit grid'])
    disp(['sum widths to: ' num2str(err(1,:))])
    disp(['Current width: ' num2str(r.blockwidth + r.canyonwidth)])
    error('Incorrect block system')
else
    disp(['Successful block network'])
    disp(['aspect ratio: ' num2str(aspectratio)])
end

if r.lflat
    
    !    bl(:,:) = [];
    
elseif r.lcastro
    
    %block starting at 1 (got errors in DALES)
    %     nrows = ie/(r.blockwidth*2);
    %     ncolumns = je/(r.blockwidth*2);
    %
    %     bl = zeros(nrows*ncolumns,13);
    %
    %     for n = 1:nrows
    %         for nn = 0:ncolumns-1
    %
    %             bl((n-1)*ncolumns+nn+1,1) = -0.5*r.blockwidth + (2*n-1) * r.blockwidth + 1;
    %             bl((n-1)*ncolumns+nn+1,2) = bl((n-1)*ncolumns+nn+1,1) + r.blockwidth - 1;
    %             bl((n-1)*ncolumns+nn+1,5) = 0;
    %             bl((n-1)*ncolumns+nn+1,6) = ceil(r.blockwidth*(h/ie)/(hz/ke));
    %
    %             if mod(n,2) == 0
    %                 bl((n-1)*ncolumns+nn+1,3) = jb + nn * r.blockwidth*2;
    %                 bl((n-1)*ncolumns+nn+1,4) = bl((n-1)*ncolumns+nn+1,3) + r.blockwidth - 1;
    %             else
    %                 bl((n-1)*ncolumns+nn+1,3) = jb + r.blockwidth + nn * r.blockwidth*2;
    %                 bl((n-1)*ncolumns+nn+1,4) = bl((n-1)*ncolumns+nn+1,3) + r.blockwidth - 1;
    %             end
    %
    %         end
    %     end
    %
    % elseif r.lcube
    %
    %     nrows = ie/(r.blockwidth*2);
    %     ncolumns = je/(r.blockwidth*2);
    %
    %     bl = zeros(nrows*ncolumns,13);
    %
    %     for n = 1:nrows
    %         for nn = 0:ncolumns-1
    %
    %             bl((n-1)*ncolumns+nn+1,1) = -0.5*r.blockwidth + (2*n-1) * r.blockwidth + 1;
    %             bl((n-1)*ncolumns+nn+1,2) = bl((n-1)*ncolumns+nn+1,1) + r.blockwidth - 1;
    %             bl((n-1)*ncolumns+nn+1,5) = 0;
    %             bl((n-1)*ncolumns+nn+1,6) = ceil(r.blockwidth*(h/ie)/(hz/ke));
    %             bl((n-1)*ncolumns+nn+1,3) = jb + nn * r.blockwidth*2;
    %             bl((n-1)*ncolumns+nn+1,4) = bl((n-1)*ncolumns+nn+1,3) + r.blockwidth - 1;
    %
    %         end
    %     end
    
    %overlap of blocks
    nrows = ie/(r.blockwidth*2);
    ncolumns = je/(r.blockwidth*2);
    
    bl = zeros(nrows*ncolumns+nrows/2,13);
    
    bl(:,5) = 0;
    bl(:,6) = r.blockheight-1; %ceil(r.blockwidth*(h/ie)/(hz/ke));
    
    for n = 1:nrows
        for nn = 0:ncolumns
            
            
            if mod(n,2) == 0
                if nn==0
                    bl(length(nonzeros(bl(:,3)))+1,3) = jb;
                    bl(length(nonzeros(bl(:,4)))+1,4) = r.blockwidth/2;
                elseif nn==ncolumns
                    bl(length(nonzeros(bl(:,3)))+1,3) = je-r.blockwidth/2+1;
                    bl(length(nonzeros(bl(:,4)))+1,4) = je;
                else
                    bl(length(nonzeros(bl(:,3)))+1,3) = jb + nn * r.blockwidth*2 - r.blockwidth/2;
                    bl(length(nonzeros(bl(:,4)))+1,4) = bl(length(nonzeros(bl(:,4)))+1,3) + r.blockwidth - 1;
                end
                bl(length(nonzeros(bl(:,1)))+1,1) = -0.5*r.blockwidth + (2*n-1) * r.blockwidth + 1;
                bl(length(nonzeros(bl(:,2)))+1,2) = bl(length(nonzeros(bl(:,2)))+1,1) + r.blockwidth - 1;
            end
        end
        for nn = 0:ncolumns-1
            if mod(n,2) ~= 0
                bl(length(nonzeros(bl(:,3)))+1,3) = jb + r.blockwidth + nn * r.blockwidth*2 - r.blockwidth/2;
                bl(length(nonzeros(bl(:,4)))+1,4) = bl(length(nonzeros(bl(:,4)))+1,3) + r.blockwidth - 1;
                bl(length(nonzeros(bl(:,1)))+1,1) = -0.5*r.blockwidth + (2*n-1) * r.blockwidth + 1;
                bl(length(nonzeros(bl(:,2)))+1,2) = bl(length(nonzeros(bl(:,2)))+1,1) + r.blockwidth - 1;
            end
            
        end
    end
    
elseif r.lcube
    
    nrows = ie/(r.blockwidth*2);
    ncolumns = je/(r.blockwidth*2);
    
    bl = zeros(nrows*ncolumns,11); %13
    
    for n = 1:nrows
        for nn = 0:ncolumns-1
            
%             bl((n-1)*ncolumns+nn+1,1) = -0.5*r.blockwidth + (2*n-1) * r.blockwidth + 1;
%             bl((n-1)*ncolumns+nn+1,2) = bl((n-1)*ncolumns+nn+1,1) + r.blockwidth - 1;
%             bl((n-1)*ncolumns+nn+1,5) = 0;
%             bl((n-1)*ncolumns+nn+1,6) = ceil(r.blockwidth*(h/ie)/(hz/ke));
%             bl((n-1)*ncolumns+nn+1,3) = jb + nn * r.blockwidth*2;
%             bl((n-1)*ncolumns+nn+1,4) = bl((n-1)*ncolumns+nn+1,3) + r.blockwidth - 1;
            bl((n-1)*ncolumns+nn+1,1) = -0.5*r.blockwidth + (2*n-1) * r.blockwidth + 1;
            bl((n-1)*ncolumns+nn+1,2) = bl((n-1)*ncolumns+nn+1,1) + r.blockwidth - 1;
            bl((n-1)*ncolumns+nn+1,5) = 0;
            bl((n-1)*ncolumns+nn+1,6) = ceil(r.blockwidth*(h/ie)/(hz/ke));
            bl((n-1)*ncolumns+nn+1,3) = jb + r.blockwidth/2 + nn * r.blockwidth*2;
            bl((n-1)*ncolumns+nn+1,4) = bl((n-1)*ncolumns+nn+1,3) + r.blockwidth - 1;
            
        end
    end
    
elseif r.lblocks
    %spanwise blocks
    
    bl = zeros(nrows, 13);
    bl(1:nrows,1) = (r.canyonwidth/2+1:r.canyonwidth+r.blockwidth:ie-r.canyonwidth/2)';
    bl(1:nrows,2) = bl(1:nrows,1) + r.blockwidth - 1;
    bl(:,3) = jb;
    bl(:,4) = je;
    bl(:,5) = 0;
    bl(1:nrows,6) = r.blockheight-1;
    %bl(nrows+1:end,6) = 1;
    nrows1 = 0;
    
    if r.ltrees %niru's code, may need adjusting for no blocks in canyons...
        if nt2
            trees= zeros(nrows*2, 6);
            for i = 1:nrows
                trees(i,1) =  bl(i,1) - tree_dx - tree_w;
                trees(i,2) =  bl(i,1) - tree_dx;
            end
            for i = 1:nrows
                trees(nrows+i,1) =  bl(i,2) + tree_dx;
                trees(nrows+i,2) =  bl(i,2) + tree_dx + tree_w;
            end
            
            trees(:,3) = 1;
            trees(:,4) = je;
            
            trees(:,5) = tree_dz;
            trees(:,6) = tree_dz + tree_h;
            
        elseif nt1 && md
            trees= zeros(nrows+1, 6);
            tree_dx=0;
            for i = 1:size(trees,1)
                if i==1 || i==size(trees,1)
                    xh(trees(i,1)) =  0.5*(bl(nrows+i,2)+bl(nrows+i,1)) - tree_dx - 0.5*tree_w/2;
                    xh(trees(i,2)+1) =  0.5*(bl(nrows+i,2)+bl(nrows+i,1)) - tree_dx + 0.5*tree_w/2;
                    
                else
                    xh(trees(i,1)) =  0.5*(bl(nrows+i,2)+bl(nrows+i,1)) - tree_dx - tree_w/2;
                    xh(trees(i,2)+1) =  0.5*(bl(nrows+i,2)+bl(nrows+i,1)) - tree_dx + tree_w/2;
                    
                end
            end
            trees(:,3) = 1;
            trees(:,4) = je;
            
            trees(:,5) = tree_dz - tree_h/2;
            trees(:,6) = tree_dz + tree_h/2;
        elseif nt1 && lw
            
            trees= zeros(nrows, 6);
            
            for i = 1:nrows
                xh(trees(i,1)) =  bl(nrows+1+i,1) + tree_dx - tree_w/2;
                xh(trees(i,2)+1) =  bl(nrows+1+i,1) + tree_dx + tree_w/2;
            end
            
            trees(:,3) = 1;
            trees(:,4) = je;
            
            trees(:,5) = tree_dz - tree_h/2;
            trees(:,6) = tree_dz + tree_h/2;
            trees(:,3) = 1;
            trees(:,4) = je;
            
            trees(:,5) = tree_dz - tree_h/2;
            trees(:,6) = tree_dz + tree_h/2;
            
        elseif nt1 && ww
            trees= zeros(nrows, 6);
            for i = 1:nrows
                xh(trees(i,1)) =  bl(nrows+i,2) - tree_dx - tree_w/2;
                xh(trees(i,2)+1) =  bl(nrows+i,2) - tree_dx + tree_w/2;
            end
            
            trees(:,3) = 1;
            trees(:,4) = je;
            
            trees(:,5) = tree_dz - tree_h/2;
            trees(:,6) = tree_dz + tree_h/2;
            
        else
            error('Need to specify tree setup (nt1, nt2 etc.)')
        end
    end
    
    if r.lpurif   % also needs updating (may be om the whole a little unnecessary (can hard code in if npurif = small))
        
        purifs = zeros(nrows*2*npurif, 7);
        
        for i = 1:nrows
            for j = 1:npurif
                purifs((i-1)*npurif+j,1) = bl(nrows+i,2) - purif_dx - purif_w;
                purifs((i-1)*npurif+j,2) = bl(nrows+i,2) - purif_dx;
                if j==1
                    purifs((i-1)*npurif+j,3) = ((je/npurif)/2);
                    purifs((i-1)*npurif+j,4) = purifs((i-1)*npurif+j,3) + purif_dy;
                else
                    purifs((i-1)*npurif+j,3) = purifs((i-1)*npurif+j-1,3) + purif_dy + purif_sp;
                    purifs((i-1)*npurif+j,4) = purifs((i-1)*npurif+j-1,4) + purif_dy + purif_sp;
                end
            end
        end
        for i = 1:nrows
            for j = 1:npurif
                purifs(nrows*npurif + (i-1)*npurif+j,1) =  bl(nrows+1+i,1) + purif_dx;
                purifs(nrows*npurif + (i-1)*npurif+j,2) =  bl(nrows+1+i,1) + purif_dx + purif_w;
                if j==1
                    purifs(nrows*npurif + (i-1)*npurif+j,3) = ((je/npurif)/2);
                    purifs(nrows*npurif + (i-1)*npurif+j,4) = purifs((i-1)*npurif+j,3) + purif_dy;
                else
                    purifs(nrows*npurif + (i-1)*npurif+j,3) = purifs((i-1)*npurif+j-1,3) + purif_dy + purif_sp;
                    purifs(nrows*npurif + (i-1)*npurif+j,4) = purifs((i-1)*npurif+j-1,4) + purif_dy + purif_sp;
                end
            end
        end
        
        purifs(:,5) = purif_dz;
        purifs(:,6) = purif_dz + purif_h;
        purifs(:,7) = purif_i;
        
    end
    
end

% bl other variables
bl(:,7) = r.z0horiz;
bl(:,8) = r.z0hhoriz;
bl(:,9) = r.Thoriz;
bl(:,10) = r.Twest;
bl(:,11) = r.Teast;
bl(:,12) = r.Tnorth;
bl(:,13) = r.Tsouth;

blocks = fopen( ['blocks.inp.' expnr], 'w');
fprintf(blocks, '%-12s\n', '# Fence location');
fprintf(blocks, '%-100s\n', '#  il  iu  jl  ju  kl  ku  z0horiz[m]  z0hhoriz[m]  Thoriz[K]  Twest  Teast  Tnorth  Tsouth');
fprintf(blocks, '%-3.0f %-3.0f %-3.0f %-3.0f %-3.0f %-3.0f %-8.6f %8.6f %-4.2f %-4.2f %-4.2f %-4.2f %-4.2f\n', bl');
fclose(blocks);
disp(['... written blocks.inp.' expnr])
disp(['nblocks : ' num2str(size(bl,1))])

if r.ltrees
    tree_write = fopen( ['trees.inp.' expnr], 'w');
    fprintf(tree_write, '%-12s\n', '# Tree location');
    fprintf(tree_write, '%-60s\n', '#  il  iu  jl  ju  kl  ku');
    fprintf(tree_write, '%-3.0f %-3.0f %-3.0f %-3.0f %-3.0f %-3.0f\n', trees');
    fclose(tree_write);
    disp(['... written trees.inp.' expnr])
end

if r.lpurif
    purif_write = fopen( ['purifs.inp.' expnr], 'w');
    fprintf(purif_write, '%-12s\n', '# Purifier location');
    fprintf(purif_write, '%-60s\n', '#  il  iu  jl  ju  kl  ku  ipu');
    fprintf(purif_write, '%-3.0f %-3.0f %-3.0f %-3.0f %-3.0f %-3.0f %-3.0f\n', purifs');
    fclose(purif_write);
    disp(['... written purifs.inp.' expnr])
end

%% Plots

figure
view(52,23)

if (r.lcastro || r.lcube || r.lblocks)
    for i = 1:size(bl,1)
        
        patch([xh(bl(i,1)) xh(bl(i,2)+1) xh(bl(i,2)+1)  xh(bl(i,1))], [yh(bl(i,3))  yh(bl(i,3)) yh(bl(i,4)+1) yh(bl(i,4)+1)], [zh(bl(i,6)+1)  zh(bl(i,6)+1) zh(bl(i,6)+1) zh(bl(i,6)+1)], [245 245 245] ./ 255) % [245 245 245] ./ 255
        patch([xh(bl(i,1)) xh(bl(i,1)) xh(bl(i,2)+1) xh(bl(i,2)+1) ], [yh(bl(i,3))  yh(bl(i,3)) yh(bl(i,3)) yh(bl(i,3))], [bl(i,5)  zh(bl(i,6)+1) zh(bl(i,6)+1) bl(i,5)], [245 245 245] ./ 255)
        patch([xh(bl(i,1)) xh(bl(i,1)) xh(bl(i,2)+1) xh(bl(i,2)+1) ], [yh(bl(i,4)+1)  yh(bl(i,4)+1) yh(bl(i,4)+1) yh(bl(i,4)+1)], [bl(i,5)  zh(bl(i,6)+1) zh(bl(i,6)+1) bl(i,5)], [245 245 245] ./ 255)
        patch([xh(bl(i,1)) xh(bl(i,1)) xh(bl(i,1)) xh(bl(i,1)) ], [yh(bl(i,4)+1)  yh(bl(i,4)+1) yh(bl(i,3)) yh(bl(i,3))], [bl(i,5)  zh(bl(i,6)+1) zh(bl(i,6)+1) bl(i,5)], [245 245 245] ./ 255)
        patch([xh(bl(i,2)+1) xh(bl(i,2)+1) xh(bl(i,2)+1) xh(bl(i,2)+1) ], [yh(bl(i,3))  yh(bl(i,3)) yh(bl(i,4)+1) yh(bl(i,4)+1)], [bl(i,5)  zh(bl(i,6)+1) zh(bl(i,6)+1) bl(i,5)], [245 245 245] ./ 255)
        % patch([xh(1) xh(end) xh(end)  xh(1)], [yh(1)  yh(1) yh(end) yh(end)], [zh(1)  zh(1) zh(1) zh(1)], [245 245 245] ./ 255)
    end
elseif r.lflat
    
else
    for i = 1:size(bl,1)
        
        patch([xh(bl(i,1)) xh(bl(i,2)+1) xh(bl(i,2)+1)  xh(bl(i,1))], [yh(bl(i,3))  yh(bl(i,3)) yh(bl(i,4)+1) yh(bl(i,4)+1)], [zh(bl(i,6)+1)  zh(bl(i,6)+1) zh(bl(i,6)+1) zh(bl(i,6)+1)], 'w')
        patch([xh(bl(i,1)) xh(bl(i,1)) xh(bl(i,2)+1) xh(bl(i,2)+1) ], [yh(bl(i,3))  yh(bl(i,3)) yh(bl(i,3)) yh(bl(i,3))], [bl(i,5)  zh(bl(i,6)+1) zh(bl(i,6)+1) bl(i,5)], 'w')
        patch([xh(bl(i,1)) xh(bl(i,1)) xh(bl(i,2)+1) xh(bl(i,2)+1) ], [yh(bl(i,4)+1)  yh(bl(i,4)+1) yh(bl(i,4)+1) yh(bl(i,4)+1)], [bl(i,5)  zh(bl(i,6)+1) zh(bl(i,6)+1) bl(i,5)], 'w')
        patch([xh(bl(i,1)) xh(bl(i,1)) xh(bl(i,1)) xh(bl(i,1)) ], [yh(bl(i,4)+1)  yh(bl(i,4)+1) yh(bl(i,3)) yh(bl(i,3))], [bl(i,5)  zh(bl(i,6)+1) zh(bl(i,6)+1) bl(i,5)], 'w')
        patch([xh(bl(i,2)+1) xh(bl(i,2)+1) xh(bl(i,2)+1) xh(bl(i,2)+1) ], [yh(bl(i,3))  yh(bl(i,3)) yh(bl(i,4)+1) yh(bl(i,4)+1)], [bl(i,5)  zh(bl(i,6)+1) zh(bl(i,6)+1) bl(i,5)], 'w')
        
    end
end

if r.ltrees  % put interms of xh, zh etc.
    for i = 1:size(trees,1)
        
        patch([xh(trees(i,1)) xh(trees(i,2)+1) xh(trees(i,2)+1)  xh(trees(i,1))], [yh(trees(i,3))  yh(trees(i,3)) yh(trees(i,4)+1) yh(trees(i,4)+1)], [zh(trees(i,6)+1)  zh(trees(i,6)+1) zh(trees(i,6)+1) zh(trees(i,6)+1)], [128 229 78] ./ 255)
        patch([xh(trees(i,1)) xh(trees(i,1)) xh(trees(i,2)+1) xh(trees(i,2)+1) ], [yh(trees(i,3))  yh(trees(i,3)) yh(trees(i,3)) yh(trees(i,3))], [zh(trees(i,5))  zh(trees(i,6)+1) zh(trees(i,6)+1) zh(trees(i,5))], [128 229 78] ./ 255)
        patch([xh(trees(i,1)) xh(trees(i,1)) xh(trees(i,2)+1) xh(trees(i,2)+1) ], [yh(trees(i,4)+1)  yh(trees(i,4)+1) yh(trees(i,4)+1) yh(trees(i,4)+1)], [zh(trees(i,5))  zh(trees(i,6)+1) zh(trees(i,6)+1) zh(trees(i,5))],[128 229 78] ./ 255 )
        patch([xh(trees(i,1)) xh(trees(i,1)) xh(trees(i,1)) xh(trees(i,1)) ], [yh(trees(i,4)+1)  yh(trees(i,4)+1) yh(trees(i,3)) yh(trees(i,3))], [zh(trees(i,5))  zh(trees(i,6)+1) zh(trees(i,6)+1) zh(trees(i,5))],[128 229 78] ./ 255 )
        patch([xh(trees(i,2)+1) xh(trees(i,2)+1) xh(trees(i,2)+1) xh(trees(i,2)+1) ], [yh(trees(i,3))  yh(trees(i,3)) yh(trees(i,4)+1) yh(trees(i,4)+1)], [zh(trees(i,5))  zh(trees(i,6)+1) zh(trees(i,6)+1) zh(trees(i,5))], [128 229 78] ./ 255)
        
    end
end

if r.lpurif
    for i = 1:size(purifs,1)
        
        patch([xh(purifs(i,1)) xh(purifs(i,2)+1) xh(purifs(i,2)+1)  xh(purifs(i,1))], [yh(purifs(i,3))  yh(purifs(i,3)) yh(purifs(i,4)+1) yh(purifs(i,4)+1)], [zh(purifs(i,6)+1)  zh(purifs(i,6)+1) zh(purifs(i,6)+1) zh(purifs(i,6)+1)], [245 245 245] ./ 255)
        patch([xh(purifs(i,1)) xh(purifs(i,1)) xh(purifs(i,2)+1) xh(purifs(i,2)+1) ], [yh(purifs(i,3))  yh(purifs(i,3)) yh(purifs(i,3)) yh(purifs(i,3))], [zh(purifs(i,5))  zh(purifs(i,6)+1) zh(purifs(i,6)+1) zh(purifs(i,5))], [245 245 245] ./ 255)
        patch([xh(purifs(i,1)) xh(purifs(i,1)) xh(purifs(i,2)+1) xh(purifs(i,2)+1) ], [yh(purifs(i,4)+1)  yh(purifs(i,4)+1) yh(purifs(i,4)+1) yh(purifs(i,4)+1)], [zh(purifs(i,5))  zh(purifs(i,6)+1) zh(purifs(i,6)+1) zh(purifs(i,5))], [245 245 245] ./ 255)
        patch([xh(purifs(i,1)) xh(purifs(i,1)) xh(purifs(i,1)) xh(purifs(i,1)) ], [yh(purifs(i,4)+1)  yh(purifs(i,4)+1) yh(purifs(i,3)) yh(purifs(i,3))], [zh(purifs(i,5))  zh(purifs(i,6)+1) zh(purifs(i,6)+1) zh(purifs(i,5))], [245 245 245] ./ 255)
        patch([xh(purifs(i,2)+1) xh(purifs(i,2)+1) xh(purifs(i,2)+1) xh(purifs(i,2)+1) ], [yh(purifs(i,3))  yh(purifs(i,3)) yh(purifs(i,4)+1) yh(purifs(i,4)+1)], [zh(purifs(i,5))  zh(purifs(i,6)+1) zh(purifs(i,6)+1) zh(purifs(i,5))], [245 245 245] ./ 255)
        
    end
end

zlim([0 zh(end)]); %/(r.blockheight-1))
xlim([0 xh(end)]); %/(r.blockheight-1))
ylim([0 yh(end)]); %/(r.blockheight-1))

set(gca,'ticklabelinterpreter','latex')
xlabel('x [m]','interpreter','latex')
ylabel('y [m]','interpreter','latex')
zlabel('z [m]','interpreter','latex')
set(gca,'BoxStyle','full','Box','on')
ax = gca;

%    ax.YTick = [0,je/2,je]; %/(r.blockheight-1);%
%    ax.ZTick = [0,ie/2,ie]; %(r.blockheight-1);
%    ax.XTick = [0,ie/4,ie/2,ie*0.75,ie]; %/(r.blockheight-1);

%    ax.YTickLabels = [0,(h/2)*(je/ie),h*(je/ie)]; %/(r.blockheight-1);
%    ax.ZTickLabels = [0,hz/2,hz]; %(r.blockheight-1);
%    ax.XTickLabels = [0,h/4,h/2,h*0.75,h]; %/(r.blockheight-1);

%    set(ax,'FontSize',14,'YTick',[0,je/2,je],'ZTick',[0,ie/2,ie],XTick,[0,ie/4,ie/2,ie*0.75,ie],'YTickLabels',[0,(h/2)*(je/ie),h*(je/ie)],'ZTickLabels',[0,hz/2,hz],'XTickLabels',[0,h/4,h/2,h*0.75,h])
%    set(ax,'FontSize',14)
daspect([1 1 1])
%pbaspect([1 1 3])

%addpath(genpath([DA_TOPDIR '/bin/b2f/']));

%main

disp('Check plots and press a key to continue...')
pause

disp([expnr ' .inp. files succesfully created'])
