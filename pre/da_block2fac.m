%% File to take old exp folder format and create facets etc.

close all
clear all

% expnr
expnr = '001';

saveasnew     = true     ; %saves all date in a new folder with current datetime

% read env
DA_EXPDIR = getenv('DA_EXPDIR');
DA_TOPDIR = getenv('DA_TOPDIR');
DA_PREDIR = getenv('DA_PREDIR');

% read namoptions
addpath([DA_PREDIR '/']);
exppath=[DA_EXPDIR '/'];
r=da_pp(expnr,exppath);

% recreate matlab variables needed for Ivo's preproc
dx = r.xsize/r.imax;
dy = r.ysize/r.jtot;
dz = r.zsize/r.kmax;
dh = r.zsize;
ni = r.imax; nj = r.jtot; nk = r.kmax;

solaz = 135; % azimuth angle
Z = 28.4066; % zenith angle
delta = 0.01; % small adjustment for position of rad corners
centerweight = 12/32;
cornerweight = (1-centerweight)/4;
I = 184.8775; % Direct solar irradiation [W/m2]
Dsk = 418.8041; % Diffuse incoming radiation [W/m2]
pad=0;

source = 1;
stretch = 'no';
lwritefiles = 1;
ltestplot = 1;
lhqplot = 1;

lradiation = 1;

%maximum size of floors and bounding walls (cells in each dimension)
if lradiation
    maxsize = 10; %ADD TO NAMOPTIONS EVENTUALY
else
    maxsize = inf;
end

if source == 1
    
    % read blocks
    blocks = dlmread([DA_EXPDIR '/' expnr '/blocks.inp.' expnr],'',2,0);
    
elseif source == 2
    
    dx = 2;
    dy = 2;
    dz = 2;
    dh = 200;
    ni = 400; nj = 400; nk = 100;
    
    % DAPPLE LIDAR
    %sourcename = '/projects/Validation/DAPPLE'; % 
    %dxinp=1;	dyinp=1;	dzinp=0.5; % 1  
    %centeri=850;   centerj=650;
    %maxh= 72; %
    %pad=5;
    
    %sourcename = '/projects/Validation/wpd.png';
    sourcename = '/projects/Validation/DAP_bld_rot.mat';
    dxinp=1;	dyinp=1;
    centeri = 550;%1500*0.2;
    centerj = 700;%1500*0.2;
    maxh = 64;
    
%    sourcename = '/projects/NY/NY_3725.png';
%    dxinp=1;	dyinp=1;	dzinp=1;
%    centeri=2000;   centerj=1365;
%    maxh = 1239/3.28084;

    % resolution of image [m/pixel]
    
%    dxinp=1;	dyinp=1;	dzinp=1;
    
    % center of area of interest in original image [pixel]
    

    % magimum height in image [m]
    
    %padding. A padding of 0 makes only sense for idealised cases. There should be no building at domain edge
    pad=5;
    
    %objects smaller than this will be deleted
    
    smallarea=round(400/(dx*dy));
    Astrel = 20; % size of structural element for morphological processing [m^2] (10*dx...)
    %smallarea=round(200/(dx*dy));
    
end

%% initialising

%run check
%checkprep
%
parentdir = '/projects/dales-u/pre-post/pre/';
cd(parentdir)

addpath(genpath(parentdir))  %add all subdirectories to path

if saveasnew
    starttime=datetime;
    outputdir = [parentdir 'output/' datestr(starttime,'yyyy-mm-dd') '-' datestr(starttime,'HH-MM')];
    inputdir = [parentdir 'input/' datestr(starttime,'yyyy-mm-dd') '-' datestr(starttime,'HH-MM')]; 
    tempdir = [parentdir 'temp/' datestr(starttime,'yyyy-mm-dd') '-' datestr(starttime,'HH-MM')];
    mkdir(outputdir)
    mkdir(tempdir)
    clear starttime
else
    outputdir = [parentdir 'output/'];
    tempdir = [parentdir 'temp/'];
end

%outputdir = '/projects/dales-u/exp/010/';

addpath(outputdir)

copyfile('default/walltypes.inp.xxx',[outputdir '/walltypes.inp.' num2str(expnr)]) %copy default walltypes


%% make grids

makexygrid
makezgrid

if source == 1
    
    topomask = zeros(nj,ni);
    topo = zeros(nj,ni);
    
    if isnan(blocks)  
    else
        for n = 1:size(blocks,1)
            topo(blocks(n,3):blocks(n,4),blocks(n,1):blocks(n,2)) = zh(blocks(n,6)+1);
            topomask(blocks(n,3):blocks(n,4),blocks(n,1):blocks(n,2)) = 1;
        end
    end
    figure
    pcolor(topo)
    shading flat
    
elseif source == 2
    
    topo2blocks
    
end

%%

if ~lradiation
    
    makeblocks
    block2fac
    nwallfcts = 0;
    createfloors
    
else
    
    makeblocks
    block2fac
    addboundingwalls
    createfloors
    vsolcalc
    vfc
    rayit
    
end

%%

% %% Manually make facetarea
% 
% fname = [outputdir '/facetarea.inp.' num2str(expnr)];
% fileID = fopen(fname,'W');
% fprintf(fileID, '# area of facets\n');
% fclose(fileID);
% dlmwrite(fname,A,'-append','delimiter',' ','precision','%4f')

%% Manually make Tfacinit

if ~lradiation

    Tnew = ones(nfcts,1)*288;

    fname = [outputdir '/Tfacinit.inp.' num2str(expnr)];
    fileID = fopen(fname,'W');
    fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium\n');
    fclose(fileID);
    dlmwrite(fname,Tnew(:,1),'-append','delimiter',' ','precision','%4f')

end