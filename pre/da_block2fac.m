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

% read blocks
blocks = dlmread([DA_EXPDIR '/' expnr '/blocks.inp.' expnr],'',2,0);

% recreate matlab variables needed for Ivo's preproc
dx = r.xsize/r.imax;
dy = r.ysize/r.jtot;
dz = r.zsize/r.kmax;
dh = r.zsize;
ni = r.imax; nj = r.jtot; nk = r.kmax;

solaz = 30; % azimuth angle
Z = 45; % zenith angle
delta = 0.01; % small adjustment for position of rad corners
centerweight = 0.5;
cornerweight = (1-centerweight)/4;
I = 800; % Direct solar irradiation [W/m2]
Dsk = 300; % Diffuse incoming radiation [W/m2]
pad=0;

source = 1;
stretch = 'no';
lwritefiles = 1;
ltestplot = 1;
lhqplot = 0;
lradiation = 0;

%maximum size of floors and bounding walls (cells in each dimension)
if lradiation
    maxsize = 10; %ADD TO NAMOPTIONS EVENTUALY
else
    maxsize = inf; 
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
    outputdir = [parentdir '/output/' datestr(starttime,'yyyy-mm-dd') '-' datestr(starttime,'HH-MM')];
    inputdir = [parentdir '/input/' datestr(starttime,'yyyy-mm-dd') '-' datestr(starttime,'HH-MM')]; 
    tempdir = [parentdir '/temp/' datestr(starttime,'yyyy-mm-dd') '-' datestr(starttime,'HH-MM')];
    mkdir(outputdir)
    mkdir(tempdir)
    clear starttime
else
    outputdir = [parentdir '/output/'];
    tempdir = [parentdir '/temp/'];
end

copyfile('default/walltypes.inp.xxx',[outputdir '/walltypes.inp.' num2str(expnr)]) %copy default walltypes


%% make grids

makexygrid
makezgrid

%% make matrix with block heights

topomask = zeros(nj,ni);
topo = zeros(nj,ni);

for n = 1:size(blocks,1)
   
    topo(blocks(n,1):blocks(n,2),blocks(n,3):blocks(n,4)) = zh(blocks(n,6)+1);
    topomask(blocks(n,1):blocks(n,2),blocks(n,3):blocks(n,4)) = 1;
    
end

%%

if ~lradiation
    
    makefacets
    block2fac
    nwallfcts = 0;
    createfloors
    
else
    
    makefacets
    block2fac
    addboundingwalls
    createfloors
    addboundingwalls
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

% %% Manually make Tfacinit
% 
% Tnew = ones(nfcts,1)*288;
% 
% fname = [outputdir '/Tfacinit.inp.' num2str(expnr)];
% fileID = fopen(fname,'W');
% fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium\n');
% fclose(fileID);
% dlmwrite(fname,Tnew(:,1),'-append','delimiter',' ','precision','%4f')
