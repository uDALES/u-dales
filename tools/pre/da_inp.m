%% da_inp

expnr = '009';
ncpus = 2;
LIDAR = 0;

DA_EXPDIR = getenv('DA_EXPDIR');
DA_PREDIR = getenv('DA_PREDIR');
addpath([DA_PREDIR '/']);
exppath = [DA_EXPDIR '/'];
cd([DA_EXPDIR '/' expnr])

copyfile([DA_PREDIR, '/default/walltypes.inp.xxx'], ['walltypes.inp.', expnr]);

r = da_pp(expnr, exppath);
da_pp.addvar(r, 'walltypes', dlmread(['walltypes.inp.', expnr],'',3,0));
da_pp.set_defaults(r, ncpus);


da_pp.generate_xygrid(r);
da_pp.write_xgrid(r)
disp(['Written xgrid.inp.', r.expnr])
da_pp.generate_zgrid(r);
da_pp.write_zgrid(r);
disp(['Written zgrid.', r.expnr])
da_pp.generate_lscale(r)
da_pp.write_lscale(r)
disp(['Written lscal.inp.', r.expnr])
da_pp.generate_prof(r);
da_pp.write_prof(r);
disp(['Written prof.inp.', r.expnr])
%da_pp.plot_profiles(r);

da_pp.generate_scalar(r);
da_pp.write_scalar(r);
disp(['Written scalar.inp.', r.expnr])

%addpath(genpath(DA_PREDIR));

if ~LIDAR
    da_pp.generate_bl_from_namoptions(r)
    da_pp.generate_topo_from_bl(r)
else
    sourcename;
    % resolution of image [m/pixel]
    dxinp; dyinp; dzinp;
    % center of area of interest in original image [pixel]
    centeri; centerj;
    % magimum height in image [m]
    maxh;
    %padding. A padding of 0 makes only sense for idealised cases. There should be no building at domain edge
    pad;
    %objects smaller than this will be deleted
    smallarea;
    
    da_pp.generate_topo_from_LIDAR(r, sourcename, dxinp, dyinp, centeri, centerj, maxh, pad, smallarea)
end

da_pp.makeblocks(r)
da_pp.generate_facets(r);

da_pp.write_blocks(r)
disp(['Written blocks.inp.', r.expnr])
da_pp.write_facets(r)
disp(['Written facets.inp.', r.expnr])
da_pp.plot_blocks(r)
da_pp.plot_facets(r)


if r.ltrees
    da_pp.generate_trees(r, true);
    disp(['Written trees.inp.', r.expnr])
end

if r.lpurif
    da_pp.generate_purifs(r, true);
    disp(['Written purifs.inp.', r.expnr])
end

if r.lEB
    da_pp.vsolc(r)
    da_pp.plot_shading(r)
    da_pp.vfc(r)
    da_pp.write_svf(r)
    disp(['Written svf.inp.', r.expnr])
    da_pp.write_vf(r)
    disp(['Written vf.nc.inp.', r.expnr])
    da_pp.write_facetarea(r)
    disp(['Written facetarea.inp.', r.expnr])
    da_pp.rayit(r)
    da_pp.write_netsw(r)
    disp(['Written netsw.inp.', r.expnr])
    da_pp.generate_Tfacinit(r)
else
    da_pp.addvar(r, 'Tfacinit', ones(r.nfcts,1) * 288);
end

da_pp.write_Tfacinit(r)
disp(['Written Tfacinit.inp.', r.expnr])

disp('Check plots and press a key to continue...')
pause
