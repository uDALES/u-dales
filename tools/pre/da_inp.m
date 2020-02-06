%% da_inp
% This script is run by the bash script da_inp. It used to generate the
% necessary input files for uDALES.
expnr = '100';
ncpus = 2;

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
da_pp.generate_scalar(r);
da_pp.write_scalar(r);
disp(['Written scalar.inp.', r.expnr])

if ~r.llidar
    if r.lflat
        da_pp.addvar(r, 'bl', [])
    else
        if (r.lcastro || r.lcube || r.lcanyons)
            disp('Generating blocks from namoptions')
            da_pp.generate_bl_from_namoptions(r)
        elseif r.lblocksfile
            disp('Generating blocks from file')
            da_pp.generate_bl_from_file(r) 
        end
        da_pp.generate_topo_from_bl(r)
    end
else
    disp('Generating blocks from LIDAR data')
    da_pp.generate_topo_from_LIDAR(r)
end

if ~r.lflat
    da_pp.makeblocks(r)
    da_pp.block2fac(r)
    da_pp.addvar(r, 'nboundingwallfacets', 0)
    if r.lEB
        da_pp.addboundingwalls(r)
    end
else
    da_pp.addvar(r, 'blocks', []);
    da_pp.addvar(r, 'buildings', []);
    da_pp.addvar(r, 'facets', []);
    da_pp.addvar(r, 'nblockfcts', 0);
    da_pp.addvar(r, 'nboundingwallfacets', 0) % lflat not currently compatible with lEB
    da_pp.addvar(r, 'boundingwallfacets', [])
end

da_pp.createfloors(r);
da_pp.write_blocks(r)
disp(['Written blocks.inp.', r.expnr])
da_pp.write_facets(r)
disp(['Written facets.inp.', r.expnr])

if r.lEB
    da_pp.vsolc(r)
    disp('Done vsolc')
    da_pp.vfc(r)
    disp('Done vfc')
    da_pp.write_svf(r)
    disp(['Written svf.inp.', r.expnr])
    da_pp.write_vf(r)
    disp(['Written vf.nc.inp.', r.expnr])
    da_pp.write_facetarea(r)
    disp(['Written facetarea.inp.', r.expnr])
    da_pp.rayit(r)
    disp('Done rayit')
    da_pp.write_netsw(r)
    disp(['Written netsw.inp.', r.expnr])
end

da_pp.generate_Tfacinit(r, r.lEB)
da_pp.write_Tfacinit(r)
disp(['Written Tfacinit.inp.', r.expnr])
