%% write_inputs

% uDALES (https://github.com/uDALES/u-dales).

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Copyright (C) 2016-2021 the uDALES Team.

% This script is run by the bash script da_inp.sh.
% It used to generate the necessary input files for uDALES.

expnr = '401';
ncpus = 2;

DA_EXPDIR = getenv('DA_EXPDIR');
DA_TOOLSDIR = getenv('DA_TOOLSDIR');
addpath([DA_TOOLSDIR '/']);
exppath = [DA_EXPDIR '/'];
cd([DA_EXPDIR '/' expnr])

r = preprocessing(expnr, exppath);
preprocessing.set_defaults(r, ncpus);

preprocessing.generate_xygrid(r);
preprocessing.write_xgrid(r)
disp(['Written xgrid.inp.', r.expnr])
preprocessing.generate_zgrid(r);
preprocessing.write_zgrid(r);
disp(['Written zgrid.', r.expnr])
preprocessing.generate_lscale(r)
preprocessing.write_lscale(r)
disp(['Written lscal.inp.', r.expnr])
preprocessing.generate_prof(r);
preprocessing.write_prof(r);
disp(['Written prof.inp.', r.expnr])

if r.nsv>0
    preprocessing.generate_scalar(r);
    preprocessing.write_scalar(r);
    disp(['Written scalar.inp.', r.expnr])
end

if ~r.lflat
    if ~r.lfloors
        if r.llidar
            disp('Generating blocks from LIDAR data')
            preprocessing.generate_topo_from_LIDAR(r)      
        elseif r.ltxtblocks
            disp('Generating blocks from text')
            preprocessing.generate_topo_from_txt(r)
        else
            if (r.lstaggered || r.lstaggeredv || r.lcube || r.lcanyons)
                disp('Generating blocks from namoptions')
                preprocessing.generate_bl_from_namoptions(r)
            elseif r.lblocksfile
                disp('Generating blocks from file')
                preprocessing.generate_bl_from_file(r) 
            end
            preprocessing.generate_topo_from_bl(r)        
        end
    
        preprocessing.makeblocks(r)
        preprocessing.block2fac(r)
    else
        preprocessing.addvar(r, 'blocks', [])
        preprocessing.addvar(r, 'facets', [])
        preprocessing.addvar(r, 'nblockfcts', 0)
    end
    
    if r.lEB
        preprocessing.addboundingwalls(r)
    else
        preprocessing.addvar(r, 'nboundingwallfacets', 0)
    end
    
    preprocessing.createfloors(r);
    preprocessing.write_blocks(r)
    disp(['Written blocks.inp.', r.expnr])
    preprocessing.write_facets(r)
    disp(['Written facets.inp.', r.expnr])
    
    if isfile(['walltypes.inp.', expnr])
        r.walltypes = dlmread(['walltypes.inp.', expnr],'',3,0);
    else
        preprocessing.write_walltypes(r)
        disp(['Written walltypes.inp.', r.expnr])
    end
    
    if r.lEB
        preprocessing.vsolc(r)
        disp('Done vsolc')
        preprocessing.vfc(r)
        disp('Done vfc')
        preprocessing.write_svf(r)
        disp(['Written svf.inp.', r.expnr])
        preprocessing.write_vf(r)
        disp(['Written vf.nc.inp.', r.expnr])
        preprocessing.write_facetarea(r)
        disp(['Written facetarea.inp.', r.expnr])
        preprocessing.rayit(r)
        disp('Done rayit')
        preprocessing.write_netsw(r)
        disp(['Written netsw.inp.', r.expnr])
    end
    
    if (r.lEB || (r.iwalltemp ~= 1))
        preprocessing.generate_Tfacinit(r, r.lEB)
        preprocessing.write_Tfacinit(r)
        disp(['Written Tfacinit.inp.', r.expnr])
    end  
end




