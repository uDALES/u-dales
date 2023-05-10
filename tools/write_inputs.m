%% write_inputs
% This script is run by the bash script da_inp.sh.
% It used to generate the necessary input files for uDALES.

expnr = '820';
ncpus = 1;

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

if (r.ltimedepsw || r.ltimedeplw || r.ltimedepnudge)
    preprocessing.readWeatherFile(r)
    if r.ltimedepnudge
       preprocessing.write_timedepnudge(r) 
    end
end

if r.nsv>0
    preprocessing.generate_scalar(r);
    preprocessing.write_scalar(r);
    disp(['Written scalar.inp.', r.expnr])
end

lgenerate_geometry = 1;

if ~r.lflat
    if lgenerate_geometry
    
    if ~r.lfloors
        if r.llidar
            disp('Generating geometry from LIDAR data')
            preprocessing.generate_topo_from_LIDAR(r)      
        elseif r.ltxtblocks
            disp('Generating geometry from text')
            preprocessing.generate_topo_from_txt(r)
            if r.ltxtmat
                disp('Generating materiality from text')
                preprocessing.generate_topomat_from_txt(r)
            end
                    
            
            
        elseif r.ltxtTG
            disp('Generating terrain and geometry from text files')
            preprocessing.generate_topo_from_txt_TG(r)
        else
            if (r.lstaggered || r.lstaggeredv || r.lcube || r.lcanyons)
                disp('Generating blocks from namoptions')
                preprocessing.generate_bl_from_namoptions(r)
                
%                 if (r.expnr == '006')
%                      r.bl((1:end/2),:) = [];
%                      r.bl((end - 12 + 1):end,:) = [];
%                 elseif (r.expnr == '007')
%                      r.bl((1:end/2),:) = [];
%                      r.bl((end - 12 + 1):end,:) = [];
%                 elseif (r.expnr == '008')
%                      r.bl((1:end/2),:) = [];
%                      r.bl((end - 12 + 1):end,:) = [];
%                 elseif (r.expnr == '009')
%                      r.bl((1:end/2),:) = [];
%                      r.bl((end - 13 + 1):end,:) = [];
%                  elseif(r.expnr == '109')
%                      r.bl((1:end/2),:) = [];
%                      r.bl((end - 12 + 1):end,:) = [];
%                  end
                
            elseif r.lblocksfile
                disp('Generating blocks from file')
                preprocessing.generate_bl_from_file(r) 
            end
            preprocessing.generate_topo_from_bl(r)        
        end
    
        if r.loverhangs
            disp('Generating overhangs')
            preprocessing.generate_overhangs(r)
        else
            preprocessing.addvar(r, 'overhangs', zeros(r.jtot, r.imax));
        end
                    if r.loverlayers
               disp('Generating overlayer')
               preprocessing.generate_overlayer(r)
            else                
                preprocessing.addvar(r, 'overlayer', zeros(r.jtot, r.imax)); % Do we need to do this?
                preprocessing.addvar(r, 'noverblocks', 0);
            end
        
        preprocessing.makeblocks(r)
        size(r.blocks)
        disp('made blocks')
        preprocessing.addvar(r, 'nunderblocks', r.nblocks);
        
        if r.loverlayers
            preprocessing.addvar(r, 'overblocks', [])
            preprocessing.addvar(r, 'noverblocks', 0)
            preprocessing.addvar(r, 'overbuildings', [])
            %preprocessing.addvar(r, 'noverbuildings', 0)
            preprocessing.addvar(r, 'layerheights', zeros(r.noverlayers,1))
            preprocessing.addvar(r, 'nblockslayers', zeros(r.noverlayers,1))
            for n=1:r.noverlayers
                preprocessing.makeblocks_overlayer(r, n)
            end
        end
        
        preprocessing.block2fac(r)
        disp('made facets')
        if r.loverlayers           
            preprocessing.addvar(r, 'overlayerfacets', [])
            preprocessing.addvar(r, 'noverlayerfacets', 0)
            for n=1:r.noverlayers
                preprocessing.block2fac_overlayer(r, n)
            end
        end
    else
        preprocessing.addvar(r, 'blocks', [])
        preprocessing.addvar(r, 'facets', [])
        preprocessing.addvar(r, 'nblockfcts', 0)
    end
    
    if (r.lEB && r.lboundingwalls)
        preprocessing.addboundingwalls(r)
    else
        preprocessing.addvar(r, 'boundingwallfacets', [])
        preprocessing.addvar(r, 'nboundingwallfacets', 0)
    end
    
    preprocessing.createfloors(r);
    disp('made floor facets');
    preprocessing.write_blocks(r)
    disp(['Written blocks.inp.', r.expnr])
       
    if r.ltxtmat
        preprocessing.set_walltypes_from_topomat(r)
    end
    
    % Temporary fix until have materials topo
    for n=1:r.nblockfcts
        if r.facets(n,end) == 0
            r.facets(n,2) = -1;
        end
    end
        
    preprocessing.write_facets(r)
    disp(['Written facets.inp.', r.expnr])
    
    %Not actually necessary if not using EB but useful for converting geometry formats.
    if ~(r.lfloors || r.lflat)
        preprocessing.generateFacetAreas(r, r.lEB)
        preprocessing.generateBlocksPhys(r);
        preprocessing.generateUniqueFacetCorners(r)
    end

    else
       load('geometry.mat')
       preprocessing.addvar(r, 'blocks', blocks);
       preprocessing.addvar(r, 'facets', facets);
       preprocessing.addvar(r, 'floorfacets', floorfacets);
       preprocessing.addvar(r, 'boundingwallfacets', boundingwallfacets);
       preprocessing.addvar(r, 'buildings', buildings);
       r.nblocks = nblocks;
       preprocessing.addvar(r, 'nblockstotal', nblockstotal);
       r.nfcts = nfcts;
       preprocessing.addvar(r, 'nfloorfacets', nfloorfacets);
       preprocessing.addvar(r, 'nboundingwallfacets', nboundingwallfacets);
       preprocessing.addvar(r, 'nbuildings', nbuildings);
       preprocessing.addvar(r, 'cornm', cornm);
       preprocessing.addvar(r, 'facetdim', facetdim)
       preprocessing.addvar(r, 'facetarea', facetarea)
       preprocessing.addvar(r, 'facetcorner', facetcorner)
       preprocessing.addvar(r, 'facetcornervis', facetcornervis)
       preprocessing.addvar(r, 'blocks_phys', blocks_phys)
       disp('Loaded geometry')
    end
    
    if isfile(['walltypes.inp.', expnr])
        r.walltypes = dlmread(['walltypes.inp.', expnr],'',3,0);
    else
        preprocessing.write_walltypes(r)
        disp(['Written walltypes.inp.', r.expnr])
    end
    
    lcalculate_direct_shortwave = 0;
    lcalculate_view_factors = 0;
    if r.lEB        
        if lcalculate_direct_shortwave
            disp('calculating direct shortwave')
            tic
            if ~(r.ltimedepsw || r.ltimedeplw || r.ltimedepnudge)
                preprocessing.vsolcAccelerated(r)
            else
                preprocessing.readWeatherFile(r)
                preprocessing.vsolcAccelerated_td(r)
            end
            toc
            disp('Done vsolc')
        else
            load('direct_shortwave.mat')
            preprocessing.addvar(r, 'Sdir_td', Sdir_td)           
            %preprocessing.addvar(r, 'asl', asl)
        end
        toc
        if lcalculate_view_factors
            preprocessing.vfcAccelerated(r)
            toc
            disp('Done vfc')
        else
            vf = ncread(['vf.nc.inp.' expnr], 'view factor');
            preprocessing.addvar(r, 'vf', vf);
            preprocessing.addvar(r, 'svf', max(1 - sum(vf, 2), 0));
        end
        
        preprocessing.write_svf(r)
        disp(['Written svf.inp.', r.expnr])
        preprocessing.write_vf(r)
        disp(['Written vf.nc.inp.', r.expnr])
        preprocessing.write_facetarea(r)
        disp(['Written facetarea.inp.', r.expnr])
        
        if ~(r.ltimedepsw)
            preprocessing.rayit(r)
            preprocessing.write_netsw(r)
        else
            preprocessing.rayit_td(r)
            preprocessing.write_timedepsw(r)
        end
        disp('Done rayit')
              
        disp(['Written netsw.inp.', r.expnr])
    else
        %preprocessing.generateFacetAreas(r, 0)
    end
    
    if (r.lEB || (r.iwalltemp ~= 1))
        preprocessing.generate_Tfacinit(r, 0)
        preprocessing.write_Tfacinit(r)
        disp(['Written Tfacinit.inp.', r.expnr])
    end  
end

toc


