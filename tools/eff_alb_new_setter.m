%% new effective albedo
% This script creates effective albedo based only on alphas and K*

%expnrs = [201:2:205,209:2:259];
expnrs = [201:2:259];
efalbs3 = [];
plan_area = 384*256;
I = 800;
zen = 45;
diff = 80; 
flat_tot_in = (I*cosd(zen)+diff)*plan_area;
for i = 1:length(expnrs)
    expnr = num2str(expnrs(i));
    fpath = ['/media/chris/Project3/uDALES2.0/experiments/original_200_259/' expnr];
    fpathout = ['/media/chris/Project4/outputs/' expnr];
    knet_file = [fpath '/netsw.inp.' expnr];
    fac_file = [fpath '/facets.inp.' expnr];
    type_file = [fpath '/factypes.inp.' expnr];
    facarea_file = [fpath '/facetarea.inp.' expnr];
    knet = dlmread(knet_file,'',1,0);
    facareas = dlmread(facarea_file,'',1,0);
    facs = dlmread(fac_file,'',1,0);
    fac_type = dlmread(type_file,'',3,0);
    albs = fac_type(:,5);
    typs = fac_type(:,1);
    tot_knet = sum(facareas.*knet);
    efalb3 = 1 - tot_knet/flat_tot_in;
    dlmwrite([fpathout '/efalb3.' expnr], efalb3);
    efalbs3 = [efalbs3, efalb3];
end 
figure()
xs = ones(size(efalbs3));
scatter(ones, efalbs3)