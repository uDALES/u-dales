%% new effective albedo
% This script creates effective albedo based only on alphas and K*

%expnrs = [201:2:205,209:2:259];
expnrs = [201:2:259];
efalbs2 = [];
for i = 1:length(expnrs)
    expnr = num2str(expnrs(i));
    fpath = ['/media/chris/Project3/uDALES2.0/experiments/' expnr];
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
    tot_knet = 0;
    tot_kin = 0;
    for j = 1:length(facs(:,1))
        type = facs(j,1);
        alb_in = find(typs==type);
        alb = albs(alb_in);
        tot_knet = tot_knet+facareas(j)*knet(j);
        tot_kin = tot_kin + facareas(j)*knet(j)/(1-alb);
    end
    efalb2 = 1 - tot_knet/tot_kin;
    dlmwrite([fpathout '/efalb2.' expnr], efalb2);
    efalbs2 = [efalbs2, efalb2];
end 
figure()
xs = ones(size(efalbs2));
scatter(ones, efalbs2)
