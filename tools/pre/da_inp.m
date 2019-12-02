%% da_inp

close all
clear all

expnr = '009';
ncpus = 2;

DA_EXPDIR = getenv('DA_EXPDIR');
DA_PREDIR = getenv('DA_PREDIR');
addpath([DA_PREDIR '/']);
exppath = [DA_EXPDIR '/'];
cd([DA_EXPDIR '/' expnr])

r = da_pp(expnr, exppath);
da_pp.set_defaults(r, ncpus);

da_pp.generate_xygrid(r, true);
da_pp.generate_zgrid(r, true);
da_pp.generate_lscale(r, true)
da_pp.generate_prof(r, true);
da_pp.generate_scalar(r, true);

addpath(genpath(DA_PREDIR));
tempdir = [DA_PREDIR '/makefacets/temp'];
da_pp.generate_blocks(r, true);
da_pp.generate_trees(r, true);
da_pp.generate_purifs(r, true);
da_pp.plot_domain(r);
% plot vertical profiles
