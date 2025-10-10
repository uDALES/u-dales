
# Using uDALES utilities in MATLAB
This tutorial demonstrates how to use the uDALES MATLAB utilities for post-processing simulation data. We'll cover: 1. Merging short-term time-averaged data to long-term,  2. Spatial filtering (coarse-graining method). Specifically, including the functions:
   -   [**merge_stat_var**](#M_6c19)**.** This method merge long-term time average for a single variable and its variance. 
   -   [**merge_stat_cov**](#M_6c19)**.** This method merge long-term time average for two variables and their covariance. 
   -   [**coarse_graining**](#M_2c34)**.** This method using FFT to spatial averaging a field using a squared filter [1]. 
# Initialising udbase and loading data
```matlab
% preamble
clear variables
close all
% add the uDALES matlab path
% addpath('path_to_udales\tools\matlab')
addpath('path_to_udales\tools\matlab')
% create an instance of the udbase class
expnr = 110;  % Experiment number (corresponds to namoptions.110)
expdir = 'path_to_experiments\110';  % Path to simulation directory
sim = udbase(expnr, expdir);
```
```text
Warning: prof.inp.110 not found. Assuming equidistant grid.
```
<a id="M_6c19"></a>
# Combining short-term time average to long-term time average 
Load short-term time-averaged data, for example the 1-D plane average
```matlab
uxyt = sim.load_stat_xyt('uxyt');       % u-velocity profile (z) [m/s]
wxyt = sim.load_stat_xyt('wxyt');       % w-velocity profile [m/s]
upupxyt = sim.load_stat_xyt('upuptxyc'); % u-velocity variance [m²/s²]
upwpxyt = sim.load_stat_xyt('upwpxyt'); % u-w velocity covariance [m²/s²]
time = sim.load_stat_xyt('time');       % Time coordinate for xyt data [s]
```
Check the averaging time interval
```matlab
time
```
```text
time = 3x1 single column vector
1.0e+03
2.0001    
4.0001    
6.0001    
```
```matlab
length(time)
```
```text
ans = 3
```
Long-term time averaging for a single variable can be obtained by using
```matlab
% over entire time series
% for u-velocity only
[uxyt_longterm, ~] = merge_stat_var(uxyt, zeros(size(uxyt)), length(time)); % the last dimension for the input short-term average must be the time.
% for u-velocity and its variance
[uxyt_longterm, upupxyt_longterm] = merge_stat_var(uxyt, upupxyt, length(time));
% check the dimension
size(uxyt)
```
```text
ans = 1x2
256    3      
```
```matlab
size(uxyt_longterm)
```
```text
ans = 1x2
256    1      
```
```matlab
% over a specific time window
Nwindow = 2;  % Number of consecutive time steps to average together
% in this example, the average are taken using the last two time-series
% data, see function merge_stat_va for details.
[uxyt_longterm, upupxyt_longterm] = merge_stat_var(uxyt, upupxyt, Nwindow);
% check the dimension:
% The last dimension(time) has been averaged
size(uxyt)
```
```text
ans = 1x2
256    3      
```
```matlab
size(uxyt_longterm)
```
```text
ans = 1x2
256    1      
```
Long-term time averaging for variables and their covariance can be obtained by using
```matlab
% over entire time series
[uxyt_longterm, wxyt_longterm, upwpxyt_longterm] = merge_stat_cov(uxyt, wxyt, upwpxyt, length(time));
size(upwpxyt)
```
```text
ans = 1x2
256    3      
```
```matlab
size(upwpxyt_longterm)
```
```text
ans = 1x2
256    1      
```
```matlab
% over a specific time window
[uxyt_longterm, wxyt_longterm, upwpxyt_longterm] = merge_stat_cov(uxyt, wxyt, upwpxyt, Nwindow);
size(upwpxyt)
```
```text
ans = 1x2
256    3      
```
```matlab
size(upwpxyt_longterm)
```
```text
ans = 1x2
256    1      
```
<a id="M_2c34"></a>
# Coarse-graining Method
```matlab
% Apply spatial filters to smooth field data at different length scales
% This is useful for multi-scale analysis
ut = sim.load_stat_t('ut'); % load 3-D a field
% Use last two time interval to merge to one long-term time average
ut_longterm = merge_stat_var(ut, zeros(size(ut)), 2);
filter_lengths = [8, 32, 128];  % Filter widths in meters
u_filtered = coarse_graining(ut_longterm, filter_lengths, sim.dx, sim.xm, sim.ym);
```
```text
Applying coarse-graining filters...
  Filter 1/3 (Lflt=7.5m) completed
  Filter 2/3 (Lflt=32.5m) completed
  Filter 3/3 (Lflt=127.5m) completed
Coarse-graining completed in 1.16 seconds
```
```matlab
% Create figure with subplots comparing original and filtered fields
figure
k_level = 10;
% Original field (unfiltered)
subplot(2,2,1)
pcolor(sim.xt, sim.yt, ut(:,:,k_level)'); 
shading flat; axis equal tight; colorbar; 
clim([-1 2])
title('Original Field (Unfiltered)')
xlabel('x [m]'); ylabel('y [m]');
% Filtered fields with increasing filter lengths
filter_titles = {
    sprintf('Filter Length = %.0f m', filter_lengths(1)),
    sprintf('Filter Length = %.0f m', filter_lengths(2)), 
    sprintf('Filter Length = %.0f m', filter_lengths(3))
};
for i = 1:3
    subplot(2,2,i+1)
    pcolor(sim.xt, sim.yt, u_filtered(:,:,k_level,i)'); 
    shading flat; axis equal tight; colorbar; 
    clim([-1 2])
    title(filter_titles{i})
    xlabel('x [m]'); ylabel('y [m]');
end
```
![figure_0.png](udales-utility-tutorial_media/figure_0.png)
# References
[1] Maarten van Reeuwijk, Jingzi Huang (2025) Multi-scale Analysis of Flow over Heterogeneous Urban Environments, *Bound-Lay. Met.* **191**, 47.
