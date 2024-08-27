clear all;

global ny nz nt dt_sig dy dz H u_H ltemp lmoist

expstr = '980';
DA_EXP = 'D:/Postdoc1/simulation/ecse1/experiments';
DA_WORK = 'D:/Postdoc1/simulation/ecse1/outputs';
DA_TOOLSDIR = 'D:/Postdoc1/simulation/ecse1/u-dales/tools';

nly = 4; nlz = 4; 
nlt = 25; nltu = 25; nltv = 25; nltw = 35;      %t_scale/dtmax

ylen = 1.8; zsize = 0.96; 
jtot = 64; ktot = 32;
dtmax = 0.02; runtime = 20.0;
ltemp = false; lmoist = false;

H = 0.16; u_H = 4.1;

ny = jtot+1; nz = ktot+1; nt = ceil(runtime/dtmax)+1; dt_sig = dtmax;
dy = ylen/jtot;
dz = zsize/ktot;

%%
addpath([DA_TOOLSDIR '/syntheticInflow']);
synInput_path = [DA_EXP '/' expstr '/syntheticInflow_inputs/'];
cd([DA_WORK '/' expstr]);
%%

fileID = fopen('z_driver.txt','r');
z_driver = fscanf(fileID,'%f');
fclose(fileID);

zm = z_driver;
zt = z_driver;

% zm = z_driver(1:end-1);
% for k = 1:length(zm)
%     zt(k) = (z_driver(k) + z_driver(k+1) ) / 2.0;
% end

%%

fileID = fopen('u_driver.txt','r');
u_driver = fscanf(fileID,'%f');
fclose(fileID);
u = (reshape(u_driver,ny,nz,nt));
clear u_driver

fileID = fopen('v_driver.txt','r');
v_driver = fscanf(fileID,'%f');
fclose(fileID);
v = (reshape(v_driver,ny,nz,nt));
clear v_driver

fileID = fopen('w_driver.txt','r');
w_driver = fscanf(fileID,'%f');
fclose(fileID);
w = (reshape(w_driver,ny,nz,nt));
clear w_driver

%%

profile_vel = readmatrix([synInput_path 'Reynolds_stress_profiles_velocity.txt'],'Range', 2);

postProcess_func.easy_statistics('upup',profile_vel(:,1),profile_vel(:,2),profile_vel(:,3),u,u,zt)
postProcess_func.easy_statistics('upvp',profile_vel(:,1),profile_vel(:,2),profile_vel(:,4),u,v,zt)
postProcess_func.easy_statistics('vpvp',profile_vel(:,1),profile_vel(:,2),profile_vel(:,5),v,v,zt)
postProcess_func.easy_statistics('upwp',profile_vel(:,1),profile_vel(:,2),profile_vel(:,6),u,w,zm)
postProcess_func.easy_statistics('vpwp',profile_vel(:,1),profile_vel(:,2),profile_vel(:,7),v,w,zm)
postProcess_func.easy_statistics('wpwp',profile_vel(:,1),profile_vel(:,2),profile_vel(:,8),w,w,zt)

%%

postProcess_func.plot_autocorr_t(nltu,u,'u');
postProcess_func.plot_autocorr_t(nltv,v,'v');
postProcess_func.plot_autocorr_t(nltw,w,'w');

%%

postProcess_func.plot_autocorr_y(nly,u,'u');
postProcess_func.plot_autocorr_y(nly,v,'v');
postProcess_func.plot_autocorr_y(nly,w,'w');

%%

postProcess_func.plot_autocorr_z(nlz,u,'u');
postProcess_func.plot_autocorr_z(nlz,v,'v');
postProcess_func.plot_autocorr_z(nlz,w,'w');

%%

if(ltemp)
    fileID = fopen('t_driver.txt','r');
    t_driver = fscanf(fileID,'%f');
    fclose(fileID);
    temp = (reshape(t_driver,ny,nz,nt));
    clear t_driver
    postProcess_func.plot_autocorr_t(nlt,temp,'t');
    postProcess_func.plot_autocorr_y(nly,temp,'t');
    postProcess_func.plot_autocorr_z(nlz,temp,'t');
end

%%

if(lmoist)
    fileID = fopen('q_driver.txt','r');
    q_driver = fscanf(fileID,'%f');
    fclose(fileID);
    qt = (reshape(q_driver,ny,nz,nt));
    clear q_driver
    postProcess_func.plot_autocorr_t(nlt,qt,'q');
    postProcess_func.plot_autocorr_y(nly,qt,'q');
    postProcess_func.plot_autocorr_z(nlz,qt,'q');
end
