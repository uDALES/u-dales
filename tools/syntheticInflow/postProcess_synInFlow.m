clear all;

global ny nz nt dt_sig dy dz H u_H ltemp lmoist

expstr = '799';
DA_EXPDIR = 'D:/Postdoc1/simulation/testu2/experiments';
DA_TOOLSDIR = 'D:/Postdoc1/simulation/testu2/u-dales/tools';

nly = 4; nlz = 4; nlt = 25; nltw = 35;
dy = 1.8/128; dz = 0.96/128;
H = 0.16; u_H = 4.3;
ltemp = true; lmoist = true;

ny = 257; nz = 129; nt = 2501; dt_sig = 0.02;  % jtot+1; ktot+1; ceil(runtime/dtmax)+1; dtmax;

%%
addpath([DA_TOOLSDIR '/syntheticInflow']);
synInput_path = [DA_EXPDIR '/' expstr '/syntheticInflow_inputs/'];
cd([DA_EXPDIR '/' expstr]);
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

postProcess_func.easy_statistics(synInput_path,'upup.txt',z_driver,u,u,zt)
postProcess_func.easy_statistics(synInput_path,'upvp.txt',z_driver,u,v,zt)
postProcess_func.easy_statistics(synInput_path,'upwp.txt',z_driver,u,w,zm)
postProcess_func.easy_statistics(synInput_path,'vpvp.txt',z_driver,v,v,zt)
postProcess_func.easy_statistics(synInput_path,'vpwp.txt',z_driver,v,w,zm)
postProcess_func.easy_statistics(synInput_path,'wpwp.txt',z_driver,w,w,zt)

%%

postProcess_func.plot_autocorr_t(nlt,u,'u');
postProcess_func.plot_autocorr_t(nlt,v,'v');
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
