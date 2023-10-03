%% This script creates the necessarey inputs files for synthetic inflow turbulence generator
% It reads the tdump file from DA_SOURCE/expnr_source and then writes the reynolds stress and
% length and time scales as txt files inside DA_TARGET/expnr_target/syntheticInflow_inputs

clear all;
global column1 column2
global iplane
global stride
global nvec len

expstr_source = '515';
DA_SOURCE = 'D:/Postdoc1/simulation/testu2/outputs';
expstr_target = '799';
DA_TARGET = 'D:/Postdoc1/simulation/testu2/experiments';
iplane = 492;
column1 = 4; column2 = 7;
stride = 1;
ltempeq = true;
lmoist = true;
compute_length_time_scales_from_simulation_tdump = true;

%%

fileID = [DA_SOURCE '/' expstr_source '/tdump.' expstr_source '.nc'];

zt = ncread(fileID, 'zt');
zm = ncread(fileID, 'zm');
ut = ncread(fileID, 'ut');
upuptc = ncread(fileID, 'upuptc');
vpvptc = ncread(fileID, 'vpvptc');
wpwptc = ncread(fileID, 'wpwptc');
upvpt = ncread(fileID, 'upvpt');
upwpt = ncread(fileID, 'upwpt');
vpwpt = ncread(fileID, 'vpwpt');

nvec = length(zm);
len = length(zm)+1;

%%

fpath = [DA_TARGET '/' expstr_target '/'];
cd(fpath);
mkdir syntheticInflow_inputs
cd('syntheticInflow_inputs');

write_me(zm,ut,zt,'ut.txt')
write_me(zm,upuptc,zt,'upup.txt')
write_me(zm,vpvptc,zt,'vpvp.txt')
write_me(zm,wpwptc,zt,'wpwp.txt')
write_me(zm,upvpt,zt,'upvp.txt')
write_me(zm,upwpt,zm,'upwp.txt')
write_me(zm,vpwpt,zm,'vpwp.txt')

if (ltempeq)
    thlt = ncread(fileID, 'thlt');
    thlpthlpt = ncread(fileID, 'thlpthlpt');
    wpthlpt = ncread(fileID, 'wpthlpt');

    write_me(zm,thlt,zt,'thlt.txt')
    write_me(zm,thlpthlpt,zt,'thlpthlpt.txt')
    write_me(zm,wpthlpt,zm,'wpthlpt.txt')
end

if (lmoist)
    qtt = ncread(fileID, 'qtt');
    write_me(zm,qtt,zt,'qtt.txt')
end

%%
if (compute_length_time_scales_from_simulation_tdump)
    filename = 'length_time_scales_u.txt';
    lengthy = 4*ones(len,1);
    lengthz = 4*ones(len,1);
    ts = 0.5*ones(len,1);
    fileID = fopen(filename,'w');
    fprintf(fileID,'%8d %8d %15.10f\n',[lengthy lengthz ts]');
    fclose(fileID);

    filename = 'length_time_scales_v.txt';
    lengthy = 4*ones(len,1);
    lengthz = 4*ones(len,1);
    ts = 0.5*ones(len,1);
    fileID = fopen(filename,'w');
    fprintf(fileID,'%8d %8d %15.10f\n',[lengthy lengthz ts]');
    fclose(fileID);

    filename = 'length_time_scales_w.txt';
    lengthy = 4*ones(len,1);
    lengthz = 4*ones(len,1);
    ts = 0.5*sqrt(2)*ones(len,1);
    fileID = fopen(filename,'w');
    fprintf(fileID,'%8d %8d %15.10f\n',[lengthy lengthz ts]');
    fclose(fileID);

    if (ltempeq)  
        filename = 'length_time_scales_temp.txt';
        lengthy = 4*ones(len,1);
        lengthz = 4*ones(len,1);
        ts = 0.5*ones(len,1);
        fileID = fopen(filename,'w');
        fprintf(fileID,'%8d %8d %15.10f\n',[lengthy lengthz ts]');
        fclose(fileID);
    end

    if (lmoist)
        filename = 'length_time_scales_qt.txt';
        lengthy = 4*ones(len,1);
        lengthz = 4*ones(len,1);
        ts = 0.5*ones(len,1);
        fileID = fopen(filename,'w');
        fprintf(fileID,'%8d %8d %15.10f\n',[lengthy lengthz ts]');
        fclose(fileID);
    end
end

%%

function write_me(zm,var,ztzm,filename)

    global column1 column2
    global iplane
    global stride
    global nvec len

    var_m = 0;
    for i=column1:column2
        var_m = var_m + var(:,:,:,i);
    end
    var_m = var_m/(column2-column1+1);
    var_y_m = squeeze(mean(var_m,2));
    vec = var_y_m(iplane,:);
    vec = vec';
    
    vec1 = zeros(len,1);
    z = zeros(len,1);
    z(1:len-1) = zm;
    z(len) = zm(end) + (zm(end)-zm(end-1));
    vec1(2:len-1) = interp1(ztzm,vec,z(2:len-1), 'linear');
    vec1(len) = vec(nvec) + (( vec(nvec)-vec(nvec-1) )/( ztzm(nvec)-ztzm(nvec-1) ))*( z(len)-ztzm(nvec) );
    fileID = fopen(filename,'w');
    fprintf(fileID,'%15.10f\n',vec1(1:stride:end)');
    fclose(fileID);
%     figure
%     plot(vec,ztzm,vec1,z)
%     plottools
end