%% This script creates the necessarey inputs files for synthetic inflow turbulence generator
% It reads the tdump file from DA_SOURCE/expnr_source and then writes the reynolds stress and
% length and time scales as txt files inside DA_TARGET/expnr_target/syntheticInflow_inputs

clear all; clc;

set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

global column1 column2
global iplane
global stride
global nvec len
global z

%%
expstr_source = '515';
DA_SOURCE = 'D:/Postdoc1/simulation/testu2/outputs';
expstr_target = '980';
DA_TARGET = 'D:/Postdoc1/simulation/ecse1/experiments';
iplane = 492;
column1 = 4; column2 = 7;
stride = 4;
ltempeq = false;
lmoist = false;
lplot_profiles = true;
lcalc_time_and_length_scale = false;

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

z = zeros(len,1);
z(1:len-1) = zm;
z(len) = zm(end) + (zm(end)-zm(end-1));

%%

fpath = [DA_TARGET '/' expstr_target '/'];
cd(fpath);
mkdir syntheticInflow_inputs
cd('syntheticInflow_inputs');

u = obtain_profile(ut,zt);
upup = obtain_profile(upuptc,zt);
vpvp = obtain_profile(vpvptc,zt);
wpwp = obtain_profile(wpwptc,zt);
upvp = obtain_profile(upvpt,zt);
upwp = obtain_profile(upwpt,zm);
vpwp = obtain_profile(vpwpt,zm);

profile_vel = [z(1:stride:end) u(1:stride:end) upup(1:stride:end) upvp(1:stride:end) vpvp(1:stride:end) upwp(1:stride:end) vpwp(1:stride:end) wpwp(1:stride:end)];

filename = fopen('Reynolds_stress_profiles_velocity.txt','w');
fprintf(filename,"z u u'u' u'v' v'v' u'w' v'w' w'w'\n");
fprintf(filename,'%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n',profile_vel');
fclose(filename);


if (ltempeq)
    thlt = ncread(fileID, 'thlt');
    thlpthlpt = ncread(fileID, 'thlpthlpt');
    wpthlpt = ncread(fileID, 'wpthlpt');

    thl = obtain_profile(thlt,zt);
    thl(1) = thl(2);
    thlpthlp = obtain_profile(thlpthlpt,zt);
    wpthlp = obtain_profile(wpthlpt,zm);

    upthlp = wpthlp;
    vpthlp = wpthlp;

    profile_temp = [z(1:stride:end) thl(1:stride:end) upthlp(1:stride:end) vpthlp(1:stride:end) wpthlp(1:stride:end) thlpthlp(1:stride:end)];

    filename = fopen('Reynolds_stress_profiles_temp.txt','w');
    fprintf(filename,"z thl u'thl' v'thl' w'thl' thl'thl'\n");
    fprintf(filename,'%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n',profile_temp');
    fclose(filename);
end

if (lmoist)
    qtt = ncread(fileID, 'qtt');
    qt = obtain_profile(qtt,zt);

    upqtp = wpthlp;
    vpqtp = upqtp;
    wpqtp = upqtp;
    thlpqtp = upqtp;
    qtpqtp = upqtp;

    profile_moist = [z(1:stride:end) qt(1:stride:end) upqtp(1:stride:end) vpqtp(1:stride:end) wpqtp(1:stride:end) thlpqtp(1:stride:end) qtpqtp(1:stride:end)];

    filename = fopen('Reynolds_stress_profiles_moist.txt','w');
    fprintf(filename,"z qt u'q' v'q' w'q' thl'q' q'q'\n");
    fprintf(filename,'%15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n',profile_moist');
    fclose(filename);
end

%%
if (lplot_profiles)
    profile_info{1} = '$z$ [m]';
    profile_info{2} = '$\overline{u}$';
    profile_info{3} = '$\overline{u^{\prime}u^{\prime}}$';
    profile_info{4} = '$\overline{u^{\prime}v^{\prime}}$';
    profile_info{5} = '$\overline{v^{\prime}v^{\prime}}$';
    profile_info{6} = '$\overline{u^{\prime}w^{\prime}}$';
    profile_info{7} = '$\overline{v^{\prime}w^{\prime}}$';
    profile_info{8} = '$\overline{w^{\prime}w^{\prime}}$';
    
    plot_profiles(profile_info,profile_vel)
    clear profile_info
    

    if (ltempeq)
        profile_info{1} = '$z$ [m]';
        profile_info{2} = '$\overline{\theta}$';
        profile_info{3} = '$\overline{u^{\prime}\theta^{\prime}}$';
        profile_info{4} = '$\overline{v^{\prime}\theta^{\prime}}$';
        profile_info{5} = '$\overline{w^{\prime}\theta^{\prime}}$';
        profile_info{6} = '$\overline{\theta^{\prime}\theta^{\prime}}$';
        plot_profiles(profile_info,profile_temp)
        clear profile_info
    end

    if (lmoist)
        profile_info{1} = '$z$ [m]';
        profile_info{2} = '$\overline{q}$';
        profile_info{3} = '$\overline{u^{\prime}q^{\prime}}$';
        profile_info{4} = '$\overline{v^{\prime}q^{\prime}}$';
        profile_info{5} = '$\overline{w^{\prime}q^{\prime}}$';
        profile_info{6} = '$\overline{\theta^{\prime}q^{\prime}}$';
        profile_info{7} = '$\overline{q^{\prime}q^{\prime}}$';
        plot_profiles(profile_info,profile_moist)
        clear profile_info
    end
end

%%
if (~lcalc_time_and_length_scale)
    filename = 'length_time_scales_u.txt';
    lengthy = 4*ones(len,1);
    lengthz = 4*ones(len,1);
    ts = 0.5*ones(len,1);
    fileID = fopen(filename,'w');
    fprintf(fileID,'z nLy nLz T\n');
    fprintf(fileID,'%15.10f %8d %8d %15.10f\n',[z(1:stride:end) lengthy(1:stride:end) lengthz(1:stride:end) ts(1:stride:end)]');
    fclose(fileID);

    filename = 'length_time_scales_v.txt';
    lengthy = 4*ones(len,1);
    lengthz = 4*ones(len,1);
    ts = 0.5*ones(len,1);
    fileID = fopen(filename,'w');
    fprintf(fileID,'z nLy nLz T\n');
    fprintf(fileID,'%15.10f %8d %8d %15.10f\n',[z(1:stride:end) lengthy(1:stride:end) lengthz(1:stride:end) ts(1:stride:end)]');
    fclose(fileID);

    filename = 'length_time_scales_w.txt';
    lengthy = 4*ones(len,1);
    lengthz = 4*ones(len,1);
    ts = 0.5*sqrt(2)*ones(len,1);
    fileID = fopen(filename,'w');
    fprintf(fileID,'z nLy nLz T\n');
    fprintf(fileID,'%15.10f %8d %8d %15.10f\n',[z(1:stride:end) lengthy(1:stride:end) lengthz(1:stride:end) ts(1:stride:end)]');
    fclose(fileID);

    if (ltempeq)  
        filename = 'length_time_scales_temp.txt';
        lengthy = 4*ones(len,1);
        lengthz = 4*ones(len,1);
        ts = 0.5*ones(len,1);
        fileID = fopen(filename,'w');
        fprintf(fileID,'z nLy nLz T\n');
        fprintf(fileID,'%15.10f %8d %8d %15.10f\n',[z(1:stride:end) lengthy(1:stride:end) lengthz(1:stride:end) ts(1:stride:end)]');
        fclose(fileID);
    end

    if (lmoist)
        filename = 'length_time_scales_qt.txt';
        lengthy = 4*ones(len,1);
        lengthz = 4*ones(len,1);
        ts = 0.5*ones(len,1);
        fileID = fopen(filename,'w');
        fprintf(fileID,'z nLy nLz T\n');
        fprintf(fileID,'%15.10f %8d %8d %15.10f\n',[z(1:stride:end) lengthy(1:stride:end) lengthz(1:stride:end) ts(1:stride:end)]');
        fclose(fileID);
    end
end

%%

function profile = obtain_profile(var,ztzm)

    global column1 column2
    global iplane
    global nvec len
    global z

    var_y_m = squeeze( mean( squeeze( mean( var(:,:,:,column1:column2), 4 ) ), 2 ) );
    vec = var_y_m(iplane,:);
    vec = vec';
    
    profile = zeros(len,1);
    profile(2:len-1) = interp1(ztzm,vec,z(2:len-1), 'linear');
    profile(len) = vec(nvec) + (( vec(nvec)-vec(nvec-1) )/( ztzm(nvec)-ztzm(nvec-1) ))*( z(len)-ztzm(nvec) );
    
    % figure
    % plot(vec,ztzm,profile,z)
    % plottools
end

function plot_profiles(profile_info,profile)
    for i=2:length(profile(1,:))
        figure
        set(gca,"Box","on")
        set(gca,'FontSize',20.0)
        set(gca,'LineWidth',2.5)
        grid('on')
        hold on
        plot(profile(:,i),profile(:,1),'LineWidth',2)
        xlabel(profile_info{i},'FontSize',30)
        ylabel(profile_info{1},'FontSize',30)
        axis tight
        plottools
    end
end