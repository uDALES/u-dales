classdef postProcess_func
    methods(Static)

        function plot_autocorr_t(n_lag_t,vel,vel_name)
            global ny nz dt_sig
            n_lag = 3*n_lag_t;
            
            set(groot,'defaultAxesTickLabelInterpreter','latex');  
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
            figure
            hold on
            
            counter = 0;
            acf_u = [];
            for j=2:ny
                for k=2:nz
                    acf_u_1 = autocorr(squeeze(vel(j,k,:)),'NumLags',n_lag);
                    if j==2 && k==2
                        acf_u = acf_u_1;
                    else
                        acf_u = acf_u + acf_u_1;
                    end
        %             plot(dt_sig*(0:1:n_lag),acf_u_1)
                    counter = counter + 1;
                end
            end
            acf_u = acf_u/counter;
            
            plot(dt_sig*(0:1:n_lag),exp(-pi*(0:1:n_lag)/(2*n_lag_t)),'LineWidth',2,'Color','blue')
            plot(dt_sig*(0:1:n_lag),acf_u,'--','LineWidth',2,'Color','black')
            set(gca,"Box","on")
            set(gca,'FontSize',18.0)
            set(gca,'LineWidth',3)
            xlabel('$t_{lag}$','FontSize',24)
            ylabel(['$R_{' vel_name vel_name '}$'],'FontSize',24)
            legend('Analytical', 'Synthetic inflow')
            title(['Time correlation : $R_{' vel_name vel_name '}$'])
        end

        
        function plot_autocorr_y(n_lag_y,vel,vel_name)
            global nz nt dy
            
            n_lag = 3*n_lag_y;
            
            set(groot,'defaultAxesTickLabelInterpreter','latex');  
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
            figure
            hold on
            
            counter = 0;
            acf_u = [];
            for j=2:nz
                for k=2:nt
                    acf_u_1 = autocorr(squeeze(vel(:,j,k)),'NumLags',n_lag);
                    if j==2 && k==2
                        acf_u = acf_u_1;
                    else
                        acf_u = acf_u + acf_u_1;
                    end
                    counter = counter + 1;
                end
            end
            acf_u = acf_u/counter;
            
            plot(dy*(0:1:n_lag),exp(-pi*(0:1:n_lag)/(2*n_lag_y)),'LineWidth',2,'Color','blue')
            plot(dy*(0:1:n_lag),acf_u,'--','LineWidth',2,'Color','black')           
            set(gca,"Box","on")
            set(gca,'FontSize',18.0)
            set(gca,'LineWidth',3)
            xlabel('$y_{lag}$','FontSize',24)
            ylabel(['$R_{' vel_name vel_name '}$'],'FontSize',24)
            legend('Analytical', 'Synthetic inflow')
            title(['Spatial correlation along y: $R_{' vel_name vel_name '}$'])
        end
        
        
        function plot_autocorr_z(n_lag_z,vel,vel_name)
            global ny nz nt dz
            
            n_lag = 3*n_lag_z;
            
            if vel_name=='u' || vel_name=='t' || vel_name=='q'
                vel_t = squeeze(mean(vel,3));
                vel_ty = squeeze(mean(vel_t,1));
                for k = 1:nz
                    vel(:,k,:) = vel(:,k,:) - vel_ty(k);
                end
            end
            
            set(groot,'defaultAxesTickLabelInterpreter','latex');  
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
            figure
            hold on
            
            counter = 0;
            acf_u = [];
            for j=2:ny
                for k=2:nt
                    acf_u_1 = autocorr(squeeze(vel(j,:,k)),'NumLags',n_lag);
                    if j==2 && k==2
                        acf_u = acf_u_1;
                    else
                        acf_u = acf_u + acf_u_1;
                    end
                    counter = counter + 1;
                end
            end
            acf_u = acf_u/counter;
            
            plot(dz*(0:1:n_lag),exp(-pi*(0:1:n_lag)/(2*n_lag_z)),'LineWidth',2,'Color','blue')
            plot(dz*(0:1:n_lag),acf_u,'--','LineWidth',2,'Color','black')           
            set(gca,"Box","on")
            set(gca,'FontSize',18.0)
            set(gca,'LineWidth',3)
            xlabel('$z_{lag}$','FontSize',24)
            ylabel(['$R_{' vel_name vel_name '}$'],'FontSize',24)
            legend('Analytical', 'Synthetic inflow')
            title(['Spatial correlation along z: $R_{' vel_name vel_name '}$'])
        end

        
        function easy_statistics(filename,zoriginal,Uoriginal,vpvp_original,vel1,vel2,z)
            global H u_H
            
            set(groot,'defaultAxesTickLabelInterpreter','latex');  
            set(groot,'defaulttextinterpreter','latex');
            set(groot,'defaultLegendInterpreter','latex');
            
            if(filename(1)=='u')
                fprintf('Make sure you have provided u data as the first velocity input.\n');
                
                vel1_t = squeeze(mean(vel1,3));
                vel1_ty = squeeze(mean(vel1_t,1));
                for k = 1:length(z)
                    vel1(:,k,:) = vel1(:,k,:) - vel1_ty(k);
                end
                
                if(filename(3)=='u')
                    
                    vel2_t = squeeze(mean(vel2,3));
                    vel2_ty = squeeze(mean(vel2_t,1));
                    for k = 1:length(z)
                        vel2(:,k,:) = vel2(:,k,:) - vel2_ty(k);
                    end

                    figure
                    hold on
                    plot(Uoriginal/u_H,zoriginal/H,'LineWidth',2,'Color','blue')
                    plot(vel2_ty/u_H,z/H,'--','LineWidth',2,'Color','black')
                    xlim([0 2])
                    ylim([0 6])
                    set(gca,"Box","on")
                    set(gca,'FontSize',18.0)
                    set(gca,'LineWidth',3)
                    xlabel(['$\overline{' filename(1) '}/U_H$'],'FontSize',24)
                    ylabel('$z/H$','FontSize',24)
                    legend('Original input', 'Synthetic inflow','Location','northwest')
                    title(['Mean inflow : $\overline{' filename(1) '}$'])
                    
                    figure
                    pcolor(vel2_t')
                    shading interp

                    figure
                    contourf(vel2_t')
                end           
            end
            
            vpvp_computed = vel1.*vel2;
            vpvp_t_computed = squeeze(mean(vpvp_computed,3));
            vpvp_ty_computed = squeeze(mean(vpvp_t_computed,1));
            
            figure
            hold on
            if filename(1) == filename(3)
                plot(vpvp_original,zoriginal/H,'LineWidth',2,'Color','blue')
                plot(vpvp_ty_computed,z/H,'--','LineWidth',2,'Color','black')
            else
                plot(-vpvp_original,zoriginal/H,'LineWidth',2,'Color','blue')
                plot(-vpvp_ty_computed,z/H,'--','LineWidth',2,'Color','black')
            end
            set(gca,"Box","on")
            set(gca,'FontSize',18.0)
            set(gca,'LineWidth',3)
            xlabel(['$\overline{' filename(1) '^{\prime}' filename(3) '^{\prime}}$'],'FontSize',24)
            ylabel('$z/H$','FontSize',24)
            legend('Original input', 'Synthetic inflow')
            title(['Easy statistics : $\overline{' filename(1) '^{\prime}' filename(3) '^{\prime}}$'])  
        end
        
    end 
end