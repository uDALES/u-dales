% clear all
% close all

%% iterate shortwave and longwave radiation
%% iterate over every facet and over entire system until energy is balanced


%% derived quantities

[fct, wall] = loadfacets(expnr);
[sortt, sorti]=sort(F(:,1));  %sort by walltype

%% shortwave
%

%Z=35;          %zenith angle of the sun (could be function of time, location etc)
Dsky=zeros(nfcts,1);
Denv=zeros(nfcts,1);

albedo=zeros(nfcts,1);
emissivity=zeros(nfcts,1);

for i=1:nfcts
    j=find(fct.wlid(i)==wall.id);
    albedo(i)=wall.al(j);
    emissivity(i)=wall.em(j);
end

%isroof
isnotroof=ones(nfcts,1);
isnotroof(find(F(:,1)==1 & F(:,4)>0))=0;

%diffuse flux from sky and other walls
Kinnew=zeros(nfcts,1);

Kininit=(1-albedo).*(Sdir+Dsk.*svf);
Koutinit=albedo.*(Sdir+Dsk.*svf);

if ltestplot
figure
a=(1-albedo).*Sdir;
b=(1-albedo).*Dsk.*svf;
plot(1:nfcts,a(sorti),1:nfcts,b(sorti))
ax1 = gca; % current axes
set(ax1,'Xtick',1:1:nfcts)
xticklabels(sortt)
ax1_pos = ax1.Position; % position of first axes
xlabel('orientation')
ylabel('S and D')
ax2 = axes('Position',ax1_pos,...
    'YAxisLocation', 'right',...
    'XAxisLocation', 'top', ...
    'Color', 'None');
set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
set(ax2,'Xtick',1:1:nfcts)
set(ax2,'Xticklabels',sorti)
set(ax2,'FontSize',8)
xlabel('facet index')
end

%total radiation absorbed and reflected
sum(Kininit+Koutinit);
%total radiation absorbed
sum(Kininit);
%total radiation reflected
sum(Koutinit);



totincrease=zeros(10,1);
increase=zeros(10,1);
Kout=zeros(10,1);
Kout(2)=sum(Koutinit);


Kin=Kininit; %Kin §adds up
Kinold=Kininit;
Koutold=Koutinit; %Kout is wiped clean with every reflection
Koutnew=zeros(nfcts,1);
Kintemp=zeros(nfcts,1);
Kouttemp=zeros(nfcts,1);
count=0;
storerad=zeros(nfcts,nfcts);
Kiterin=zeros(300,1);
Kiterout=zeros(300,1);
itermaxdiff=zeros(300,1);
itermaxdiffloc=zeros(300,1);
itermaxrelloc=zeros(300,1);
moep=zeros(nfcts,20);
blub=zeros(nfcts,20);
moep(:,1)=Kininit;
while true
    count=count+1;
    for i=1:nfcts
        Kintemp(i)=0;
        Kouttemp(i)=0;
        for j=1:nfcts %sum all the incoming radiation on i, originally reflected from all the other j facets ("radiation reflected on j" x "what perecentage does i take of j's vision")
            inc=Koutold(j)*A(j)/A(i)*vf(j,i);  %[W/m2]
            storerad(i,j)=inc;
            Kintemp(i)=Kintemp(i)+(1-albedo(i))*inc;
            Kouttemp(i)=Kouttemp(i)+albedo(i)*inc;
        end
        Kin(i)=Kin(i)+Kintemp(i); %add newly absorbed radiation to already existing one
        Koutnew(i)=Kouttemp(i); %save newly reflected radiation for next iteration
    end
   % Kiterin(count)=Kin(1);
   % Kiterout(count)=Koutnew(1);
   % [itermaxdiff(count),itermaxdiffloc(count)]=max(Kintemp(i));
   % [~,itermaxrelloc(count)]=max(Koutnew./Koutold);
   % max(Koutnew./Koutold)
   % if all(Koutnew./Koutold<0.01)
   %     break
   % end
   moep(:,count+1)=Kin;
   if (max((Kin-Kinold)./Kinold)<0.01)
       disp(['reached relative tolerance after ' num2str(count) ' iterations'])
       break
   end
%    if (max(Koutnew-Koutold)<0.01 )
%        disp('reached absolute tolerance')
%        break
%    end
   Kinold=Kin;
   Koutold=Koutnew; %overwrite reflected radiation with new value 
end
%%
if lhqplot
scale=2;
scalef=1.5;
h= figure;
set(gcf,'units','centimeters','position',[0 0 14.5*scale 14.5*scale]);
set(h,'PaperPosition',[0 0 14.5*scale 14.5*scale]);
set(h,'PaperUnits','centimeters');
set(h,'renderer','painters');
select=140:-1:6;
select=[1:5 select];
hs1=subplot(1,2,1)
hold on
stylez={':x';':o';':d';':s';':+'};
cm=hsv;
colint=floor(length(cm)/5);
for k=1:length(select)
h=select(k);    
plot(0:1:count,moep(h,1:count+1),stylez{F(h,1)},'color',cm(1+F(h,1)*colint,:),'Linewidth',2)
end
xlim([0 count])
xlabel('iteration $n$','Interpreter','latex','FontSize',18)
ylabel('$K^{\downarrow}$ [Wm$^{-2}$]','Interpreter','latex','FontSize',18)
xticks(0:count)
xticklabels({'init' '1' '2' '3' '4' '5'})
set(hs1,'TickLabelInterpreter','latex')
set(hs1,'FontSize',12)
hs=subplot(1,2,2);
for r=2:count+1
blub(:,r)=((moep(:,r)-moep(:,r-1))./moep(:,r-1));
end
blub(:,1)=NaN;
hold on
for k=1:length(select)
h=select(k);  
plot(0:1:count,blub(h,1:count+1),stylez{F(h,1)},'color',cm(1+F(h,1)*colint,:),'Linewidth',2)
end
plot([0 count],[0.01 0.01],'k:','Linewidth',2)
xlim([0 count])
xlabel('iteration $n$','Interpreter','latex','FontSize',18)
xticks(0:count)
xticklabels({'init' '1' '2' '3' '4' '5'})
ylabel('$(K^{\downarrow}_n-K^{\downarrow}_{n-1})/K^{\downarrow}_{n-1}$ [-]','Interpreter','latex','FontSize',18)
hs.YScale='log';
ylim([0.0001 10])
yticklabels({'0.0001' '0.001' '0.01' '0.1' '1' '10'})
set(hs,'TickLabelInterpreter','latex')
set(hs,'FontSize',12)
set(gcf, 'Color', 'w');
legend('roof/road','west','east','north','south','Location','NorthEast')
export_fig rayswcov.eps
end
%%
if lwritefiles
    fname = [outputdir '/netsw.inp.' num2str(expnr)];
    fileID = fopen(fname,'w');
    fprintf(fileID,'# %4s\n','net shortwave on facets [W/m2] (including reflections and diffusive)');
    fprintf(fileID,'%6d\n', Kin);
    fclose(fileID);
end

if lhqplot
    scale=2;
    scalef=1.5;
    h= figure;
    set(gcf,'units','centimeters','position',[0 0 14.5*scale 14.5*scale]);
    set(h,'PaperPosition',[0 0 14.5*scale 14.5*scale]);
    set(h,'PaperUnits','centimeters');
    set(h,'renderer','painters');
    a=(1-albedo).*Sdir;
    b=(1-albedo).*Dsk.*svf;
    plot(1:nfcts,Kin(sorti),'x',1:nfcts,Kininit(sorti),'+',1:nfcts,a(sorti),'o',1:nfcts,b(sorti),'d')
    legend('Total K','Initial K','Initial Sdir','Initial D')
    ax1 = gca; % current axes
    set(ax1,'Xtick',1:1:nfcts)
    set(ax1,'FontSize',9)
    xticklabels(sortt)
    ax1_pos = ax1.Position; % position of first axes
    xlabel('orientation','Interpreter','latex','FontSize',12)
    ylabel('Radiation [$W/m^2$]','Interpreter','latex','FontSize',12)
    ax2 = axes('Position',ax1_pos,...
        'YAxisLocation', 'right',...
        'XAxisLocation', 'top', ...
        'Color', 'None');
    set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
    set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
    set(ax2,'YTickLabel','')
    set(ax2,'Xtick',1:3:nfcts)
    set(ax2,'Xticklabels',sorti(1:3:nfcts),'TickLabelInterpreter','latex')
    set(ax2,'FontSize',9)
    xl2=xlabel('facet index','FontSize',12);
    xl2.Position(2)=xl2.Position(2)+15;
    offset = repmat(ax2.YTick(end)+20,1,numel(ax2.XTick));
    % create new lables:
   % text(ax1.XTick(2:3:nfcts),offset,num2str(sorti(2:3:nfcts)),'HorizontalAlign','center','FontSize',9,'Interpreter','latex')
    
    offset = repmat(ax2.YTick(end)+30,1,numel(ax2.XTick)-1);
    % create new lables:
 %   text(ax1.XTick(3:3:nfcts),offset,num2str(sorti(3:3:nfcts)),'HorizontalAlign','center','FontSize',9,'Interpreter','latex')
    set(gcf, 'Color', 'w');
    export_fig shortwaveit.eps
    
    
    scale=2;
    scalef=1.5;
    h= figure;
    set(gcf,'units','centimeters','position',[0 0 14.5*scale 14.5*scale]);
    set(h,'PaperPosition',[0 0 14.5*scale 14.5*scale]);
    set(h,'PaperUnits','centimeters');
    set(h,'renderer','painters');
    a=(1-albedo).*Sdir;
    b=(1-albedo).*Dsk.*svf;
    plot(1:nfcts,blub(sorti,2),'x',1:nfcts,blub(sorti,3),'+',1:nfcts,blub(sorti,4),'o',1:nfcts,blub(sorti,5),'d',1:nfcts,blub(sorti,6),'s')
    legend('Total K','Initial K','Initial Sdir','Initial D')
    ax1 = gca; % current axes
    ax1.YScale='log';
    set(ax1,'Xtick',1:1:nfcts)
    set(ax1,'FontSize',9)
    xticklabels(sortt)
    ax1_pos = ax1.Position; % position of first axes
    xlabel('orientation','Interpreter','latex','FontSize',12)
    ylabel('Radiation [$W/m^2$]','Interpreter','latex','FontSize',12)
    ax2 = axes('Position',ax1_pos,...
        'YAxisLocation', 'right',...
        'XAxisLocation', 'top', ...
        'Color', 'None');
    set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
    set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
    set(ax2,'YTickLabel','')
    %set(ax2,'Xtick',1:3:nfcts)
    %set(ax2,'Xticklabels',sorti(1:3:nfcts),'TickLabelInterpreter','latex')
    set(ax2,'Xticklabels','','TickLabelInterpreter','latex')
    set(ax2,'FontSize',9)
    xl2=xlabel('facet index','FontSize',12);
    xl2.Position(2)=xl2.Position(2)+15;
    offset = repmat(ax2.YTick(end)+20,1,numel(ax2.XTick));
    % create new lables:
   % text(ax1.XTick(2:3:nfcts),offset,num2str(sorti(2:3:nfcts)),'HorizontalAlign','center','FontSize',9,'Interpreter','latex')
    
    offset = repmat(ax2.YTick(end)+30,1,numel(ax2.XTick)-1);
    % create new lables:
 %   text(ax1.XTick(3:3:nfcts),offset,num2str(sorti(3:3:nfcts)),'HorizontalAlign','center','FontSize',9,'Interpreter','latex')
    set(gcf, 'Color', 'w');
    export_fig shortwaveit2.eps

    % ax3 = axes('Position',ax1_pos,...
    %     'YAxisLocation', 'right',...
    %     'XAxisLocation', 'top', ...
    %     'Color', 'None');
    % set (ax3, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
    % set (ax3, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
    % ax3.YTickLabel='';
    % ax3.YTick=''
    % set(ax3,'Xtick',2:2:nfcts)
    % set(ax3,'Xticklabels',sorti(2:2:nfcts))
    % set(ax3,'FontSize',8)
    % box(ax3(1),'off')
    % ax3.TickLength=[0;0]
    %ax3.Position(4)=ax3.Position(4)*1.02;
    
    
    
    
    
    figure
    plot(Kin(sorti)./Kininit(sorti))
    ax1 = gca; % current axes
    set(ax1,'Xtick',1:1:nfcts)
    xticklabels(sortt)
    ax1_pos = ax1.Position; % position of first axes
    xlabel('orientation')
    ylabel('K/Kinit')
    ax2 = axes('Position',ax1_pos,...
        'YAxisLocation', 'right',...
        'XAxisLocation', 'top', ...
        'Color', 'None');
    set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
    set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
    set(ax2,'Xtick',1:1:nfcts)
    set(ax2,'Xticklabels',sorti)
    set(ax2,'FontSize',8)
    xlabel('facet index')
    
    figure
    plot(Kin(sorti)-Kininit(sorti))
    ax1 = gca; % current axes
    set(ax1,'Xtick',1:1:nfcts)
    xticklabels(sortt)
    ax1_pos = ax1.Position; % position of first axes
    xlabel('orientation')
    ylabel('K-Kinit')
    ax2 = axes('Position',ax1_pos,...
        'YAxisLocation', 'right',...
        'XAxisLocation', 'top', ...
        'Color', 'None');
    set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
    set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
    set(ax2,'Xtick',1:1:nfcts)
    set(ax2,'Xticklabels',sorti)
    set(ax2,'FontSize',8)
    xlabel('facet index')
    
end
%%





%% Inital Temperature and longwave
%solve energy budget equation in an initial steady state
%assume 0 wall heatflux
%assume 0 latent heatflux
%assume constant air temperature and heat transfer coefficient
%assume constant longwave
%solve for Tinit and Linit
%assume absorbtivity is equal to emissivity
%=>
%K+L=H
Tair=300;  %K, air temperature
Tinitial=300; %K, initial facet temperature

% Tinc=1*ones(nfcts,1); %incremental temperature change of facets if not in equilibrium
% tolerance=2.5;  %W/m2, if below this threshold, change will be made to facet temperature

Tinc=2*ones(nfcts,1); %incremental temperature change of facets if not in equilibrium
tolerance=2.5*Tinc(1);  %W/m2, if below this threshold, change will be made to facet temperature
%the two are somewhat related, an increase of 1K will result in a longwave
%change of approximately 5W/m2, e.g.:
%if we are 4.56K off we correct approximately  (1-5/22.8)*2 = 1.5616
%now we are 3K off and correct approximately for (1-5/15)*2 = 1.3333
%now we are 1.6666 off and we correct for 0.8 (had we corrected for 2 at
%this step, we would have overshot and the solution would oscilate)
%now we are satisfied at 0.866 off which is less than 1degree
%if for any


Tterminate=0; %Terminate if the absolute temperature change between iterations is equal or below this value
%somewhat redundant since this basically means all the facets are within the tolerance


absorptivity=emissivity;


sigma=5.67e-8;
Lsk=350;
hc=0;  %(rho*cp)/R    ~(1.2*1000)/100=12   R~100

Tinit=ones(nfcts,1)*Tinitial;
Lsky=Lsk.*svf;

%Loutinit=emissivity.*sigma.*Tinit.^4;
Told=Tinit;
Tnew=zeros(nfcts,100);
Lin=zeros(nfcts,1);
b
k=0;
figure
plot(1:nfcts,Told)
hold on
while true  %what happens to the reflected longwave? i.e. (1-absorptivity)*Lin
    %for k=1:10
    %count=count+1;
    
    k=k+1;
    %maxchange=0;
    change=0;
    Lin(:)=0;
    
    for i=1:nfcts
        for j=1:nfcts %sum all the incoming radiation on i, originally reflected from all the other j facets ("radiation reflected on j" x "what perecentage does i take of j's vision")
            inc=Told(j)^4*emissivity(j)*sigma*vf(j,i)*A(j)/A(i);
            Lin(i)=Lin(i)+inc;
        end
        
        %calculate energy balance
        eb=Kin(i)+absorptivity(i)*(Lin(i)+Lsky(i))-hc*(Told(i)-Tair)-emissivity(i)*sigma*Told(i)^4;
        
        if eb<-tolerance %if energy balance is negative, facet is too hot
            
            Tnew(i,k)=Told(i)-Tinc(i)*(1+tolerance/eb); %remove incremental temperature, scale by deviation from accepted tolerance in W/m2
            
        elseif eb>tolerance
            Tnew(i,k)=Told(i)+Tinc(i)*(1-tolerance/eb);
        else %do nothing, within tolerance
            Tnew(i,k)=Told(i);
        end
        
        % maxchange=max(maxchange,abs(Tnew(i,k)-Told(i)));
        change=change+abs(Tnew(i,k)-Told(i));
        
    end
    
    Told=Tnew(:,k);
    
    plot(1:nfcts,Tnew(:,k),'-x')
    
    if change<=Tterminate%maxchange<=Tterminate
        change
        k
        break
    elseif k>99
        change
        k
        break
    end
end
xlabel('facet NR')
ylabel('facet temperature')

fname = [outputdir '/Tfacinit.inp.' num2str(expnr)];
fileID = fopen(fname,'W');
fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium\n');
fclose(fileID);
dlmwrite(fname,Tnew(:,k),'-append','delimiter',' ','precision','%4f')

 fname = [outputdir '/Tfacinitnudged.inp.' num2str(expnr)];
 fileID = fopen(fname,'W');
 fprintf(fileID, '# Initial facet tempereatures in radiative equilibrium nudged to 288\n');
 fclose(fileID);
 blub=Tnew(Tnew>0);
 blublub=mean(blub(:));
 blub=Tnew-(Tnew-blublub)*0.5;
 blub(Tnew==0)=0;
 dlmwrite(fname,blub,'-append','delimiter',' ','precision','%4f')


if ltestplot
    figure
    hold on
    plot(Tnew(sorti,k)-Tinit(sorti),'r-x')
    plot(blub(sorti,k)-Tinit(sorti),'b-x')
    xlabel('facet NR')
    ylabel('T_{end}-T_{init}')
    ax1 = gca; % current axes
    set(ax1,'Xtick',1:1:nfcts)
    xticklabels(sortt)
    ax1_pos = ax1.Position; % position of first axes
    xlabel('orientation')
    ylabel('T_{new}-T_{init}')
    ax2 = axes('Position',ax1_pos,...
        'YAxisLocation', 'right',...
        'XAxisLocation', 'top', ...
        'Color', 'None');
    set (ax2, 'XLim', get (ax1, 'XLim'), 'Layer', 'top');
    set (ax2, 'YLim', get (ax1, 'YLim'), 'Layer', 'top');
    set(ax2,'Xtick',1:1:nfcts)
    set(ax2,'Xticklabels',sorti)
    set(ax2,'FontSize',8)
    xlabel('facet index')
    
    
    
    figure
    hold on
    plot(Tnew(15,1:k),'-x')
    plot(Tnew(41,1:k),'-x')
    plot(Tnew(51,1:k),'-x')
    plot(Tnew(58,1:k),'-x')
    plot(Tnew(68,1:k),'-x')
    xlabel('iteration')
    ylabel('facet temperature')
    legend('15','41','51','58','68')
end
