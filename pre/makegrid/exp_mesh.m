% Routine to produce exponentially streched mesh:

   
il=round(maxh/dz);
ir=nk-il;

zf   =  zeros(nk,1);
dzf  =  zeros(nk,1); 
zh   =  zeros(nk+1,1);

zf(1:il)=dz/2:dz:maxh;
zh(1:il+1)=0:dz:maxh;

gf=stretchconst;    

while true
zh(il+1:end)=zh(il+1)+(dh-zh(il+1))*(exp(gf*(0:1:ir)/(ir))-1)/(exp(gf)-1);

if (zh(il+2)-zh(il+1))<dz
gf=gf-0.01; %make sufficiently small steps to avoid an initial bump in dz
disp(['decreasing stretchconst to:' num2str(gf)])

else
    if (zh(end)-zh(end-1))>3*dz
    disp('WARNING: final grid spacing large! Consider reducing domain height') 
    end
    break
end

end

for i=1:nk
zf(i)=(zh(i)+zh(i+1))/2 ;
dzf(i)=zh(i+1)-zh(i);
end

disp(['growth factor ~' num2str((dzf(il+1)/dzf(il)-1)*100) '-' num2str((dzf(end)/dzf(end-1)-1)*100) '%'])

if ltestplot
figure
subplot(1,2,1)
hold on
plot(zf,'x-')
plot([0 nk],[maxh maxh],'k-')
axis tight
ylabel('z_f [m]')
xlabel('n')
subplot(1,2,2)
hold on
plot(dzf,'d-')
ylabel('dz_f [m]')
xlabel('n')
axis tight
end
