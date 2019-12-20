blockfile='blocks4.inp';
ncolumns=6;  %13 unless blocks.inp has changed its format

 
ni=400
nj=320
nk=160;
xsize=1200;
ysize=960;
zsize=480

dx=3;
dy=3;
dz=3;


%read blocks
fid = fopen(blockfile)
data = textscan(fid,['%f'],'HeaderLines',3); %include bottom block as header
data2=cell2mat(data)
l=length(data2)
data=reshape(data2,ncolumns,l/ncolumns)
fid=fclose(fid)
 
%%
%3d-patch indeces
figure
for i=1:l/ncolumns
vert = [data(1,i) data(3,i) data(5,i); data(2,i)+0.99 data(3,i) data(5,i); data(2,i)+0.99 data(4,i)+0.99 data(5,i); data(1,i) data(4,i)+0.99 data(5,i); ...
        data(1,i) data(3,i) data(6,i)+0.99; data(2,i)+0.99 data(3,i) data(6,i)+0.99; data(2,i)+0.99 data(4,i)+0.99 data(6,i)+0.99; data(1,i) data(4,i)+0.99 data(6,i)+0.99];
fac = [1 2 3 4; ...
       2 6 7 3; ...
       4 3 7 8; ...
       1 5 8 4; ...
       1 2 6 5; ...
       5 6 7 8];

p=patch('Faces',fac,'Vertices',vert,'FaceColor',[0.4 0.4 0.4],'FaceLighting','none'); 
  
end
axis equal
xlim([0 ni])
ylim([0 nj])
zlim([0 nk])
 %%
%3d-patch units
figure
for i=1:l/ncolumns
vert = [data(1,i)*dx data(3,i)*dy data(5,i)*dz; data(2,i)*dx+dx-eps(1) data(3,i)*dy data(5,i)*dz; data(2,i)*dx+dx-eps(1) data(4,i)*dy+dy-eps(1) data(5,i)*dz; data(1,i)*dx data(4,i)*dy+dy-eps(1) data(5,i)*dz; ...
        data(1,i)*dx data(3,i)*dy data(6,i)*dz+dz-eps(1); data(2,i)*dx+dx-eps(1) data(3,i)*dy data(6,i)*dz+dz-eps(1); data(2,i)*dx+dx-eps(1) data(4,i)*dy+dy-eps(1) data(6,i)*dz+dz-eps(1); data(1,i)*dx data(4,i)*dy+dy-eps(1) data(6,i)*dz+dz-eps(1)];
fac = [1 2 3 4; ...
       2 6 7 3; ...
       4 3 7 8; ...
       1 5 8 4; ...
       1 2 6 5; ...
       5 6 7 8];

p=patch('Faces',fac,'Vertices',vert,'FaceColor',[0.4 0.4 0.4],'FaceLighting','none'); 
  
end
xlim([0 xsize])
ylim([0 ysize])
zlim([0 zsize])

%%
blub=zeros(ni,nj,nk);
for i=1:l/ncolumns
blub((data(1,i):data(2,i)),(data(3,i):data(4,i)),(1:data(6,i)+1))=1;
end

[X Y Z] = meshgrid(1:ni,1:nj,1:nk) ;
moep=permute(blub,[2 1 3]);

p = patch(isosurface(X,Y,Z,moep,1));