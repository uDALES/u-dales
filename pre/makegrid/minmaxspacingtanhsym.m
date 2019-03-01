lmakegrid=1;
llinbot=1;
n=480;

%desired spacing
spacing1=2.534
spacing2=2
%available points
ngridpoints=28

hmin=spacing1*ngridpoints; %=>delta=0.0000001 would work


%% to determine stretching coefficients delta  (=d)

gp1=linspace(0,1,ngridpoints);
ming=zeros(1000,1);
maxg=zeros(1000,1);
di=linspace(0.01,2,1000);
for i=1:1000
d=di(i);
g=1/2*(1+tanh(d.*(gp1/1-1/2))/tanh(d/2))*1;    
gg=gradient(g);
ming(i)=min(gg);
maxg(i)=max(gg);
end
list=maxg./ming;

%% 
[index]=find(list<=spacing2/spacing1,1,'last') %choose delta corresponding to desired spacing 
d=di(index);

h=spacing2/maxg(index); %determine maximum height with given delta and spacings
gp=linspace(0,h,ngridpoints); %create grid
g=1/2*(1+tanh(d.*(gp/h-1/2))/tanh(d/2))*h;  
gg=gradient(g);

%gain due to tanh-mesh instead of linear
sprintf('gain due to tanh-mesh instead of linear: %2.1f%s.',g(end)-hmin,'m')


%% add linear spacing at bottom etc..
if lmakegrid
grid=zeros(n,1);
if llinbot
   ne=(n-ngridpoints+1);
   np=1:ne;
   grid(np)=spacing1/2+(np-1)*spacing1;
end
   grid(ne:end)=grid(ne)+g;
end
format long g
round(grid,4)

figure
hold on
for i=1:2:n-1
y1=0;
y2=1;
x1=grid(i);
x2=grid(i+1);
x = [x1, x2, x2, x1, x1];
y = [y1, y1, y2, y2, y1];
plot(x, y, 'k-', 'LineWidth', 0.5);

end
xlim([0 grid(end)])
plot(grid,gradient(grid))
hold off

