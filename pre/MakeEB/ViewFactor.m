function [F12,F21,A1,A2] = ViewFactor(C1,C2,A1,A2,GP,vcorner)
%ViewFactor calculates view factors between two planar surfaces in 3D.
%   [F12,F21,A1,A2] = ViewFactor(C1,C2) returns the view factors between
%   planer surfaces whose vertices are stored in coordinate arrays 
%   C1 and C2.  Cn are of the form [x1,y1,z1; x2,y2,z2; x3,y3,z3;....],
%   where (xn,yn,zn) is the coordinate of the nth vertex.  The vertices 
%   should be in order as  encountered on a trip around the perimeter, and
%   each vertex should be counted only once. The integration is carried 
%   out with 7-point Gauss-Legendre quadrature.

%
% adapted from  Matt Fig, Date: 09/20/2016

if nargin < 3
     GP = 7;
end
[A,W] = Gauss(GP); % Weights and abscissa for Gauss-Legendre quadrature



if vcorner(1)==0
L1 = 4; %size(C1,1);   % Number of vertices.
L2=[4 4 4 4]; 
%L2=[size(C2,1) size(C2,1) size(C2,1) size(C2,1)];
% Close the boundary.
C1(5,:) = C1(1,:);
C2(5,:) = C2(1,:);
c=0;
else
dorder=[1 2 3 4];
L1 = 4; %size(C1,1);   % Number of vertices.
L2=[4 4 4 3]; 
order1=circshift(dorder,4-vcorner(1));
order2=circshift(dorder,4-vcorner(2));
C1=C1(order1,:);
C2=C2(order2,:);
C1(5,:) = C1(1,:);
C2(5,:) = C2(1,:);
Lc=norm(C1(4,:)-C1(1,:));
c=Lc^2*(1.5-log(Lc))*4*sign(vcorner(3));
end

% A1 = rectarea3d(C1);
% A2 = rectarea3d(C2);
% A1 = polyarea3d(C1);
% A2 = polyarea3d(C2);
S=0;

for ii = 1:L1 % Loop over segment pairs.
    P1 = C1(ii,:);
    P2 = C1(ii+1,:);
    
    for jj = 1:L2(ii)       
        SM = 0;
        P3 = C2(jj,:);
        P4 = C2(jj+1,:);
        % Next perform the Gauss-Legendre quadrature
        for kk = 1:GP
            SM = SM + sum(W(kk)*W.*F(A(kk),A,P1,P2,P3,P4));
        end
        
        S = S + SM;
    end
end

%Calculation of the view factors
%cp = cross(C1(2,:) - C1(1,:),C1(3,:) - C1(1,:));
%A1 = polyarea3d(C1,cp/(norm(cp)));
%cp = cross(C2(2,:)-C2(1,:),C2(3,:)-C2(1,:));
%A2 = polyarea3d(C2,cp/(norm(cp)));

F12 = (abs(S+c))/(8*pi*A1);
F21 = (abs(S+c))/(8*pi*A2);


function INT  = F(s,t,P1,P2,P3,P4)
% Integrand for contour integral.  We parameterize the linesegments on 
% s,t=[-1,1] for easy 2D Gaussian quadrature..
% Parametric equations for the first line segment.
T = (t+1)/2;
x1 = P1(1)+(P2(1)-P1(1))*T;
y1 = P1(2)+(P2(2)-P1(2))*T;
z1 = P1(3)+(P2(3)-P1(3))*T;
% Parametric equations for the second line segment.
S = (s+1)/2;
x2 = P3(1)+(P4(1)-P3(1))*S;
y2 = P3(2)+(P4(2)-P3(2))*S;
z2 = P3(3)+(P4(3)-P3(3))*S;
% Distance between seg 1 and seg 2.
R = sqrt((x2-x1).^2 + (y2-y1).^2 + (z2-z1).^2);
% The integrand
INT = log(R)*((P2(1)-P1(1))*(P4(1)-P3(1))+(P2(2)-P1(2))*...
              (P4(2)-P3(2))+(P2(3)-P1(3))*(P4(3)-P3(3)));
             

% function [area] = polyarea3d(V,N)
% % Calculates the area of the polygon given by the matrix of vertices and
% % the normal vector N.  This assumes the vertices are closed, i.e., the
% % last row in V is equal to the first row.  Algorithm developed by 
% % Daniel Sunday and available in c++ here:
% % http://geomalgorithms.com/a01-_area.html
% n = size(V,1)-1;
% % select largest abs coordinate to ignore for projection
% ax = abs(N(1));
% ay = abs(N(2));
% az = abs(N(3));  
% 
% if (ax > ay)&&(ax > az)
%     area = (sum(V(2:n,2) .* (V(3:n+1,3) - V(1:n-1,3))) + ...
%            (V(n+1,2) .* (V(2,3) - V(n,3))))/ (2 * N(1));
% elseif (ay > az)&&(ay > ax)
%     area = (sum((V(2:n,3) .* (V(3:n+1,1) - V(1:n-1,1)))) + ...
%            (V(n+1,3) * (V(2,1) - V(n,1))))/ (2 * N(2));
% else
%     area = (sum((V(2:n,1) .* (V(3:n+1,2) - V(1:n-1,2)))) + ...
%            (V(n+1,1) * (V(2,2) - V(n,2))))/(2 * N(3));
% end

% function [area] = rectarea3d(V)
% % if only rectangles are used
% %area=(0.2+norm(V(1,:)-V(2,:)))*(0.2+norm(V(2,:)-V(3,:)));
% l1=norm(V(1,:)-V(2,:));
% l2=norm(V(2,:)-V(3,:));
% area=l1*l2;



function [A,W] = Gauss(n)
% Generates the abscissa and weights for a Gauss-Legendre quadrature.
% Reference:  Numerical Recipes in Fortran 77, Cornell press.
A = zeros(n,1);  % Preallocations.
W = A;
m = (n+1)/2;
for ii=1:m
    z = cos(pi*(ii-.25)/(n+.5)); % Initial estimate.
    z1 = z+1;
    
    while abs(z-z1)>eps
        p1 = 1;
        p2 = 0;
        
        for jj = 1:n
            p3 = p2;
            p2 = p1;
            p1 = ((2*jj-1)*z*p2-(jj-1)*p3)/jj; % The Legendre polynomial.
        end
        
        pp = n*(z*p1-p2)/(z^2-1);   % The L.P. derivative.
        z1 = z;
        z = z1-p1/pp;
    end
    
    A(ii) = -z;      % Build up the abscissas.
    A(n+1-ii) = z;
    W(ii) = 2/((1-z^2)*(pp^2));  % Build up the weights.
    W(n+1-ii) = W(ii);
end
