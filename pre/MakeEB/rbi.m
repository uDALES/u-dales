function [flag ,tmin] = rbi(origin, direction, vec)
%  Ray/box intersection using the Smits' algorithm
%
% Input:
%    origin.
%    direction.
%    vec = [xl,xu,yl,yu,zl,zu]  (=Block)
% Output:
%    flag: (0) Reject, (1) Intersect.
%    tmin: distance from the ray origin.

vmin=[vec(1) vec(3) vec(5)];
vmax=[vec(2) vec(4) vec(6)];


if ((direction(1) == 0) && (((vmin(1) - origin(1))==0) || ((vmax(1) - origin(1))==0))) %comparing floats!  %0/0 is undefined whereas 1/0 = inf, -1/0=-inf works
    %this will only be the case if the ray is tangential -> no intersection
    %(or infinite intersections...)
    tmin=0;
    tmax=inf;
elseif (direction(1) >= 0)
    tmin = (vmin(1) - origin(1)) / direction(1);
    tmax = (vmax(1) - origin(1)) / direction(1);
else
    tmin = (vmax(1) - origin(1)) / direction(1);
    tmax = (vmin(1) - origin(1)) / direction(1);
end

if ((direction(2) == 0) && (((vmin(2) - origin(2))==0) || ((vmax(2) - origin(2))==0))) %comparing floats!
    tymin=0;
    tymax=inf;
elseif (direction(2) >= 0)
    tymin = (vmin(2) - origin(2)) / direction(2);
    tymax = (vmax(2) - origin(2)) / direction(2);
else
    tymin = (vmax(2) - origin(2)) / direction(2);
    tymax = (vmin(2) - origin(2)) / direction(2);
end

if ( (tmin > tymax) || (tymin > tmax) )
    flag = 0;
    tmin = NaN;  %originally -1
    return;
end

if (tymin > tmin)
    tmin = tymin;
end

if (tymax < tmax)
    tmax = tymax;
end

if ((direction(3) == 0) && (((vmin(3) - origin(3))==0) || ((vmax(3) - origin(3))==0))) %comparing floats!
    tzmin=0;
    tzmax=inf;
elseif (direction(3) >= 0)
    tzmin = (vmin(3) - origin(3)) / direction(3);
    tzmax = (vmax(3) - origin(3)) / direction(3);
else
    tzmin = (vmax(3) - origin(3)) / direction(3);
    tzmax = (vmin(3) - origin(3)) / direction(3);
end


if ((tmin > tzmax) || (tzmin > tmax))
    flag = 0;
    tmin = -1;
    return;
end

if (tzmin > tmin)
    tmin = tzmin;
end

if (tzmax < tmax)
    tmax = tzmax;
end

% if( (tmin < t1) && (tmax > t0) )
flag = 1;
% else
%    flag = 0;
%    tmin = -1;
% end;
end
