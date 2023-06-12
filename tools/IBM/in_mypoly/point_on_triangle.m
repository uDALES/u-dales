% function flag = point_on_triangle(Origin,x_tri,y_tri,z_tri,incenters,faceNormals)
function flag = point_on_triangle(Origin,edge1,edge2,rhs,incenters,faceNormals,tol)
% clear all; clc;
% 
% Origin = [5 4 5]';
% x_tri = [5 5 5]';
% y_tri = [0 5 5]';
% z_tri = [0 0 5];
% incenters = [5 3 2]';
% faceNormals = [1 0 0]';

V = Origin-incenters;
for i=1:3
    if (abs(V(i)) < tol)
        V(i) = 0;
    end
end

if(abs(dot(V,faceNormals))<tol)
% if(dot(V,faceNormals)==0)

    % edge1 = [(x_tri(2)-x_tri(1)) (y_tri(2)-y_tri(1)) (z_tri(2)-z_tri(1))];
    % edge2 = [(x_tri(3)-x_tri(1)) (y_tri(3)-y_tri(1)) (z_tri(3)-z_tri(1))];

    % rhs = Origin - [x_tri(1) y_tri(1) z_tri(1)]';

    deno(1) = edge1(1)*edge2(2) - edge1(2)*edge2(1);
    deno(2) = edge1(1)*edge2(3) - edge1(3)*edge2(1);
    deno(3) = edge1(2)*edge2(3) - edge1(3)*edge2(2);

    [max_deno,max_i] = max(abs(deno));
    max_deno = max_deno*sign(deno(max_i));

    switch max_i
            case 1 
                neu_1 = -edge2(1)*rhs(2)+edge2(2)*rhs(1);
                neu_2 = -rhs(1)*edge1(2) + rhs(2)*edge1(1);
            case 2
                neu_1 = -edge2(1)*rhs(3)+edge2(3)*rhs(1);
                neu_2 = -rhs(1)*edge1(3)+rhs(3)*edge1(1);
            case 3
                neu_1 = -edge2(2)*rhs(3)+edge2(3)*rhs(2);
                neu_2 = -rhs(2)*edge1(3)+rhs(3)*edge1(2);
            otherwise
                error('This function is valid only for 3D inputs.')
    end

    if (max_deno~=0)
        a = neu_1/max_deno;
        b = neu_2/max_deno;
        if(a>=0 && b>=0 && a+b<=1)
            flag = true;
        else
            flag = false;
        end

    else
        flag = false;
    end
else
    flag = false;
end