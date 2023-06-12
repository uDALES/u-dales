% clear all; clc;

function flag = point_line_segment_intersect(O,d,V1,V2)

% O = [5 5 -5]';
% d = -5*[0 0 1]';
% V1 = [5 5 0]';
% V2 = [5 5 5]';
% size_V1 = size(V1);
% size_V2 = size(V2);
% size_d = size(d);
% size_O = size(O);
% if(size_V1(1)~=3 || size_V1(2)~=1 || size_V2(1)~=3 || size_V2(2)~=1 || size_d(1)~=3 || size_d(2)~=1 || size_O(1)~=3 || size_O(2)~=1)
%     error('All inputs must be in 3x1 format.')
% end

edge = V2-V1; 
rhs = O-V1;

deno(1) = -d(1)*edge(2) + d(2)*edge(1);
deno(2) = -d(1)*edge(3) + d(3)*edge(1);
deno(3) = -d(2)*edge(3) + d(3)*edge(2);

[max_deno,max_i] = max(abs(deno));
max_deno = max_deno*sign(deno(max_i));


    switch max_i
        case 1 
            neu_1 = -edge(1)*rhs(2)+edge(2)*rhs(1);
            neu_2 = rhs(1)*d(2)-rhs(2)*d(1);
        case 2
            neu_1 = -edge(1)*rhs(3)+edge(3)*rhs(1);
            neu_2 = rhs(1)*d(3)-rhs(3)*d(1);
        case 3
            neu_1 = -edge(2)*rhs(3)+edge(3)*rhs(2);
            neu_2 = rhs(2)*d(3)-rhs(3)*d(2);
        otherwise
            error('This function is valid only for 3D inputs.')
    end

    if (max_deno~=0)
        
        t_1 = neu_1/max_deno;
        t_2 = neu_2/max_deno;
        if(t_1>=0 && t_2>=0 && t_2<=1)
            flag = true;
        else
            flag = false;
        end

    elseif (max_deno==0 && neu_1==0 && neu_2==0)

        if (edge(1)==0 && edge(2)==0 && edge(3)==0)
            error('Points A and B mast have different coordinates')
        elseif(d(1)==0 && d(2)==0 && d(3)==0)
            error('Magnitude of the direction vector cannot be zero')
        elseif (norm(cross(rhs,edge)) == 0)
            [max_edge,max_edge_i] = max(abs(edge));
            if(dot(d,edge)>0 && rhs(max_edge_i)/edge(max_edge_i) <= 1 )
                flag = true;
            elseif (dot(d,edge)<0 && rhs(max_edge_i)/edge(max_edge_i) >= 0 )
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
    
% figure
% hold on
% plot3([V1(1) V2(1)],[V1(2) V2(2)],[V1(3) V2(3)],'-o')
% quiver3(O(1),O(2),O(3),d(1),d(2),d(3))
% plottools