function [face_mean_nodes,face_normals]=getFaceCenterAndNormals(faces,nodes)
% determines the incenter of a face and the normal direction to the face    
face_mean_nodes=zeros(size(faces));
    face_normals=zeros(size(faces));
    for count_face=1:size(faces,1)
        current_nodes=nodes(faces(count_face,:),:);
        face_mean_nodes(count_face,:)=getTriInCenter(current_nodes);

        vec1=current_nodes(2,:)-current_nodes(1,:);
        vec2=current_nodes(3,:)-current_nodes(2,:);
        vec3=cross(vec1,vec2);
        vec3=vec3/norm(vec3);

        face_normals(count_face,:)=vec3;
    end
end