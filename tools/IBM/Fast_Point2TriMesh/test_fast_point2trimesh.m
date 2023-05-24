%% clear
clear
close all
clc

%% Load geometry
[men_geom]=stlread('meniscus_2.stl');
faces=men_geom.ConnectivityList;
nodes=men_geom.Points;
%% get face means and normals (STEP 1)
% do this once and store
[face_mean_nodes,face_normals]=getFaceCenterAndNormals(faces,nodes);

%% create KD tree (STEP 2)
% do this once and store
tree_model=KDTreeSearcher(face_mean_nodes);

%% create test points (these are the query or test points)
bb.min=min(nodes)-10;
bb.max=max(nodes)+10;
num_pts=100;
pts=rand(num_pts,3).*(bb.max-bb.min)+bb.min;


%% determine distance and projection point (STEP 3)
tic
inputs.faces=faces;
inputs.nodes=nodes;
inputs.face_mean_nodes=face_mean_nodes;
inputs.face_normals=face_normals;
inputs.tree_model=tree_model;
[distances,project_pts,outside]=fastPoint2TriMesh(inputs,pts,1);
toc

%% plot test data
figure()
patch('Faces',faces,'Vertices',nodes,'FaceColor','r','FaceAlpha',.75,'EdgeAlpha',.2);
hold on
nearest_direction=project_pts-pts;
scatter3(pts(:,1),pts(:,2),pts(:,3),'bx');
quiver3(pts(:,1),pts(:,2),pts(:,3),nearest_direction(:,1),nearest_direction(:,2),nearest_direction(:,3),'k','AutoScale','off');


legend({'Object','Points','Direction to Surface'})
view([1,1,1]);