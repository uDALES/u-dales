# Fast_Point2TriMesh

![](https://github.com/thor-andreassen/Fast_Point2TriMesh/blob/main/fastPoint2TriMesh_Image.png)

Determine the nearest point between an abitrary point in space and a given triangulated surface
This code takes in a triangulated surface of faces and nodes/vertices such as those from an STL file. As well as, a list of arbitrary points in space
and determines the nearest projection point and the distance to the surface for every point.

The main function is "fastPoint2TriMesh.m". An example of the function being used for a random set of points to a triangulated surface form a mensicus is 
shown in "test_fast_point2trimesh.m"

 This code takes in a triangulated surface of faces and nodes/vertices
     such as those present in STLs. It then takes a set of arbitrary
     points in space and determines the nearest point on the triangulated
     surface to this point, as well as the distance away. This code is
     similar to an existing algoirthm by Daniel Frisch called
     "point2TriMesh", but is written more efficiently and provides
     significant speed increases.
    
     For example, to determine the projection
     of 10,000 points to the included example surface, the time in
     Point2TriMesh is approximately 100 seconds. Using this implementation
     the exact same points can be obtained in approximately 0.3 seconds
     without parallel computing and 0.18 seconds with parallel computing.
     Additionally, the proposed algorithm allows values to be stored
     allowing for subsequent computations to take less time.
    
     The algorithm works by either taking in or calculating the circle
     incenter of every triangle, as well as the normal direction to each
     face or calculating it. Additionally, it either takes in or
     calculates an initial KDTree for the original face incenters. The
     code then determines the nearest average face to every one of the
     inputted test/query points using a KNN Search and the previously
     created KDTree. These points are used to determine the initial
     distance and projection. The code then checks every projection point,
     and first ensures that it is inside the given triangle. If the
     calcualted projection point is outside of the triangle, the code
     determines the minimum distance to each of the face's 3 edges to the
     arbitrary point and then uses the minimum distance and correpsonding
     point as the projection. This is repeated for every point to get the
     distance, and projection point for every query point. The face
     normals are used to determine a positive or negative sign depending
     on whether the point is nearest the positive normal direction
     (outside) or the negative normal direction (inside)
    
     INPUTS:
     inputs.faces = (N x 3) the faces of the original triangulated surface (Required)
     inputs.nodes = (M x 3) the nodes or vertices of the original triangulated surface (Required)
     inputs.face_mean_nodes = (N x 3) the location of the incenter of each
                         triangle. This can either be calculated using "getTriInCenter.m" ahead
                         of time, or calculated herein.
     inputs.face_normals = (N x 3) the unit normal direction to each face.
                         This can either be calculated using "getFaceCenterAndNormals.m" ahead
                         of time, or calculated herein.
     inputs.tree_model = (model) the KDTree trained to each triangle
                         incenter using "KDTreeSearcher.m" This can either be calculated using "KDTreeSearcher.m"
                         ahead of time, or calculated herein.
     pts = (Q x 3) Arbitrary points that we are trying to project onto the
                         surface (Required)
     use_parallel = (0 or 1) A binary that determines whether to use (1) the parallel computing
                          or not use (0) the parallel computing. NOTE this requires the parallel computing toolbox(Required)
    
     OUTPUTS:
     distances= (Q x 1) The signed distance for each arbitary point (pts)
                         to the nearest point on the triangulated surface. Positive means the
                         point is off of the positive face normal, where as negative means the
                         point is nearest the opposite direction to the
                         face normal.
     project_pts= (Q x 3) The location of the projected point for each arbitary point (pts)
                         to the nearest point on the triangulated surface.
     outside= (Q x 1) A boolean array representing if the given point is on the positive direction
                         (1) or the negative direction (0) of the corresponding surface
    
    
     Written by Thor Andreassen
     University of Denver
     5/9/2023
