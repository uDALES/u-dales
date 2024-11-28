function geom = createFlatSurface(xsize, ysize, edgelength)
% createFlatSurface    creates flat surface consisting of triangular facets
%
%    geom = createFlatSurface(xsize, ysize, edgelength) returns a geom
%    instance that can be saved to an stl file.
%        xsize:  length of the domain in x-direction
%        ysize:  length of the domain in y-direction
%        edgelength: the length of individual facets. Best taken as xsize
%        (or ysize) divided by an integer number.

TR = generateGround([], xsize, ysize, edgelength);
geom = udgeom.udgeom(TR);