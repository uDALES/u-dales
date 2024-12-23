function TR_out = unitCube()
% create a unit cube with open bottom
P = [; ...
   -0.5000,   -0.5000,   -0.5000; ...
   -0.5000,    0.5000,   -0.5000; ...
    0.5000,    0.5000,   -0.5000; ...
    0.5000,   -0.5000,   -0.5000; ...
    0.5000,   -0.5000,    0.5000; ...
   -0.5000,   -0.5000,    0.5000; ...
   -0.5000,    0.5000,    0.5000; ...
    0.5000,    0.5000,    0.5000; ...
];

T = [ ...

     1,     2,     3; ...
     3,     4,     1; ...
     1,     4,     5; ...
     5,     6,     1; ...
     1,     6,     7; ...
     7,     2,     1; ...
     8,     5,     4; ...
     4,     3,     8; ...
     8,     3,     2; ...
     2,     7,     8; ...
     8,     7,     6; ...
     6,     5,     8; ...
];

TR_out = triangulation(T,P);