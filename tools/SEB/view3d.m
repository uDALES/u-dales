function vf = view3d(view3d_exe, infile, outfile)
% Wrapper for view3d (view factor calculator).
% view3d_exe is the path to the view3d executable
% infile is the path to the input file
% outfile is the path to the output file
% vf: view factors.

view3d_execution_command = [view3d_exe ' ' infile ' ' outfile];
system(view3d_execution_command);
vf = dlmread(outfile, ' ', 2, 0);
vf(end,:) = [];
end