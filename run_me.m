%% This script loads and solves a static finite element problem.
% The core of this program is written in Rust for multithreading and lapack
% integration. We need to load these libraries first. I don't think the
% loading of these libraries should count towards my time?

addpath('examples', 'matlab');
library_path = 'target\debug\ferust';
header_path = 'src\ferust.h';
solver = Solver(library_path, header_path);

%% Solve for displacements
% First we need to parse our input file
filename = "Biaxial_Q4_2x2.txt";
[node, element, elemType, nel, nen, nIntPts, nnd, ps, nu, E, ...
    Force_Node, bforce, disp_BC] = Read_input(filename);

% Now we can run our solver
handle = @() solver.solve(node, element, elemType, nel, nen, nIntPts, ...
    nnd, ps, nu, E, Force_Node, bforce, disp_BC);
timeit(handle)
disp("If it ran, we're probably good :)");