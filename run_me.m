%% This script loads and solves a static finite element problem.
% The core of this program is written in Rust for multithreading and lapack
% integration. We need to load these libraries first. I don't think the
% loading of these libraries should count towards my time?

%addpath('examples', 'matlab');
%library_path = 'target\debug\ferust';
%header_path = 'src\ferust.h';
%solver = Solver(library_path, header_path);

%% Solve for displacements
% First we need to parse our input file
filename = "Biaxial_Q4_2x2.txt";
[nodes, element, elemType, nel, nen, nIntPts, nnd, ps, nu, E, ...
    Force_Node, bforce, disp_BC] = Read_input(filename);

u = E / (2*(1 + nu));
l = E * nu / ((1 + nu) * (1 - 2 * nu));
l = (2 * l * u / (l + 2 * u));

dims = cast(2, 'uint16');

bc_nodes = cast(disp_BC(:, 1), 'uint16');
M = containers.Map(bc_nodes, ...
        mat2cell( ...
            NaN([length(bc_nodes) dims]), ...
            ones([1 length(bc_nodes)])));

for row = disp_BC'
    m = M(row(1));
    m(row(2)) = row(3);
    M(row(1)) = m;
end

ID = zeros([length(nodes) dims], 'uint16');
id = 1;
for node = 1:length(nodes)
    if ~M.isKey(node)
        ID(node, 1:dims) = id:id + dims - 1;
        id = id + dims;
        continue;
    end
    row = M(node);
    for dim = 1:dims
        if isnan(row(dim))
            ID(node, dim) = id;
            id = id + 1;
        end
    end
end

equations = nnd * 2 - length(disp_BC);
K = zeros([equations, equations]);

for e=element
    k = local_stiffness(nodes(e, 1), nodes(e, 2), u, l, nen);
    gx = ID(e, 1);
    gy = ID(e, 2);
    for ex = 1:nen
        gx0 = gx(ex) > 0;
        if gx0
            K(gx(ex), gx(ex)) = K(gx(ex), gx(ex)) + k(2 * ex - 1, 2 * ex - 1);
        end
        for ey = 1:nen
            if gy(ey) > 0
                if ex == 1
                    K(gy(ey), gy(ey)) = K(gy(ey), gy(ey)) + k(2 * ey, 2 * ey);
                end
                if gx0
                    K(gx(ex), gy(ey)) = K(gx(ex), gy(ey)) + k(2 * ex, 2 * ey - 1);
                    K(gy(ey), gx(ex)) = K(gy(ey), gx(ex)) + k(2 * ex - 1, 2 * ey);
                    %K(gx(ex), gy(ey) - 1) = K(gx(ex), gy(ey)) + k(2 * ex, 2 * ey - 1);
                    %K(gx(ex) - 1, gy(ey)) = K(gy(ey), gx(ex)) + k(2 * ex - 1, 2 * ey);
                end
            end
        end
    end
end

%R = chol(K);
%d = R\(R'\F);

% Now we can run our solver
%handle = @() solver.solve(node, element, elemType, nel, nen, nIntPts, ...
%    nnd, ps, nu, E, Force_Node, bforce, disp_BC);
%timeit(handle)
%disp("If it ran, we're probably good :)");