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
if ps == 2 
    l = (2 * l * u / (l + 2 * u));
end

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

equations = nnd * dim - length(disp_BC);
K = zeros([equations, equations]);
gK = zeros([nnd * dim, nnd * dim]);
F = zeros([equations, 1]);
gF = zeros([nnd * dim, 1]);
Mass = zeros([nnd * 3, nnd * 3]);
for e=element'
    [k, f, m] = local_values(nodes(e, 1), nodes(e, 2), u, l, bforce, nen);
    gx = ID(e, 1);
    gy = ID(e, 2);
    gs = reshape([gx' ; gy'], [], 2 * nen)';
    es = [(2*e -1)'; (2 * e)']';

    interleave = gs > 0;
    F(gs(interleave)) = F(gs(interleave)) + f(interleave);
    gF(es) = gF(es) + f(interleave);

    for ex = 1:nen
        gx0 = gx(ex) > 0;
        for ex1 = ex:nen
            if gx0 && gx(ex1) > 0
                K(gx(ex), gx(ex1)) = K(gx(ex), gx(ex1)) + k(2 * ex - 1, 2 * ex1 - 1);
                K(gx(ex1), gx(ex)) = K(gx(ex), gx(ex1));
            end
            Mass(3*e(ex)-2:3*e(ex),3*e(ex1)-2:3*e(ex1)) = ...
                Mass(3*e(ex)-2:3*e(ex),3*e(ex1)-2:3*e(ex1)) + ...
                m(3*ex-2:3*ex, 3*ex1-2:3*ex1);
            Mass(3*e(ex1)-2:3*e(ex1),3*e(ex)-2:3*e(ex)) = ...
                Mass(3*e(ex)-2:3*e(ex),3*e(ex1)-2:3*e(ex1));
            gK(2*e(ex) - 1, 2*e(ex1)-1) = gK(2*e(ex) - 1, 2*e(ex1) - 1) + k(2 * ex - 1, 2 * ex1 - 1);
            gK(2*e(ex1), 2*e(ex)-1) = gK(2*e(ex) - 1, 2*e(ex1) - 1);
        end
        % Partially assembled projection matrix
        %P(3*e(ex)-2:3*e(ex), :) = P(3*e(ex)-2:3*e(ex), :) + p(3*ex-2:3*ex, :);
        % Add the force contribution from the constrained x
        if ~gx0
            bc = M(e(ex));
            F(gs(interleave)) = F(gs(interleave)) ... % Previous force contributions
                - k(interleave, ex * 2 - 1) * bc(1); % discount displacement contribution
        end

        for ey = 1:nen
            gy0 = gy(ey) > 0;
            if ex == 1
                for ey1 = ey:nen
                    if gy0 && gy(ey1) > 0
                        K(gy(ey), gy(ey1)) = K(gy(ey), gy(ey1)) + k(2 * ey, 2 * ey1);
                        K(gy(ey1), gy(ey)) = K(gy(ey), gy(ey1));
                    end
                    gK(2*e(ey), 2*e(ey1)) = gK(2*e(ey), 2*e(ey1)) + k(2 * ey, 2 * ey1);
                    gK(2*e(ey1), 2*e(ey)) = gK(2*e(ey), 2*e(ey1));
                end
                if ~gy0
                    % Add the force contribution from the constrained y
                    bc = M(e(ey));
                    F(gs(interleave)) = F(gs(interleave)) ... % Previous force contributions
                        - k(interleave, ey * 2) * bc(2); % discount displacement contribution
                end
            end
            if gx0 && gy0
                K(gx(ex), gy(ey)) = K(gx(ex), gy(ey)) + k(2 * ex, 2 * ey - 1);
                K(gy(ey), gx(ex)) = K(gy(ey), gx(ex)) + k(2 * ex - 1, 2 * ey);
            end
            gK(2*e(ex)-1, 2*e(ey)) = gK(2*e(ex)-1, 2*e(ey)) + k(2 * ex, 2 * ey - 1);
            gK(2*e(ey), 2*e(ex)-1) = gK(2*e(ey), 2*e(ex)-1) + k(2 * ex - 1, 2 * ey);
        end
    end
end

R = chol(K);
d = R\(R'\F);

disp = zeros([nnd, dim]);
disp(ID > 0) = d(ID(ID > 0));
for n = 1:nnd
    for ds = 1:dim
        if ID(n, ds) == 0
            m = M(n);
            disp(n, ds) = m(ds);
        end
    end
end

for e=element'
    P(reshape(((3*e-3) + uint32(1:3))', [], 1)) = P(reshape(((3*e-3) + uint32(1:3))', [], 1)) + ... 
        local_projection(nodes(e,1), nodes(e,2), disp(e, :), u, l, nen)';
end

stresses = Mass\P;
reaction = gK * disp - gF;