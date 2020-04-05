classdef Solver < handle
    properties (SetAccess = private)
        nodes         % [n x dim] List of nodes taken from the input file
        elements      % [elements x elements per node] Elements taken from the input file
        integration   % [1] Order of Gaussian integration
        bforce        % [1] body force in the y direction
        mu            % [1] Calculated shear modulus
        lame          % [1] Calculated lame constant
        K             % [dof x dof] System Stiffness Matrix
        gK            % [dim * nodes x dim * nodes] Global Stiffness Matrix
        force         % [dof x 1] System Force Vector
        bForce        % [dim * nodes x 1] Global body force vector
        Mass          % [stress_components * nodes x stress_components * nodes] Global Mass Matrix
        projection    % [stress_components * nodes x 1] Global projection vector
        displacements % [nodes x dim] Calculated displacement of solver.nodes
        stresses      % [stress_components * nodes x 1] Calculated stresses at nodes
        reaction      % [nodes x dim] Calculated reaction forces at nodes.
        node_reaction % Calculated reaction forces at constrained nodes.
        reaction_sum  % [dim x 1] Of the sum of the reaction forces in X and Y
        local_k       % Debug value to see what the local k is
    end
    properties (Constant = true)
        dims = cast(2, 'uint16');
        stress_components = cast(3, 'uint16');
    end
    methods
        function solver = Solver(filename, overrides)
            %Solver: Class to solve a finite element problem, given an input
            % file.

            % First we need to parse our input file
            [nodes, element, ~, ~, nen, integration, nnd, ps, nu, E, ...
                Force_Node, bforce, disp_BC] = Read_input(filename);
            % We allow for programmatic overrides from the input file.
            if nargin > 1
                if overrides.isKey('integration')
                    integration = overrides('integration');
                end
                if overrides.isKey('nu')
                    nu = overrides('nu');
                end
                if overrides.isKey('E')
                    E = overrides('E');
                end
                if overrides.isKey('ps')
                    ps = overrides('ps');
                end
                if overrides.isKey('bforce')
                    bforce = overrides('bforce');
                end
                if overrides.isKey('fixelements')
                    element = solver.fix_elements(nodes, element, nen);
                end
            end

            % We set our plane strain condition
            u = E / (2*(1 + nu));
            l = E * nu / ((1 + nu) * (1 - 2 * nu));
            if ps == 2
                l = (2 * l * u / (l + 2 * u));
            end

            % Bind provided to instance
            solver.nodes = nodes;
            solver.elements = element;
            solver.integration = integration;
            solver.bforce = bforce;
            solver.mu = u;
            solver.lame = l;

            % Initialize our data with 0s
            equations = nnd * solver.dims - length(disp_BC);
            solver.K = zeros([equations, equations]);
            solver.gK = zeros([nnd * solver.dims, nnd * solver.dims]);
            solver.force = zeros([equations, 1]);
            solver.bForce = zeros([nnd * solver.dims, 1]);
            solver.Mass = zeros([nnd * solver.stress_components, ...
                nnd * solver.stress_components]);
            solver.projection = zeros([nnd * solver.stress_components, 1]);
            solver.displacements = zeros([nnd, solver.dims]);

            % Build useful data structures. Note, we opt not to use a LM
            % array, but make extensive usage of the ID array, in addition
            % to a map M for boundary condition lookup.
            M = solver.build_map(solver.dims, disp_BC);
            ID = solver.build_id(M, solver.dims, nnd);

            solver.crunch_local_data(ID, M, nen);

            % With a set K matrix and forces, we can solve for our
            % displacements. Since K is garunteed to be positive definite
            % we solve using Choleksy as a speed up method. There are
            % method for banded solves,s including from our own department!
            % See Hang Lui, R. Mitall 2013, IEEE, GPU-accelerated
            % scalable solver for banded linear systems
            % This is also a great sanity check, since this line fails if K
            % is not positive def.
            %R = chol(solver.K);
            %d = R\(R'\solver.force);
            d = solver.K\solver.force;

            % We now need to massage our local d to a global context.
            solver.contextulize_d(d, ID, M);

            % For stresses, we now need a projection vector, dependent on
            % the displacements just calculated.
            for e=element'
                % Add up the projection contribution from each element
                solver.projection(reshape(((3*e-3) + uint32(1:3))', [], 1)) = ...
                    solver.projection(reshape(((3*e-3) + uint32(1:3))', [], 1)) + ...
                    local_projection(nodes(e,1), nodes(e,2), solver.displacements(e, :), ...
                    u, l, nen, solver.integration)';
            end

            % We can finally solve for our stresses
            %R = chol(solver.Mass);
            %solver.stresses = R\(R'\solver.projection);
            solver.stresses = solver.Mass\solver.projection;

            % And thus our reaction forces.
            solver.reaction = solver.gK * reshape(solver.displacements', [], 1) ...
                - solver.bForce;
            % Including the reactions specifically asked for.
            solver.node_reaction = [solver.reaction(2 * Force_Node - 1) ...
                solver.reaction(2 * Force_Node)];
            solver.reaction_sum = [sum(solver.node_reaction(:, 1)) ...
                sum(solver.node_reaction(:, 1))]';
        end
    end
    methods(Access = private, Hidden = true, Sealed = true)
        function crunch_local_data(solver, ID, M, nen)
            %crunch_local_data: Function to determine local stiffness, force,
            % and mass contributions and assemble them.
            for e=solver.elements'
                % Compute Stifffness, Local Force and local Mass in one
                % Fell swoop.
                [k, f, m] = local_values(solver.nodes(e, 1), ...
                    solver.nodes(e, 2), ...
                    solver.mu, solver.lame, ...
                    solver.bforce, nen, ...
                    solver.integration);
                % Provide indices for the solution K matrix.
                gx = ID(e, 1);
                gy = ID(e, 2);
                gs = reshape([gx' ; gy'], [], 2 * nen)';
                % Provide indices for the global K matrix.
                es = reshape([(2 * e -1)'; (2 * e)'], [], 2 * nen)';

                % Aggregate contribution to the solvable equations.
                interleave = gs > 0;
                solver.force(gs(interleave)) = ...
                    solver.force(gs(interleave)) + f(interleave);
                % Record the global body force felt.
                solver.bForce(es) = solver.bForce(es) + f;

                % We examine every x vs x, y vs y and x vs y of the element
                % and record the pertinent data.
                for ex = 1:nen
                    % In the case gx > 0 we have a boundary condition. As such
                    % this displacement is not solvable.
                    gx0 = gx(ex) > 0;
                    % x vs x
                    for ex1 = ex:nen
                        if gx0 && gx(ex1) > 0
                            solver.K(gx(ex), gx(ex1)) = ...
                                solver.K(gx(ex), gx(ex1)) ...
                                + k(2 * ex - 1, 2 * ex1 - 1);
                            % Symmetry of K lets us just flip the result.
                            solver.K(gx(ex1), gx(ex)) = ...
                                solver.K(gx(ex), gx(ex1));
                        end
                        % Compute Mass x vs x contribution
                        solver.Mass(3*e(ex)-2:3*e(ex),3*e(ex1)-2:3*e(ex1)) = ...
                            solver.Mass(3*e(ex)-2:3*e(ex),3*e(ex1)-2:3*e(ex1)) + ...
                            m(3*ex-2:3*ex, 3*ex1-2:3*ex1);
                        % Mass is also symmetric.
                        solver.Mass(3*e(ex1)-2:3*e(ex1),3*e(ex)-2:3*e(ex)) = ...
                            solver.Mass(3*e(ex)-2:3*e(ex),3*e(ex1)-2:3*e(ex1));

                        % Record global stiffness contribution, regardless of
                        % BC. Since this will be used for our stress
                        % calculations.
                        solver.gK(2*e(ex) - 1, 2*e(ex1) - 1) = ...
                            solver.gK(2*e(ex) - 1, 2*e(ex1) - 1) + k(2 * ex - 1, 2 * ex1 - 1);
                        solver.gK(2*e(ex1)- 1, 2*e(ex) - 1) = solver.gK(2*e(ex) - 1, 2*e(ex1) - 1);
                    end
                    % Add the force contribution from the constrained x
                    if ~gx0
                        bc = M(e(ex));
                        solver.force(gs(interleave)) = ...
                            solver.force(gs(interleave)) ... % Previous force contributions
                            - k(interleave, ex * 2 - 1) * bc(1); % discount displacement contribution
                    end

                    % Inner loop to compute x vs y and y vs y contributions
                    for ey = 1:nen
                        gy0 = gy(ey) > 0;
                        % Check to make sure we only capture y v y once.
                        if ex == 1
                            for ey1 = ey:nen
                                if gy0 && gy(ey1) > 0
                                    solver.K(gy(ey), gy(ey1)) = ...
                                        solver.K(gy(ey), gy(ey1)) + k(2 * ey, 2 * ey1);
                                    solver.K(gy(ey1), gy(ey)) = ...
                                        solver.K(gy(ey), gy(ey1));
                                end
                                solver.gK(2*e(ey), 2*e(ey1)) = ...
                                    solver.gK(2*e(ey), 2*e(ey1)) + k(2 * ey, 2 * ey1);
                                solver.gK(2*e(ey1), 2*e(ey)) = ...
                                    solver.gK(2*e(ey), 2*e(ey1));
                            end
                            if ~gy0
                                % Add the force contribution from the constrained y
                                bc = M(e(ey));
                                solver.force(gs(interleave)) = ...
                                    solver.force(gs(interleave)) ... % Previous force contributions
                                    - k(interleave, ey * 2) * bc(2); % discount displacement contribution
                            end
                        end
                        % Capture cross contributions, x vs. y and y vs. x
                        if gx0 && gy0
                            solver.K(gx(ex), gy(ey)) = ...
                                solver.K(gx(ex), gy(ey)) + k(2 * ex - 1, 2 * ey);
                            solver.K(gy(ey), gx(ex)) = ...
                                solver.K(gx(ex), gy(ey));
                        end
                        solver.gK(2*e(ex)-1, 2*e(ey)) = ...
                            solver.gK(2*e(ex)-1, 2*e(ey)) + k(2 * ex - 1, 2 * ey);
                        solver.gK(2*e(ey), 2*e(ex)-1) = ...
                            solver.gK(2*e(ex)-1, 2*e(ey));
                    end
                end
            end
            solver.local_k = k;
        end
        function contextulize_d(solver, d, ID, M)
            solver.displacements(ID > 0) = d(ID(ID > 0));
            for n = 1:length(solver.nodes)
                for ds = 1:solver.dims
                    if ID(n, ds) == 0
                        m = M(n);
                        solver.displacements(n, ds) = m(ds);
                    end
                end
            end
        end
    end
    methods(Static, Access = private, Hidden = true, Sealed = true)
        function elements = fix_elements(nodes, old_elements, nen)
            %fix_elements: Ensure that convention for node numbering is
            % respected. I don't think this is violated in our examples, but
            % something Adyota said left me paranoid.
            elements = old_elements * 0;
            index = 1;
            for e = old_elements'
                c = convhull(nodes(e, 1), nodes(e, 2));
                [~, new_order] = ismember(...
                    polyshape(nodes(e(c), :), ...
                    'SolidBoundaryOrientation', 'ccw', ...
                    'Simplify', false).Vertices, ...
                    nodes(e, :),'rows');
                if nen == 4
                    elements(index, :) = e(new_order);
                else
                    elements(index, 1:4) = e(new_order(1:2:8));
                    elements(index, 5:8) = e(new_order(2:2:8));
                    if nen == 9
                        mask = true([1 9]);
                        mask(new_order) = 0;
                        elements(index, 9) = e(mask);
                    end
                end
                index = index + 1;
            end
        end
        function M = build_map(dims, disp_BC)
            %build_map: Creates hashmap for boundary condition lookup.
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
        end
        function ID = build_id(M, dims, nnd)
            %build_id: Creates the ID matrix for a given
            % system.
            ID = zeros([nnd dims], 'uint16');
            id = 1;
            for node = 1:nnd
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
        end
    end
    methods(Access = public, Sealed = true)
        function stress_fn = get_interpolate_stress_fn(solver, element)
            %get_interpolate_stress_fn: returns a function that can calculate
            % the stress at a given point for an element.
            e = solver.elements(element, :)';
            xx = local_interpolation(solver.nodes(e, 1), ...
                solver.nodes(e, 2), ...
                solver.stresses(3 * e - 2), length(e));
            yy = local_interpolation(solver.nodes(e, 1), ...
                solver.nodes(e, 2), ...
                solver.stresses(3 * e - 1), length(e));
            xy = local_interpolation(solver.nodes(e, 1), ...
                solver.nodes(e, 2), ...
                solver.stresses(3 * e), length(e));
            stress_fn = @(y, x) [xx(x ,y) yy(x ,y) xy(x ,y)];
        end
        function stress = interpolate_stress(solver, x, y)
            %interpolate_stress: Finds the element that contains x, y and
            % calculates the stresses at that point in the element.
            i = 1;
            for e = solver.elements'
                c = convhull(solver.nodes(e, 1), solver.nodes(e, 2));
                if inpolygon(x, y, solver.nodes(e(c), 1), solver.nodes(e(c), 2))
                    stress_fn = solver.get_interpolate_stress_fn(i);
                    stress = stress_fn(x, y);
                    return;
                end
                i = i + 1;
            end
            stress = NaN([1 solver.stress_components]);
        end
        function displacement_fn = get_interpolate_displacement_fn(solver, element)
            %get_interpolate_displacement_fn: returns a function that can calculate
            % the displacement at a given point for an element.
            e = solver.elements(element, :)';
            X = local_interpolation(solver.nodes(e, 1), ...
                solver.nodes(e, 2), ...
                solver.displacements(e, 1), length(e));
            Y = local_interpolation(solver.nodes(e, 1), ...
                solver.nodes(e, 2), ...
                solver.displacements(e, 2), length(e));
            displacement_fn = @(y, x) [X(x ,y) Y(x ,y)];
        end
        function displacement = interpolate_displacement(solver, x, y)
            %interpolate_displacement: Finds the element that contains x, y and
            % calculates the displacement at that point in the element.
            i = 1;
            for e = solver.elements'
                c = convhull(solver.nodes(e, 1), solver.nodes(e, 2));
                if inpolygon(x, y, solver.nodes(e(c), 1), solver.nodes(e(c), 2))
                    displacement_fn = solver.get_interpolate_displacement_fn(i);
                    displacement = displacement_fn(y, x);
                    return;
                end
                i = i + 1;
            end
            displacement = NaN([1 solver.dims]);
        end
        function [varargout] = contour_stress(solver, resolution)
            %contour_stress: Creates a mesh grid of stress results plotted at
            % the respective displacments. Called with no return, this will
            % plot.
            if nargin == 1
                resolution = 10;
            end
            dt = 1/resolution;
            resolution = uint32(resolution);
            rx = max(solver.nodes(:, 1)) * resolution + 1;
            ry = max(solver.nodes(:, 2)) * resolution + 1;
            stress = zeros([rx, ry, 3]);
            displacement = zeros([rx, ry, 2]);
            for i = 1:rx
                dx = double(i - 1) * dt;
                for j = 1:ry
                    dy = double(j - 1) * dt;
                    stress(i, j, :) = solver.interpolate_stress(dx, dy);
                    displacement(i,j,:) = [dx dy] + solver.interpolate_displacement(dx, dy);
                end
            end
            if nargout == 0
              fn = @(i) contourf(displacement(:, :, 1), displacement(:, :, 2), stress(:, :, i));
              subplot(3,1,1);
              fn(1);
              subplot(3,1,2);
              fn(2);
              subplot(3,1,3);
              fn(2);
            else
              varargout{1} = stress;
              if nargout > 1
                varargout{2} = displacement;
              end
            end
        end
        function plot_nodes_displaced(solver, displacements)
            %plot_nodes_displaced: Plots the deformed body.
            if nargin == 1
                displacements = solver.displacements;
            end
            nx = solver.nodes(:, 1) + displacements(:, 1);
            ny = solver.nodes(:, 2) + displacements(:, 2);
            scatter(nx, ny, ...
                10, 'filled');
            hold on;
            for e = solver.elements'
                c = convhull(nx(e), ny(e));
                plot(polyshape(nx(e(c)), ny(e(c)), ...
                    'Simplify', false));
            end
        end
        function plot_nodes(solver)
            %plot_nodes_displaced: Plots the undeformed body.
            solver.plot_nodes_displaced(0 * solver.nodes);
        end
    end
end
