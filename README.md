# Finite Element Solver

This project is primarily run by the Solver class. A summary of its properties and attributes are as follows:

```matlab
Solver
    Properties
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

    Constant properties
        dims = cast(2, 'uint16');
        stress_components = cast(3, 'uint16');

    Constructor
        function solver = Solver(filename, overrides)
            %Solver: Class to solve a finite element problem, given an input
            % file.

    Helpers
        function crunch_local_data(solver, ID, M, nen)
            %crunch_local_data: Function to determine local stiffness, force,
            % and mass contributions and assemble them.
        function contextulize_d(solver, d, ID, M)
        function elements = fix_elements(nodes, old_elements, nen)
            %fix_elements: Ensure that convention for node numbering is
            % respected. I don't think this is violated in our examples, but
            % something Adyota said left me paranoid.
            elements = old_elements * 0;
        function M = build_map(dims, disp_BC)
            %build_map: Creates hashmap for boundary condition lookup.
        function ID = build_id(M, dims, nnd)
            %build_id: Creates the ID matrix for a given
            % system.

    Plotting functions
        function [varargout] = contour_stress(solver, resolution)
            %contour_stress: Creates a mesh grid of stress results plotted at
            % the respective displacments. Called with no return, this will
            % plot.
        function plot_nodes_displaced(solver, displacements)
            %plot_nodes_displaced: Plots the deformed body.
        function plot_nodes(solver)
            %plot_nodes_displaced: Plots the undeformed body.
        function plot_axis(solver)
            %plot_axis: Find and plot the y displacemnt at the center
            % line. (For comparison with Euler Bernoulli)
        function plot_centerline(solver)
            %plot_centerline: Find and plot the displacment at the center
            % line.
        function plot_midline_xx(solver)
            %plot_midline_xx: plots the xx stress along the center of the
            %beam vertically.
        function plot_midline_xy(solver)
            %plot_midline_xy: plots the xy stress along the center of the
            %beam vertically.
        function val = get_nearest_node_val(solver, branch, integration, A)
            %get_nearest_node_val: Branch determines which value to return.
            % 0 for displacement, 1 for stress. Integration dictates what
            % grid points to look at. A is the point to search for.

    Plotting Helpers
        function stress_fn = get_interpolate_stress_fn(solver, element)
            %get_interpolate_stress_fn: returns a function that can calculate
            % the stress at a given point for an element.
        function stress = interpolate_stress(solver, x, y)
            %interpolate_stress: Finds the element that contains x, y and
            % calculates the stresses at that point in the element.
        function displacement_fn = get_interpolate_displacement_fn(solver, element)
            %get_interpolate_displacement_fn: returns a function that can calculate
            % the displacement at a given point for an element.
        function displacement = interpolate_displacement(solver, x, y)
            %interpolate_displacement: Finds the element that contains x, y and
            % calculates the displacement at that point in the element.

    Other:
        function write_files(solver)
            %write_files Writes out important details to files
```

Usage would look something like:
```matlab
% Sample case
biaxial_Q4_2x2 = Solver("Biaxial_Q4_2x2.txt");
biaxial_Q4_2x2.plot_nodes_displaced();
biaxial_Q4_2x2.contour_stress();
biaxial_Q4_2x2.write_files();
```

where `"Biaxial_Q4_2x2.txt"` is the input file.

NOTE: Output files may be produced with the `write_files` function.

*-node files are csv files ordered by node number where each line is
`[posx posy] [dispx dispy] [sigma_xx sigma_yy sigma_xy]`

*-points files are a list of integration points in order of each element.
Each line contains
`[posx posy] [dispx dispy] [sigma_xx sigma_yy sigma_xy]`