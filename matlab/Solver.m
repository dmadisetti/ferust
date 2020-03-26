classdef Solver < handle
    properties (SetAccess = private)
        displacements
        forces
    end
    methods
        function solver = Solver(library_path, header_path)
            persistent loaded library_cached header_cached;
            if nargin < 2
                header_path = 'src\ferust.h';
                if nargin == 0
                    library_path = 'target\debug\ferust.dll';
                end
            end
            if ~isempty(loaded) && ...
                    ~(strcmp(library_cached, library_path) && ...
                      strcmp(header_cached, header_path))
                    unloadlibrary libferust
            end
            loadlibrary(library_path, header_path, 'alias', 'libferust');
            loaded = 1;
            library_cached = library_path;
            header_cached = header_path;
        end
        function solve(solver, nodes, elements, ~, element_count, ...
                element_nodes, integration_points, node_count, ...
                plane_stress, nu, E, reactions, body_force, ...
                displacements)
            solver.displacements = zeros(size(nodes));
            solver.forces = zeros(size(nodes));
            displacement_nodes = cast(displacements(:, 1:2) - 1, 'uint16');
            displacements = displacements(:, 3);

            d = libpointer('doublePtr', solver.displacements);
            f = libpointer('doublePtr', solver.forces);
            % Note that we scramble the arguments into something a bit more
            % sane to our rust function.
            calllib('libferust', 'solve', element_count, element_nodes, ...
                integration_points, node_count, plane_stress, ...
                length(displacements), length(reactions), nu, E, ...
                body_force, nodes, elements - 1, reactions - 1, ...
                displacement_nodes, displacements, d, f);
        end
        function delete(solver)
            unloadlibrary libferust
        end
    end
end