function [P] = local_projection(X, Y, node, u, l, local_nodes, int_nodes)
  %local_stiffness Calculates the local projection vector for a given element
  % and its type, using an appropriate integration method.
  % NOTE: THIS IS GENERATED CODE. REFER TO local_projection.m.tmpl and generate.py
  if int_nodes == 4
    [P] = local_projection_4(X, Y, node, u, l, local_nodes);
  else
    [P] = local_projection_9(X, Y, node, u, l, local_nodes);
  end
end
