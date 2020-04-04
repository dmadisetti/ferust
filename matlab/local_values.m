function [K, F, M] = local_values(X, Y, u, l, f, local_nodes, int_nodes)
  %local_stiffness Calculates the local stiffness, forces, mass matrix and
  % projection vector for a given element and its type, using an appropriate
  % integration method.
  % NOTE: THIS IS GENERATED CODE. REFER TO local_values.m.tmpl and generate.py
  if int_nodes == 4
    [K, F, M] = local_values_4(X, Y, u, l, f, local_nodes);
  else
    [K, F, M] = local_values_9(X, Y, u, l, f, local_nodes);
  end
end
