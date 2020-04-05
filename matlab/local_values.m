function [K, F, M] = local_values(X, Y, u, l, f, local_nodes, int_nodes)
  %local_values Calculates the local stiffness, forces, and mass matrix and
  %for a given element and its type, using an appropriate integration method.
  if int_nodes == 4
    [K, F, M] = local_values_4(X, Y, u, l, f, local_nodes);
  else
    [K, F, M] = local_values_9(X, Y, u, l, f, local_nodes);
  end
end
