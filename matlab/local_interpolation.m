function [fn] = local_interpolation(X, Y, interp, local_nodes)
  %local_interpolation Returns a function that calculates the interpolation
  % value in a given element provided solved values at the nodes using the shape
  % functions.
  % NOTE: THIS IS GENERATED CODE. REFER TO local_stress.m.tmpl and generate.py
  if local_nodes == 4
    [fn] = crunch_quad(X(1), X(2), X(3), X(4), ...
                      Y(1), Y(2), Y(3), Y(4), ...
                      interp(1), interp(2), interp(3), interp(4));
  elseif local_nodes == 8
    [fn] = crunch_serindipity(X(1), X(2), X(3), X(4), ...
                           X(5), X(6), X(7), X(8), ...
                           Y(1), Y(2), Y(3), Y(4), ...
                           Y(5), Y(6), Y(7), Y(8), ...
                           interp(1), interp(2), interp(3), interp(4), ...
                           interp(5), interp(6), interp(7), interp(8));
  else
    [fn] = crunch_bubble(X(1), X(2), X(3), X(4), ...
                           X(5), X(6), X(7), X(8), X(9), ...
                           Y(1), Y(2), Y(3), Y(4), Y(5), ...
                           Y(6), Y(7), Y(8), Y(9), ...
                           interp(1), interp(2), interp(3), interp(4), ...
                           interp(5), interp(6), interp(7), interp(8), ...
                           interp(9));
  end
end

function [fn] = crunch_quad(A, B, C, D, ...
                           R, S, T, U, ...
                           d0, d1, d2, d3)
  
    fn = @(n, e) (d0*(1 - e)*(1 - n)/4 + d1*(1 - n)*(e + 1)/4 + d2*(e + 1)*(n + 1)/4 + d3*(1 - ...
     e)*(n + 1)/4);
end

function [fn] = crunch_serindipity(A, B, C, D, E, F, G, H, ...
                                  R, S, T, U, V, W, X, Y, ...
                                  d0, d1, d2, d3, d4, d5, d6, d7)
  
    fn = @(n, e) (d0*((1 - e)*(1 - n)/4 - (1 - e)*(1 - n^2)/4 - (1 - e^2)*(1 - n)/4) + d1*(- ...
    (1 - e^2)*(1 - n)/4 + (1 - n)*(e + 1)/4 - (1 - n^2)*(e + 1)/4) + d2*(-(1 - e^2 ...
    )*(n + 1)/4 - (1 - n^2)*(e + 1)/4 + (e + 1)*(n + 1)/4) + d3*(-(1 - e)*(1 - n^2 ...
    )/4 + (1 - e)*(n + 1)/4 - (1 - e^2)*(n + 1)/4) + d4*(1 - e^2)*(1 - n)/2 + d5*(1 ...
     - n^2)*(e + 1)/2 + d6*(1 - e^2)*(n + 1)/2 + d7*(1 - e)*(1 - n^2)/2)*(-3*A*S*e^2 ...
    *n^2/8 + 5*A*S*e^2*n/8 - A*S*e^2/4 - A*S*n^3/4 + A*S*n^2/4 + A*T*e^3/4 + 3*A* ...
    T*e^2*n/8 - 3*A*T*e*n^2/8 - A*T*n^3/4 + A*U*e^3/4 + 3*A*U*e^2*n^2/8 - A*U*e^2 ...
    /4 - 5*A*U*e*n^2/8 + A*U*n^2/4 + 3*A*V*e^2*n^2/8 - 5*A*V*e^2*n/8 + A*V*e^2/4  ...
    - A*V*e*n^2/2 + 3*A*V*e*n/4 - A*V*e/4 + A*V*n^2/8 - A*V*n/8 + 3*A*W*e^2*n^2/8 ...
     - A*W*e^2*n/2 + A*W*e^2/8 + 3*A*W*e*n^2/8 - A*W*e*n/4 - A*W*e/8 + A*W*n^3/2  ...
    - A*W*n^2/4 - A*W*n/4 - A*X*e^3/2 - 3*A*X*e^2*n^2/8 - 3*A*X*e^2*n/8 + A*X*e^2 ...
    /4 + A*X*e*n^2/2 + A*X*e*n/4 + A*X*e/4 - A*X*n^2/8 + A*X*n/8 - 3*A*Y*e^2*n^2/ ...
    8 + A*Y*e^2*n/2 - A*Y*e^2/8 + 5*A*Y*e*n^2/8 - 3*A*Y*e*n/4 + A*Y*e/8 - A*Y*n^2 ...
    /4 + A*Y*n/4 + 3*B*R*e^2*n^2/8 - 5*B*R*e^2*n/8 + B*R*e^2/4 + B*R*n^3/4 - B*R* ...
    n^2/4 + B*T*e^3/4 - 3*B*T*e^2*n^2/8 + B*T*e^2/4 - 5*B*T*e*n^2/8 - B*T*n^2/4 + ...
     B*U*e^3/4 - 3*B*U*e^2*n/8 - 3*B*U*e*n^2/8 + B*U*n^3/4 - 3*B*V*e^2*n^2/8 + 5* ...
    B*V*e^2*n/8 - B*V*e^2/4 - B*V*e*n^2/2 + 3*B*V*e*n/4 - B*V*e/4 - B*V*n^2/8 + B ...
    *V*n/8 + 3*B*W*e^2*n^2/8 - B*W*e^2*n/2 + B*W*e^2/8 + 5*B*W*e*n^2/8 - 3*B*W*e* ...
    n/4 + B*W*e/8 + B*W*n^2/4 - B*W*n/4 - B*X*e^3/2 + 3*B*X*e^2*n^2/8 + 3*B*X*e^2 ...
    *n/8 - B*X*e^2/4 + B*X*e*n^2/2 + B*X*e*n/4 + B*X*e/4 + B*X*n^2/8 - B*X*n/8 -  ...
    3*B*Y*e^2*n^2/8 + B*Y*e^2*n/2 - B*Y*e^2/8 + 3*B*Y*e*n^2/8 - B*Y*e*n/4 - B*Y*e ...
    /8 - B*Y*n^3/2 + B*Y*n^2/4 + B*Y*n/4 - C*R*e^3/4 - 3*C*R*e^2*n/8 + 3*C*R*e*n^ ...
    2/8 + C*R*n^3/4 - C*S*e^3/4 + 3*C*S*e^2*n^2/8 - C*S*e^2/4 + 5*C*S*e*n^2/8 + C ...
    *S*n^2/4 - 3*C*U*e^2*n^2/8 - 5*C*U*e^2*n/8 - C*U*e^2/4 + C*U*n^3/4 + C*U*n^2/ ...
    4 + C*V*e^3/2 - 3*C*V*e^2*n^2/8 + 3*C*V*e^2*n/8 + C*V*e^2/4 - C*V*e*n^2/2 + C ...
    *V*e*n/4 - C*V*e/4 - C*V*n^2/8 - C*V*n/8 - 3*C*W*e^2*n^2/8 - C*W*e^2*n/2 - C* ...
    W*e^2/8 - 5*C*W*e*n^2/8 - 3*C*W*e*n/4 - C*W*e/8 - C*W*n^2/4 - C*W*n/4 + 3*C*X ...
    *e^2*n^2/8 + 5*C*X*e^2*n/8 + C*X*e^2/4 + C*X*e*n^2/2 + 3*C*X*e*n/4 + C*X*e/4  ...
    + C*X*n^2/8 + C*X*n/8 + 3*C*Y*e^2*n^2/8 + C*Y*e^2*n/2 + C*Y*e^2/8 - 3*C*Y*e*n ...
    ^2/8 - C*Y*e*n/4 + C*Y*e/8 - C*Y*n^3/2 - C*Y*n^2/4 + C*Y*n/4 - D*R*e^3/4 - 3* ...
    D*R*e^2*n^2/8 + D*R*e^2/4 + 5*D*R*e*n^2/8 - D*R*n^2/4 - D*S*e^3/4 + 3*D*S*e^2 ...
    *n/8 + 3*D*S*e*n^2/8 - D*S*n^3/4 + 3*D*T*e^2*n^2/8 + 5*D*T*e^2*n/8 + D*T*e^2/ ...
    4 - D*T*n^3/4 - D*T*n^2/4 + D*V*e^3/2 + 3*D*V*e^2*n^2/8 - 3*D*V*e^2*n/8 - D*V ...
    *e^2/4 - D*V*e*n^2/2 + D*V*e*n/4 - D*V*e/4 + D*V*n^2/8 + D*V*n/8 - 3*D*W*e^2* ...
    n^2/8 - D*W*e^2*n/2 - D*W*e^2/8 - 3*D*W*e*n^2/8 - D*W*e*n/4 + D*W*e/8 + D*W*n ...
    ^3/2 + D*W*n^2/4 - D*W*n/4 - 3*D*X*e^2*n^2/8 - 5*D*X*e^2*n/8 - D*X*e^2/4 + D* ...
    X*e*n^2/2 + 3*D*X*e*n/4 + D*X*e/4 - D*X*n^2/8 - D*X*n/8 + 3*D*Y*e^2*n^2/8 + D ...
    *Y*e^2*n/2 + D*Y*e^2/8 - 5*D*Y*e*n^2/8 - 3*D*Y*e*n/4 - D*Y*e/8 + D*Y*n^2/4 +  ...
    D*Y*n/4 - 3*E*R*e^2*n^2/8 + 5*E*R*e^2*n/8 - E*R*e^2/4 + E*R*e*n^2/2 - 3*E*R*e ...
    *n/4 + E*R*e/4 - E*R*n^2/8 + E*R*n/8 + 3*E*S*e^2*n^2/8 - 5*E*S*e^2*n/8 + E*S* ...
    e^2/4 + E*S*e*n^2/2 - 3*E*S*e*n/4 + E*S*e/4 + E*S*n^2/8 - E*S*n/8 - E*T*e^3/2 ...
     + 3*E*T*e^2*n^2/8 - 3*E*T*e^2*n/8 - E*T*e^2/4 + E*T*e*n^2/2 - E*T*e*n/4 + E* ...
    T*e/4 + E*T*n^2/8 + E*T*n/8 - E*U*e^3/2 - 3*E*U*e^2*n^2/8 + 3*E*U*e^2*n/8 + E ...
    *U*e^2/4 + E*U*e*n^2/2 - E*U*e*n/4 + E*U*e/4 - E*U*n^2/8 - E*U*n/8 - 3*E*W*e^ ...
    2*n^2/4 + E*W*e^2*n - E*W*e^2/4 - E*W*e*n^2 + E*W*e*n - E*W*n^2/4 + E*W/4 + E ...
    *X*e^3 - E*X*e + 3*E*Y*e^2*n^2/4 - E*Y*e^2*n + E*Y*e^2/4 - E*Y*e*n^2 + E*Y*e* ...
    n + E*Y*n^2/4 - E*Y/4 - 3*F*R*e^2*n^2/8 + F*R*e^2*n/2 - F*R*e^2/8 - 3*F*R*e*n ...
    ^2/8 + F*R*e*n/4 + F*R*e/8 - F*R*n^3/2 + F*R*n^2/4 + F*R*n/4 - 3*F*S*e^2*n^2/ ...
    8 + F*S*e^2*n/2 - F*S*e^2/8 - 5*F*S*e*n^2/8 + 3*F*S*e*n/4 - F*S*e/8 - F*S*n^2 ...
    /4 + F*S*n/4 + 3*F*T*e^2*n^2/8 + F*T*e^2*n/2 + F*T*e^2/8 + 5*F*T*e*n^2/8 + 3* ...
    F*T*e*n/4 + F*T*e/8 + F*T*n^2/4 + F*T*n/4 + 3*F*U*e^2*n^2/8 + F*U*e^2*n/2 + F ...
    *U*e^2/8 + 3*F*U*e*n^2/8 + F*U*e*n/4 - F*U*e/8 - F*U*n^3/2 - F*U*n^2/4 + F*U* ...
    n/4 + 3*F*V*e^2*n^2/4 - F*V*e^2*n + F*V*e^2/4 + F*V*e*n^2 - F*V*e*n + F*V*n^2 ...
    /4 - F*V/4 - 3*F*X*e^2*n^2/4 - F*X*e^2*n - F*X*e^2/4 - F*X*e*n^2 - F*X*e*n -  ...
    F*X*n^2/4 + F*X/4 + F*Y*n^3 - F*Y*n + G*R*e^3/2 + 3*G*R*e^2*n^2/8 + 3*G*R*e^2 ...
    *n/8 - G*R*e^2/4 - G*R*e*n^2/2 - G*R*e*n/4 - G*R*e/4 + G*R*n^2/8 - G*R*n/8 +  ...
    G*S*e^3/2 - 3*G*S*e^2*n^2/8 - 3*G*S*e^2*n/8 + G*S*e^2/4 - G*S*e*n^2/2 - G*S*e ...
    *n/4 - G*S*e/4 - G*S*n^2/8 + G*S*n/8 - 3*G*T*e^2*n^2/8 - 5*G*T*e^2*n/8 - G*T* ...
    e^2/4 - G*T*e*n^2/2 - 3*G*T*e*n/4 - G*T*e/4 - G*T*n^2/8 - G*T*n/8 + 3*G*U*e^2 ...
    *n^2/8 + 5*G*U*e^2*n/8 + G*U*e^2/4 - G*U*e*n^2/2 - 3*G*U*e*n/4 - G*U*e/4 + G* ...
    U*n^2/8 + G*U*n/8 - G*V*e^3 + G*V*e + 3*G*W*e^2*n^2/4 + G*W*e^2*n + G*W*e^2/4 ...
     + G*W*e*n^2 + G*W*e*n + G*W*n^2/4 - G*W/4 - 3*G*Y*e^2*n^2/4 - G*Y*e^2*n - G* ...
    Y*e^2/4 + G*Y*e*n^2 + G*Y*e*n - G*Y*n^2/4 + G*Y/4 + 3*H*R*e^2*n^2/8 - H*R*e^2 ...
    *n/2 + H*R*e^2/8 - 5*H*R*e*n^2/8 + 3*H*R*e*n/4 - H*R*e/8 + H*R*n^2/4 - H*R*n/ ...
    4 + 3*H*S*e^2*n^2/8 - H*S*e^2*n/2 + H*S*e^2/8 - 3*H*S*e*n^2/8 + H*S*e*n/4 + H ...
    *S*e/8 + H*S*n^3/2 - H*S*n^2/4 - H*S*n/4 - 3*H*T*e^2*n^2/8 - H*T*e^2*n/2 - H* ...
    T*e^2/8 + 3*H*T*e*n^2/8 + H*T*e*n/4 - H*T*e/8 + H*T*n^3/2 + H*T*n^2/4 - H*T*n ...
    /4 - 3*H*U*e^2*n^2/8 - H*U*e^2*n/2 - H*U*e^2/8 + 5*H*U*e*n^2/8 + 3*H*U*e*n/4  ...
    + H*U*e/8 - H*U*n^2/4 - H*U*n/4 - 3*H*V*e^2*n^2/4 + H*V*e^2*n - H*V*e^2/4 + H ...
    *V*e*n^2 - H*V*e*n - H*V*n^2/4 + H*V/4 - H*W*n^3 + H*W*n + 3*H*X*e^2*n^2/4 +  ...
    H*X*e^2*n + H*X*e^2/4 - H*X*e*n^2 - H*X*e*n + H*X*n^2/4 - H*X/4);;
end

function [fn] = crunch_bubble(A, B, C, D, E, F, G, H, I, ...
                             R, S, T, U, V, W, X, Y, Z, ...
                             d0, d1, d2, d3, d4, d5, d6, d7, d8)
  
    fn = @(n, e) (d0*((1 - e)*(1 - n)/4 - (1 - e)*(1 - n^2)/4 - (1 - e^2)*(1 - n)/4 + (1 -  ...
    e^2)*(1 - n^2)/4) + d1*(-(1 - e^2)*(1 - n)/4 + (1 - e^2)*(1 - n^2)/4 + (1 - n) ...
    *(e + 1)/4 - (1 - n^2)*(e + 1)/4) + d2*((1 - e^2)*(1 - n^2)/4 - (1 - e^2)*(n + ...
     1)/4 - (1 - n^2)*(e + 1)/4 + (e + 1)*(n + 1)/4) + d3*(-(1 - e)*(1 - n^2)/4 +  ...
    (1 - e)*(n + 1)/4 + (1 - e^2)*(1 - n^2)/4 - (1 - e^2)*(n + 1)/4) + d4*((1 - e^ ...
    2)*(1 - n)/2 - (1 - e^2)*(1 - n^2)/2) + d5*(-(1 - e^2)*(1 - n^2)/2 + (1 - n^2) ...
    *(e + 1)/2) + d6*(-(1 - e^2)*(1 - n^2)/2 + (1 - e^2)*(n + 1)/2) + d7*((1 - e)*( ...
    1 - n^2)/2 - (1 - e^2)*(1 - n^2)/2) + d8*(1 - e^2)*(1 - n^2))*(A*S*e^2*n^3/4 - ...
     3*A*S*e^2*n^2/8 + A*S*e^2*n/8 - A*T*e^3*n^2/4 + A*T*e^2*n^3/4 - A*T*e^2*n/8  ...
    + A*T*e*n^2/8 - A*U*e^3*n^2/4 + 3*A*U*e^2*n^2/8 - A*U*e*n^2/8 - A*V*e^2*n^3/4 ...
     + 3*A*V*e^2*n^2/8 - A*V*e^2*n/8 + A*V*e*n^3/2 - 3*A*V*e*n^2/4 + A*V*e*n/4 -  ...
    A*V*n^3/4 + 3*A*V*n^2/8 - A*V*n/8 + A*W*e^3*n^2/4 - A*W*e^3*n/2 + A*W*e^3/4 - ...
     A*W*e^2*n^3/2 + 3*A*W*e^2*n^2/8 + A*W*e^2*n/4 - A*W*e^2/8 - A*W*e*n^2/8 + A* ...
    W*e*n/4 - A*W*e/8 + A*X*e^3*n^2/2 - A*X*e^2*n^3/4 - 3*A*X*e^2*n^2/8 + A*X*e^2 ...
    *n/8 + A*X*e*n^3/2 - A*X*e*n^2/4 - A*X*e*n/4 - A*X*n^3/4 + A*X*n^2/8 + A*X*n/ ...
    8 + A*Y*e^3*n^2/4 - A*Y*e^3*n/2 + A*Y*e^3/4 - 3*A*Y*e^2*n^2/8 + 3*A*Y*e^2*n/4 ...
     - 3*A*Y*e^2/8 + A*Y*e*n^2/8 - A*Y*e*n/4 + A*Y*e/8 - A*Z*e^3*n^2/2 + A*Z*e^3* ...
    n - A*Z*e^3/2 + A*Z*e^2*n^3/2 - A*Z*e^2*n + A*Z*e^2/2 - A*Z*e*n^3 + A*Z*e*n^2 ...
     + A*Z*n^3/2 - A*Z*n^2/2 - B*R*e^2*n^3/4 + 3*B*R*e^2*n^2/8 - B*R*e^2*n/8 - B* ...
    T*e^3*n^2/4 - 3*B*T*e^2*n^2/8 - B*T*e*n^2/8 - B*U*e^3*n^2/4 - B*U*e^2*n^3/4 + ...
     B*U*e^2*n/8 + B*U*e*n^2/8 + B*V*e^2*n^3/4 - 3*B*V*e^2*n^2/8 + B*V*e^2*n/8 +  ...
    B*V*e*n^3/2 - 3*B*V*e*n^2/4 + B*V*e*n/4 + B*V*n^3/4 - 3*B*V*n^2/8 + B*V*n/8 + ...
     B*W*e^3*n^2/4 - B*W*e^3*n/2 + B*W*e^3/4 + 3*B*W*e^2*n^2/8 - 3*B*W*e^2*n/4 +  ...
    3*B*W*e^2/8 + B*W*e*n^2/8 - B*W*e*n/4 + B*W*e/8 + B*X*e^3*n^2/2 + B*X*e^2*n^3 ...
    /4 + 3*B*X*e^2*n^2/8 - B*X*e^2*n/8 + B*X*e*n^3/2 - B*X*e*n^2/4 - B*X*e*n/4 +  ...
    B*X*n^3/4 - B*X*n^2/8 - B*X*n/8 + B*Y*e^3*n^2/4 - B*Y*e^3*n/2 + B*Y*e^3/4 + B ...
    *Y*e^2*n^3/2 - 3*B*Y*e^2*n^2/8 - B*Y*e^2*n/4 + B*Y*e^2/8 - B*Y*e*n^2/8 + B*Y* ...
    e*n/4 - B*Y*e/8 - B*Z*e^3*n^2/2 + B*Z*e^3*n - B*Z*e^3/2 - B*Z*e^2*n^3/2 + B*Z ...
    *e^2*n - B*Z*e^2/2 - B*Z*e*n^3 + B*Z*e*n^2 - B*Z*n^3/2 + B*Z*n^2/2 + C*R*e^3* ...
    n^2/4 - C*R*e^2*n^3/4 + C*R*e^2*n/8 - C*R*e*n^2/8 + C*S*e^3*n^2/4 + 3*C*S*e^2 ...
    *n^2/8 + C*S*e*n^2/8 - C*U*e^2*n^3/4 - 3*C*U*e^2*n^2/8 - C*U*e^2*n/8 - C*V*e^ ...
    3*n^2/2 + C*V*e^2*n^3/4 - 3*C*V*e^2*n^2/8 - C*V*e^2*n/8 + C*V*e*n^3/2 + C*V*e ...
    *n^2/4 - C*V*e*n/4 + C*V*n^3/4 + C*V*n^2/8 - C*V*n/8 - C*W*e^3*n^2/4 - C*W*e^ ...
    3*n/2 - C*W*e^3/4 - 3*C*W*e^2*n^2/8 - 3*C*W*e^2*n/4 - 3*C*W*e^2/8 - C*W*e*n^2 ...
    /8 - C*W*e*n/4 - C*W*e/8 + C*X*e^2*n^3/4 + 3*C*X*e^2*n^2/8 + C*X*e^2*n/8 + C* ...
    X*e*n^3/2 + 3*C*X*e*n^2/4 + C*X*e*n/4 + C*X*n^3/4 + 3*C*X*n^2/8 + C*X*n/8 - C ...
    *Y*e^3*n^2/4 - C*Y*e^3*n/2 - C*Y*e^3/4 + C*Y*e^2*n^3/2 + 3*C*Y*e^2*n^2/8 - C* ...
    Y*e^2*n/4 - C*Y*e^2/8 + C*Y*e*n^2/8 + C*Y*e*n/4 + C*Y*e/8 + C*Z*e^3*n^2/2 + C ...
    *Z*e^3*n + C*Z*e^3/2 - C*Z*e^2*n^3/2 + C*Z*e^2*n + C*Z*e^2/2 - C*Z*e*n^3 - C* ...
    Z*e*n^2 - C*Z*n^3/2 - C*Z*n^2/2 + D*R*e^3*n^2/4 - 3*D*R*e^2*n^2/8 + D*R*e*n^2 ...
    /8 + D*S*e^3*n^2/4 + D*S*e^2*n^3/4 - D*S*e^2*n/8 - D*S*e*n^2/8 + D*T*e^2*n^3/ ...
    4 + 3*D*T*e^2*n^2/8 + D*T*e^2*n/8 - D*V*e^3*n^2/2 - D*V*e^2*n^3/4 + 3*D*V*e^2 ...
    *n^2/8 + D*V*e^2*n/8 + D*V*e*n^3/2 + D*V*e*n^2/4 - D*V*e*n/4 - D*V*n^3/4 - D* ...
    V*n^2/8 + D*V*n/8 - D*W*e^3*n^2/4 - D*W*e^3*n/2 - D*W*e^3/4 - D*W*e^2*n^3/2 - ...
     3*D*W*e^2*n^2/8 + D*W*e^2*n/4 + D*W*e^2/8 + D*W*e*n^2/8 + D*W*e*n/4 + D*W*e/ ...
    8 - D*X*e^2*n^3/4 - 3*D*X*e^2*n^2/8 - D*X*e^2*n/8 + D*X*e*n^3/2 + 3*D*X*e*n^2 ...
    /4 + D*X*e*n/4 - D*X*n^3/4 - 3*D*X*n^2/8 - D*X*n/8 - D*Y*e^3*n^2/4 - D*Y*e^3* ...
    n/2 - D*Y*e^3/4 + 3*D*Y*e^2*n^2/8 + 3*D*Y*e^2*n/4 + 3*D*Y*e^2/8 - D*Y*e*n^2/8 ...
     - D*Y*e*n/4 - D*Y*e/8 + D*Z*e^3*n^2/2 + D*Z*e^3*n + D*Z*e^3/2 + D*Z*e^2*n^3/ ...
    2 - D*Z*e^2*n - D*Z*e^2/2 - D*Z*e*n^3 - D*Z*e*n^2 + D*Z*n^3/2 + D*Z*n^2/2 + E ...
    *R*e^2*n^3/4 - 3*E*R*e^2*n^2/8 + E*R*e^2*n/8 - E*R*e*n^3/2 + 3*E*R*e*n^2/4 -  ...
    E*R*e*n/4 + E*R*n^3/4 - 3*E*R*n^2/8 + E*R*n/8 - E*S*e^2*n^3/4 + 3*E*S*e^2*n^2 ...
    /8 - E*S*e^2*n/8 - E*S*e*n^3/2 + 3*E*S*e*n^2/4 - E*S*e*n/4 - E*S*n^3/4 + 3*E* ...
    S*n^2/8 - E*S*n/8 + E*T*e^3*n^2/2 - E*T*e^2*n^3/4 + 3*E*T*e^2*n^2/8 + E*T*e^2 ...
    *n/8 - E*T*e*n^3/2 - E*T*e*n^2/4 + E*T*e*n/4 - E*T*n^3/4 - E*T*n^2/8 + E*T*n/ ...
    8 + E*U*e^3*n^2/2 + E*U*e^2*n^3/4 - 3*E*U*e^2*n^2/8 - E*U*e^2*n/8 - E*U*e*n^3 ...
    /2 - E*U*e*n^2/4 + E*U*e*n/4 + E*U*n^3/4 + E*U*n^2/8 - E*U*n/8 - E*W*e^3*n^2/ ...
    2 + E*W*e^3*n - E*W*e^3/2 + E*W*e^2*n^3/2 - 3*E*W*e^2*n^2/4 + E*W*e^2*n/2 - E ...
    *W*e^2/4 + E*W*e*n^3 - E*W*e*n^2/2 - E*W*e*n + E*W*e/2 + E*W*n^3/2 - E*W*n^2/ ...
    4 - E*W*n/2 + E*W/4 - E*X*e^3*n^2 + E*X*e*n^2 - E*Y*e^3*n^2/2 + E*Y*e^3*n - E ...
    *Y*e^3/2 - E*Y*e^2*n^3/2 + 3*E*Y*e^2*n^2/4 - E*Y*e^2*n/2 + E*Y*e^2/4 + E*Y*e* ...
    n^3 - E*Y*e*n^2/2 - E*Y*e*n + E*Y*e/2 - E*Y*n^3/2 + E*Y*n^2/4 + E*Y*n/2 - E*Y ...
    /4 + E*Z*e^3*n^2 - 2*E*Z*e^3*n + E*Z*e^3 - E*Z*e*n^2 + 2*E*Z*e*n - E*Z*e - F* ...
    R*e^3*n^2/4 + F*R*e^3*n/2 - F*R*e^3/4 + F*R*e^2*n^3/2 - 3*F*R*e^2*n^2/8 - F*R ...
    *e^2*n/4 + F*R*e^2/8 + F*R*e*n^2/8 - F*R*e*n/4 + F*R*e/8 - F*S*e^3*n^2/4 + F* ...
    S*e^3*n/2 - F*S*e^3/4 - 3*F*S*e^2*n^2/8 + 3*F*S*e^2*n/4 - 3*F*S*e^2/8 - F*S*e ...
    *n^2/8 + F*S*e*n/4 - F*S*e/8 + F*T*e^3*n^2/4 + F*T*e^3*n/2 + F*T*e^3/4 + 3*F* ...
    T*e^2*n^2/8 + 3*F*T*e^2*n/4 + 3*F*T*e^2/8 + F*T*e*n^2/8 + F*T*e*n/4 + F*T*e/8 ...
     + F*U*e^3*n^2/4 + F*U*e^3*n/2 + F*U*e^3/4 + F*U*e^2*n^3/2 + 3*F*U*e^2*n^2/8  ...
    - F*U*e^2*n/4 - F*U*e^2/8 - F*U*e*n^2/8 - F*U*e*n/4 - F*U*e/8 + F*V*e^3*n^2/2 ...
     - F*V*e^3*n + F*V*e^3/2 - F*V*e^2*n^3/2 + 3*F*V*e^2*n^2/4 - F*V*e^2*n/2 + F* ...
    V*e^2/4 - F*V*e*n^3 + F*V*e*n^2/2 + F*V*e*n - F*V*e/2 - F*V*n^3/2 + F*V*n^2/4 ...
     + F*V*n/2 - F*V/4 - F*X*e^3*n^2/2 - F*X*e^3*n - F*X*e^3/2 - F*X*e^2*n^3/2 -  ...
    3*F*X*e^2*n^2/4 - F*X*e^2*n/2 - F*X*e^2/4 - F*X*e*n^3 - F*X*e*n^2/2 + F*X*e*n ...
     + F*X*e/2 - F*X*n^3/2 - F*X*n^2/4 + F*X*n/2 + F*X/4 - F*Y*e^2*n^3 + F*Y*e^2* ...
    n + F*Z*e^2*n^3 - F*Z*e^2*n + 2*F*Z*e*n^3 - 2*F*Z*e*n + F*Z*n^3 - F*Z*n - G*R ...
    *e^3*n^2/2 + G*R*e^2*n^3/4 + 3*G*R*e^2*n^2/8 - G*R*e^2*n/8 - G*R*e*n^3/2 + G* ...
    R*e*n^2/4 + G*R*e*n/4 + G*R*n^3/4 - G*R*n^2/8 - G*R*n/8 - G*S*e^3*n^2/2 - G*S ...
    *e^2*n^3/4 - 3*G*S*e^2*n^2/8 + G*S*e^2*n/8 - G*S*e*n^3/2 + G*S*e*n^2/4 + G*S* ...
    e*n/4 - G*S*n^3/4 + G*S*n^2/8 + G*S*n/8 - G*T*e^2*n^3/4 - 3*G*T*e^2*n^2/8 - G ...
    *T*e^2*n/8 - G*T*e*n^3/2 - 3*G*T*e*n^2/4 - G*T*e*n/4 - G*T*n^3/4 - 3*G*T*n^2/ ...
    8 - G*T*n/8 + G*U*e^2*n^3/4 + 3*G*U*e^2*n^2/8 + G*U*e^2*n/8 - G*U*e*n^3/2 - 3 ...
    *G*U*e*n^2/4 - G*U*e*n/4 + G*U*n^3/4 + 3*G*U*n^2/8 + G*U*n/8 + G*V*e^3*n^2 -  ...
    G*V*e*n^2 + G*W*e^3*n^2/2 + G*W*e^3*n + G*W*e^3/2 + G*W*e^2*n^3/2 + 3*G*W*e^2 ...
    *n^2/4 + G*W*e^2*n/2 + G*W*e^2/4 + G*W*e*n^3 + G*W*e*n^2/2 - G*W*e*n - G*W*e/ ...
    2 + G*W*n^3/2 + G*W*n^2/4 - G*W*n/2 - G*W/4 + G*Y*e^3*n^2/2 + G*Y*e^3*n + G*Y ...
    *e^3/2 - G*Y*e^2*n^3/2 - 3*G*Y*e^2*n^2/4 - G*Y*e^2*n/2 - G*Y*e^2/4 + G*Y*e*n^ ...
    3 + G*Y*e*n^2/2 - G*Y*e*n - G*Y*e/2 - G*Y*n^3/2 - G*Y*n^2/4 + G*Y*n/2 + G*Y/4 ...
     - G*Z*e^3*n^2 - 2*G*Z*e^3*n - G*Z*e^3 + G*Z*e*n^2 + 2*G*Z*e*n + G*Z*e - H*R* ...
    e^3*n^2/4 + H*R*e^3*n/2 - H*R*e^3/4 + 3*H*R*e^2*n^2/8 - 3*H*R*e^2*n/4 + 3*H*R ...
    *e^2/8 - H*R*e*n^2/8 + H*R*e*n/4 - H*R*e/8 - H*S*e^3*n^2/4 + H*S*e^3*n/2 - H* ...
    S*e^3/4 - H*S*e^2*n^3/2 + 3*H*S*e^2*n^2/8 + H*S*e^2*n/4 - H*S*e^2/8 + H*S*e*n ...
    ^2/8 - H*S*e*n/4 + H*S*e/8 + H*T*e^3*n^2/4 + H*T*e^3*n/2 + H*T*e^3/4 - H*T*e^ ...
    2*n^3/2 - 3*H*T*e^2*n^2/8 + H*T*e^2*n/4 + H*T*e^2/8 - H*T*e*n^2/8 - H*T*e*n/4 ...
     - H*T*e/8 + H*U*e^3*n^2/4 + H*U*e^3*n/2 + H*U*e^3/4 - 3*H*U*e^2*n^2/8 - 3*H* ...
    U*e^2*n/4 - 3*H*U*e^2/8 + H*U*e*n^2/8 + H*U*e*n/4 + H*U*e/8 + H*V*e^3*n^2/2 - ...
     H*V*e^3*n + H*V*e^3/2 + H*V*e^2*n^3/2 - 3*H*V*e^2*n^2/4 + H*V*e^2*n/2 - H*V* ...
    e^2/4 - H*V*e*n^3 + H*V*e*n^2/2 + H*V*e*n - H*V*e/2 + H*V*n^3/2 - H*V*n^2/4 - ...
     H*V*n/2 + H*V/4 + H*W*e^2*n^3 - H*W*e^2*n - H*X*e^3*n^2/2 - H*X*e^3*n - H*X* ...
    e^3/2 + H*X*e^2*n^3/2 + 3*H*X*e^2*n^2/4 + H*X*e^2*n/2 + H*X*e^2/4 - H*X*e*n^3 ...
     - H*X*e*n^2/2 + H*X*e*n + H*X*e/2 + H*X*n^3/2 + H*X*n^2/4 - H*X*n/2 - H*X/4  ...
    - H*Z*e^2*n^3 + H*Z*e^2*n + 2*H*Z*e*n^3 - 2*H*Z*e*n - H*Z*n^3 + H*Z*n + I*R*e ...
    ^3*n^2/2 - I*R*e^3*n + I*R*e^3/2 - I*R*e^2*n^3/2 + I*R*e^2*n - I*R*e^2/2 + I* ...
    R*e*n^3 - I*R*e*n^2 - I*R*n^3/2 + I*R*n^2/2 + I*S*e^3*n^2/2 - I*S*e^3*n + I*S ...
    *e^3/2 + I*S*e^2*n^3/2 - I*S*e^2*n + I*S*e^2/2 + I*S*e*n^3 - I*S*e*n^2 + I*S* ...
    n^3/2 - I*S*n^2/2 - I*T*e^3*n^2/2 - I*T*e^3*n - I*T*e^3/2 + I*T*e^2*n^3/2 - I ...
    *T*e^2*n - I*T*e^2/2 + I*T*e*n^3 + I*T*e*n^2 + I*T*n^3/2 + I*T*n^2/2 - I*U*e^ ...
    3*n^2/2 - I*U*e^3*n - I*U*e^3/2 - I*U*e^2*n^3/2 + I*U*e^2*n + I*U*e^2/2 + I*U ...
    *e*n^3 + I*U*e*n^2 - I*U*n^3/2 - I*U*n^2/2 - I*V*e^3*n^2 + 2*I*V*e^3*n - I*V* ...
    e^3 + I*V*e*n^2 - 2*I*V*e*n + I*V*e - I*W*e^2*n^3 + I*W*e^2*n - 2*I*W*e*n^3 + ...
     2*I*W*e*n - I*W*n^3 + I*W*n + I*X*e^3*n^2 + 2*I*X*e^3*n + I*X*e^3 - I*X*e*n^ ...
    2 - 2*I*X*e*n - I*X*e + I*Y*e^2*n^3 - I*Y*e^2*n - 2*I*Y*e*n^3 + 2*I*Y*e*n + I ...
    *Y*n^3 - I*Y*n);;
end

