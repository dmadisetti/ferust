clc; clear;
H = 1e-3;
B = 1e-3;
rho = 2710;
L = 12e-3;
E = 70e9;
P=0.5e-3;
g = 9.82;

I = B * H^3/12;

q = rho * g * H * B;

% From M(L) = 0
c2 = @(c1) -q*L^2/(2*E*I) - c1*L

% From u(L) = -P
c1 = -(5*q*L^4 - 24*E*I*P)/(8*E*I*L^3)

c2 = c2(c1)
c2 =(q*L^4 - 24*E*I*P)/(8*E*I*L^2)
u = @(x) q*x^4/(24*E*I) + c1 * x^3/6 + c2 * x^2/2

shear = @(x) q * x + c1*E*I 
M = @(x) q*x^2/(2*E*I) + c1 * x + c2
u = @(x) q*x.^4/(24*E*I) + c1 * x.^3/6 + c2 * x.^2/2

%assert(M(L) == 0)
%assert(u(L) == -P)

r = shear(0)

% Consider at point A
sigma_xy = -shear(L/2)
sigma_xx = -(M(L/2) * (H/2)/I)
u(L)
%%
clc; clear;
syms c1 c2 q L E I H P q
M = @(x) q*x^2/(2*E*I) + c1 * x + c2
u = @(x) q*x^4/(24*E*I) + c1 * x^3/6 + c2 * x^2/2

[M(L), u(L) + P]
[c1, c2] = solve([M(L), u(L) + P], [c1 c2])
shear = @(x) q * x + c1*E*I 
shear(0)