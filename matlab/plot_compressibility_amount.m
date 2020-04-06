reaction = 1e3;

override = containers.Map({'nu'}, 0.4);
override1 = containers.Map({'nu'}, 0.49);
override2 = containers.Map({'nu'}, 0.4999);

[stress_4, disp_4] = Solver("Beam_Bending_Q4_16x4_PU.txt", override).contour_stress();
[stress_49, disp_49] = Solver("Beam_Bending_Q4_16x4_PU.txt", override1).contour_stress();
[stress_4999, disp_4999] = Solver("Beam_Bending_Q4_16x4_PU.txt", override2).contour_stress();

%%
figure
hold on
plot((0:0.1:1) - 0.5, stress_4(61, :, 1))
plot((0:0.1:1) - 0.5, stress_49(61, :, 1))
plot((0:0.1:1) - 0.5, stress_4999(61, :, 1))
title("\sigma_{xx} along A -> A' for various \nu values.")
ylabel("\sigma_{xx}")
xlabel("Distance from central axis (mm)")
legend("\nu=0.4", "\nu=0.49", "\nu=0.4999")
%%
figure
hold on
plot(0:0.1:12, disp_4(:, 6, 2) -0.5)
plot(0:0.1:12, disp_49(:, 6, 2) -0.5)
plot(0:0.1:12, disp_4999(:, 6, 2) -0.5)
title("Displacement along central axis for various \nu values.")
xlabel("Distance from wall (mm)")
ylabel("Displacement (mm)")
legend("\nu=0.4", "\nu=0.49", "\nu=0.4999")