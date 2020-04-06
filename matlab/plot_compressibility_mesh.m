% override = containers.Map({'nu'}, 0.4999);
% 
% [stress_44, disp_44] = Solver("Beam_Bending_Q4_16x4_PU.txt", override).contour_stress();
% [stress_48, disp_48] = Solver("Beam_Bending_Q4_16x8_PU.txt", override1).contour_stress();
% [stress_94, disp_94] = Solver("Beam_Bending_Q9_16x4_PU.txt", override2).contour_stress();

%%
figure
hold on
plot((0:0.1:1) - 0.5, stress_44(61, :, 1))
plot((0:0.1:1) - 0.5, stress_48(61, :, 1))
plot((0:0.1:1) - 0.5, stress_94(61, :, 1))
title("\sigma_{xx} along A -> A' for various PU meshes.")
ylabel("\sigma_{xx}")
xlabel("Distance from central axis (mm)")
legend("Q4-4", "Q4-8", "Q9")

%%
figure
hold on
plot(0:0.1:12, disp_44(:, 6, 2) - 0.5)
plot(0:0.1:12, disp_48(:, 6, 2) - 0.5)
plot(0:0.1:12, disp_94(:, 6, 2) - 0.5)
title("Displacement along central axis for various PU meshes.")
xlabel("Distance from wall (mm)")
ylabel("Displacement (mm)")
legend("Q4-4", "Q4-8", "Q9")