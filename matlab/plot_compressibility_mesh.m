override = containers.Map({'nu'}, 0.4999);

[stress_44, disp_44] = Solver("Beam_Bending_Q4_16x4_PU.txt", override).contour_stress();
[stress_48, disp_48] = Solver("Beam_Bending_Q4_16x8_PU.txt", override1).contour_stress();
[stress_94, disp_94] = Solver("Beam_Bending_Q9_16x4_PU.txt", override2).contour_stress();

%%
figure
hold on
plot(0:0.1:12, stress_44(6, :, 1))
plot(0:0.1:12, stress_48(6, :, 1))
plot(0:0.1:12, stress_94(6, :, 1))

figure
hold on
plot(0:0.1:1, disp_44(:, 61, 1))
plot(0:0.1:1, disp_48(:, 61, 1))
plot(0:0.1:1, disp_94(:, 61, 1))