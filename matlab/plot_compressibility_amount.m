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
plot(0:0.1:12, stress_4(6, :, 1))
plot(0:0.1:12, stress_49(6, :, 1))
plot(0:0.1:12, stress_4999(6, :, 1))

figure
hold on
plot(0:0.1:1, disp_4(:, 61, 1))
plot(0:0.1:1, disp_49(:, 61, 1))
plot(0:0.1:1, disp_4999(:, 61, 1))