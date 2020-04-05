reaction = 1e3;

override = containers.Map({'integration'}, 4);

[stress_Q4_4x1_Al, displacement_Q4_4x1_Al] = Solver("Beam_Bending_Q4_4x1_Al.txt").contour_stress();
[stress_Q4_8x2_Al, displacement_Q4_8x2_Al] = Solver("Beam_Bending_Q4_8x2_Al.txt").contour_stress();
[stress_Q4_16x4_Al, displacement_Q4_16x4_Al] = Solver("Beam_Bending_Q4_16x4_Al.txt").contour_stress();
[stress_Q8_4x1_Al, displacement_Q8_4x1_Al] = Solver("Beam_Bending_Q8_4x1_Al.txt").contour_stress();
[stress_Q9_4x1_Al, displacement_Q9_4x1_Al] = Solver("Beam_Bending_Q9_4x1_Al.txt").contour_stress();
[stress_Q8_4x1_Al_under, displacement_Q8_4x1_Al_under] = Solver("Beam_Bending_Q8_4x1_Al.txt", override).contour_stress();
[stress_Q9_4x1_Al_under, displacement_Q9_4x1_Al_under] = Solver("Beam_Bending_Q9_4x1_Al.txt", override).contour_stress();

%%
figure
plot()
beam_Bending_Q4_4x1_Al.contour_stress()