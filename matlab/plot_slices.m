reaction = 1e3;

override = containers.Map({'integration'}, 4);

[stress_Q4_4x1_Al, displacement_Q4_4x1_Al] = Solver("Beam_Bending_Q4_4x1_Al.txt").contour_stress();
[stress_Q4_8x2_Al, displacement_Q4_8x2_Al] = Solver("Beam_Bending_Q4_8x2_Al.txt").contour_stress();
[stress_Q4_16x4_Al, displacement_Q4_16x4_Al] = Solver("Beam_Bending_Q4_16x4_Al.txt").contour_stress();
[stress_Q8_4x1_Al, displacement_Q8_4x1_Al] = Solver("Beam_Bending_Q8_4x1_Al.txt").contour_stress();
[stress_Q9_4x1_Al, displacement_Q9_4x1_Al] = Solver("Beam_Bending_Q9_4x1_Al.txt").contour_stress();
[stress_Q8_4x1_Al_under, displacement_Q8_4x1_Al_under] = Solver("Beam_Bending_Q8_4x1_Al.txt", override).contour_stress();
[stress_Q9_4x1_Al_under, displacement_Q9_4x1_Al_under] = Solver("Beam_Bending_Q9_4x1_Al.txt", override).contour_stress();

names = ['Q4-4', 'Q4-8', 'Q4-16', 'Q8', 'Q9' ,'Q8 under-integrated', 'Q9 under-integrated'];
stresses = [stress_Q4_4x1_Al; stress_Q4_8x2_Al; stress_Q4_16x4_Al; stress_Q8_4x1_Al; stress_Q9_4x1_Al; stress_Q8_4x1_Al_under; stress_Q9_4x1_Al_under];
displacement = [displacement_Q4_4x1_Al displacement_Q4_8x2_Al displacement_Q4_16x4_Al displacement_Q8_4x1_Al displacement_Q9_4x1_Al displacement_Q8_4x1_Al_under displacement_Q9_4x1_Al_under];
%%
names = {'Q4-4', 'Q4-8', 'Q4-16', 'Q8', 'Q9' ,'Q8 under-integrated', 'Q9 under-integrated'}
figure
hold on
count = 1;
for name = names
    hold on
    plot(0:0.1:12, displacement(:, 6 + 11 * (count - 1), 2) - 0.5)
    count = count + 1;
end
plot(0:0.1:12, (0:0.1:12)*.5/12)
names{8} = 'Beam Theory';
legend(names')
title("Displacement along central axis for various meshes.")
xlabel("Distance from wall (mm)")
ylabel("Displacement (mm)")
%%
figure
subplot(2, 1, 1)
names = {'Q4-4', 'Q4-8', 'Q4-16', 'Q8', 'Q9' ,'Q8 under-integrated', 'Q9 under-integrated'}
hold on
count = 1;
for name = names
    hold on
    plot((0:0.1:1) - 0.5, stresses(61 + 121 * (count - 1), :, 1))
    count = count + 1;
end
legend(names')
title("\sigma_{xx} along A -> A' for various meshes.")
ylabel("\sigma_{xx}")
hold on
%%
subplot(2, 1, 2)
names = {'Q4-4', 'Q4-8', 'Q4-16', 'Q8', 'Q9' ,'Q8 under-integrated', 'Q9 under-integrated'}
count = 1;
for name = names
    hold on
    61 + 121 * (count - 1)
    plot((0:0.1:1) - 0.5, stresses(61 + 121 * (count - 1), :, 2))
    count = count + 1;
end
% legend(names')
title("\sigma_{xy} along A -> A' for various meshes.")
xlabel("Distance from central axis (mm)")
ylabel("\sigma_{xy}")