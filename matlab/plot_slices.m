reaction = 1e3;

override = containers.Map({'integration'}, 4);

beam_Bending_Q4_4x1_Al = Solver("Beam_Bending_Q4_4x1_Al.txt");
beam_Bending_Q4_8x2_Al = Solver("Beam_Bending_Q4_8x2_Al.txt");
beam_Bending_Q4_16x4_Al = Solver("Beam_Bending_Q4_16x4_Al.txt");
beam_Bending_Q8_4x1_Al = Solver("Beam_Bending_Q8_4x1_Al.txt");
beam_Bending_Q9_4x1_Al = Solver("Beam_Bending_Q9_4x1_Al.txt");
beam_Bending_Q8_4x1_Al_under = Solver("Beam_Bending_Q8_4x1_Al.txt", override);
beam_Bending_Q9_4x1_Al_under = Solver("Beam_Bending_Q9_4x1_Al.txt", override);

names = ['Q4-4', 'Q4-8', 'Q4-16', 'Q8', 'Q9' ,'Q8 under-integrated', 'Q9 under-integrated'];
beams = [beam_Bending_Q4_4x1_Al; beam_Bending_Q4_8x2_Al; beam_Bending_Q4_16x4_Al; beam_Bending_Q8_4x1_Al; beam_Bending_Q9_4x1_Al; beam_Bending_Q8_4x1_Al_under; beam_Bending_Q9_4x1_Al_under];
%%
names = {'Q4-4', 'Q4-8', 'Q4-16', 'Q8', 'Q9' ,'Q8 under-integrated', 'Q9 under-integrated'}
figure
hold on
count = 1;
for name = names
    hold on
    beams(count).plot_axis();
    count = count + 1;
end
plot(0:0.1:12, 1e3*u((0:0.1:12) * 1e-3))
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
    beams(count).plot_midline_xx();
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
    beams(count).plot_midline_xy();
    count = count + 1;
end
% legend(names')
title("\sigma_{xy} along A -> A' for various meshes.")
xlabel("Distance from central axis (mm)")
ylabel("\sigma_{xy}")