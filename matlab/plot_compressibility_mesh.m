override = containers.Map({'nu'}, 0.4999);
 
PU_44 = Solver("Beam_Bending_Q4_16x4_PU.txt", override);
PU_48 = Solver("Beam_Bending_Q4_16x8_PU.txt", override);
PU_94 = Solver("Beam_Bending_Q9_16x4_PU.txt", override);

%%
figure
subplot(1,2,1)
hold on
PU_44.plot_midline_xx();
PU_48.plot_midline_xx();
PU_94.plot_midline_xx();
title("\sigma_{xx} along A -> A' for various PU meshes.")
ylabel("\sigma_{xx}")
xlabel("Distance from central axis (mm)")
legend("Q4-4", "Q4-8", "Q9")
%%
subplot(1,2,2)
hold on
PU_44.plot_midline_xy();
PU_48.plot_midline_xy();
PU_94.plot_midline_xy();
title("\sigma_{xy} along A -> A' for various PU meshes.")
ylabel("\sigma_{xy}")
xlabel("Distance from central axis (mm)")
legend("Q4-4", "Q4-8", "Q9")


%%
figure
hold on
PU_44.plot_centerline();
PU_48.plot_centerline();
PU_94.plot_centerline();
title("Displacement along central axis for various PU meshes.")
xlabel("Distance from wall (mm)")
ylabel("Displacement (mm)")
legend("Q4-4", "Q4-8", "Q9")