reaction = 1e3;

override = containers.Map({'nu'}, 0.4);
override1 = containers.Map({'nu'}, 0.49);
override2 = containers.Map({'nu'}, 0.4999);

nu4 = Solver("Beam_Bending_Q4_16x4_PU.txt", override);
nu49 = Solver("Beam_Bending_Q4_16x4_PU.txt", override1);
nu4999 = Solver("Beam_Bending_Q4_16x4_PU.txt", override2);

%%
figure
hold on
nu4.plot_midline_xx()
nu49.plot_midline_xx()
nu4999.plot_midline_xx()
title("\sigma_{xx} along A -> A' for various \nu values.")
ylabel("\sigma_{xx}")
xlabel("Distance from central axis (mm)")
legend("\nu=0.4", "\nu=0.49", "\nu=0.4999")
%%
figure
hold on
nu4.plot_centerline()
nu49.plot_centerline()
nu4999.plot_centerline()
title("Displacement along central axis for various \nu values.")
xlabel("Distance from wall (mm)")
ylabel("Distance from center line (mm)")
legend("\nu=0.4", "\nu=0.49", "\nu=0.4999")