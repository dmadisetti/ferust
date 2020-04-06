%% This script loads and solves a static finite element problem.

% You may need to add the following directories to your path.
%addpath('examples', 'matlab');
%% Solve for displacements
% Sample case
biaxial_Q4_2x2 = Solver("Biaxial_Q4_2x2.txt");

% AL Beams
beam_Bending_Q4_4x1_Al = Solver("Beam_Bending_Q4_4x1_Al.txt");
beam_Bending_Q4_4x1_Al.plot_centerline()
beam_Bending_Q4_8x2_Al = Solver("Beam_Bending_Q4_8x2_Al.txt");
beam_Bending_Q4_16x4_Al = Solver("Beam_Bending_Q4_16x4_Al.txt");
beam_Bending_Q8_4x1_Al = Solver("Beam_Bending_Q8_4x1_Al.txt");
beam_Bending_Q9_4x1_Al = Solver("Beam_Bending_Q9_4x1_Al.txt");

% PU Beams
beam_Bending_Q9_16x4_PU = Solver("Beam_Bending_Q9_16x4_PU.txt");
beam_Bending_Q4_16x4_PU = Solver("Beam_Bending_Q4_16x4_PU.txt");
beam_Bending_Q4_16x8_PU = Solver("Beam_Bending_Q4_16x8_PU.txt");

%% Examine results
figure
subplot(3,3,1)
beam_Bending_Q4_4x1_Al.plot_nodes_displaced()
xlim([0, 12.5])
ylim([-1, 1.5])

subplot(3,3,2)
beam_Bending_Q4_8x2_Al.plot_nodes_displaced()
xlim([0, 12.5])
ylim([-1, 1.5])

subplot(3,3,3)
beam_Bending_Q4_16x4_Al.plot_nodes_displaced()
xlim([0, 12.5])
ylim([-1, 1.5])

subplot(3,3,4)
beam_Bending_Q8_4x1_Al.plot_nodes_displaced()
xlim([0, 12.5])
ylim([-1, 1.5])

subplot(3,3,5)
biaxial_Q4_2x2.plot_nodes_displaced()
xlim([-1, 3])
ylim([-1, 3])

subplot(3,3,6)
beam_Bending_Q9_4x1_Al.plot_nodes_displaced()
xlim([0, 12.5])
ylim([-1, 1.5])

subplot(3,3,7)
beam_Bending_Q9_16x4_PU.plot_nodes_displaced()
xlim([0, 12.5])
ylim([-1, 1.5])

subplot(3,3,8)
beam_Bending_Q4_16x4_PU.plot_nodes_displaced()
xlim([0, 12.5])
ylim([-1, 1.5])

subplot(3,3,9)
beam_Bending_Q4_16x8_PU.plot_nodes_displaced()
xlim([0, 12.5])
ylim([-1, 1.5])