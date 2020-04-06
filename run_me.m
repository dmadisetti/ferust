%% This script loads and solves a static finite element problem.

% You may need to add the following directories to your path.
%addpath('examples', 'matlab');
%% Solve for displacements
% Sample case
biaxial_Q4_2x2 = Solver("Biaxial_Q4_2x2.txt");

% AL Beams
beam_Bending_Q4_4x1_Al = Solver("Beam_Bending_Q4_4x1_Al.txt");
beam_Bending_Q4_8x2_Al = Solver("Beam_Bending_Q4_8x2_Al.txt");
beam_Bending_Q4_16x4_Al = Solver("Beam_Bending_Q4_16x4_Al.txt");
beam_Bending_Q8_4x1_Al = Solver("Beam_Bending_Q8_4x1_Al.txt");
beam_Bending_Q9_4x1_Al = Solver("Beam_Bending_Q9_4x1_Al.txt");

% PU Beams
beam_Bending_Q9_16x4_PU = Solver("Beam_Bending_Q9_16x4_PU.txt");
beam_Bending_Q4_16x4_PU = Solver("Beam_Bending_Q4_16x4_PU.txt");
beam_Bending_Q4_16x8_PU = Solver("Beam_Bending_Q4_16x8_PU.txt");