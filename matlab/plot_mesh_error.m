reaction = 1e3;

override = containers.Map({'integration'}, 4);

beam_Bending_Q4_4x1_Al = Solver("Beam_Bending_Q4_4x1_Al.txt");
beam_Bending_Q4_8x2_Al = Solver("Beam_Bending_Q4_8x2_Al.txt");
beam_Bending_Q4_16x4_Al = Solver("Beam_Bending_Q4_16x4_Al.txt");
beam_Bending_Q8_4x1_Al = Solver("Beam_Bending_Q8_4x1_Al.txt");
beam_Bending_Q9_4x1_Al = Solver("Beam_Bending_Q9_4x1_Al.txt");
beam_Bending_Q8_4x1_Al_under = Solver("Beam_Bending_Q8_4x1_Al.txt", override);
beam_Bending_Q9_4x1_Al_under = Solver("Beam_Bending_Q9_4x1_Al.txt", override);


log_lengths = [log(12/4), log(12/8) log(12/16)];
figure
plot(log_lengths, [...
        log(abs(beam_Bending_Q4_4x1_Al.reaction_sum(2) - reaction)/abs(reaction)), ...
        log(abs(beam_Bending_Q4_8x2_Al.reaction_sum(2) - reaction)/abs(reaction)), ...
        log(abs(beam_Bending_Q4_16x4_Al.reaction_sum(2) - reaction)/abs(reaction)), ...
    ], '.-', 'MarkerSize', 20);
hold on
scatter([log(12/4)],[...
        log(abs(beam_Bending_Q8_4x1_Al.reaction_sum(2) - reaction)/abs(reaction)), ...
    ], 30, 	[0.8500, 0.3250, 0.0980], 'filled');
scatter([log(12/4)],[...
        log(abs(beam_Bending_Q9_4x1_Al.reaction_sum(2) - reaction)/abs(reaction)), ...
    ], 30, [0.4660, 0.6740, 0.1880], 'filled');
legend()