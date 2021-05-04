% main to generate all figures in one go

global talkative 
talkative = 0;

try
    lorenz84_cont_2
catch
    warning('lorenz48 continuation from Hopf to Hopf did not complete')
end
try
    figure_rychkov_saddle
catch
    warning('figure_rychkov continuation through saddle did not complete')
end
try
    lorenz84_validated_cont
catch
    warning('lorenz84_validated_cont continuation from Hopf did not complete')
end
try
    main_hyper
catch
    warning('main_hyper continuation from Hopf did not complete')
end
try
    main_hamiltonian
catch
    warning('main_hamiltonian continuation from Hopf did not complete')
end