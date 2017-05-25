function [nr sigmar nl sigmal]=GaussianAnalysisStimators(matrix_1,matrix_2)

    nr=mean(matrix_1);
    nl=mean(matrix_2);
    sigmar=var(matrix_1);
    sigmal=var(matrix_2);
end
