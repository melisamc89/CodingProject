function [p_theta p_theta_n pn vec phase_mean kappa]=GaussianVonMisesDistribution2(n_matrix,matrixa,matrixb)

    [p_theta_n pn p_theta vec phase_mean kappa]=ConditionalVonMisesProbability2(n_matrix,matrixa,matrixb);
    
end