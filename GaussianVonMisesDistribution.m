%%Revisar esta funcion

function [p_thetat p_theta pt]=GaussianVonMisesDistribution(n_matrix,matrixa,matrixb)
    

    [p_thetat p_theta pt]=ConditionalVonMisesProbability(n_matrix,matrixa,matrixb);
    
    for i=1:length(pn)
        pntheta(i,:)=pn(i)*p_theta_n(i,:);
    end    
    pntheta=pntheta/(sum(sum(pntheta)));
    
    [nbins,pbins,timebins]=size(p_theta_nt);
    for t=1:timebins
        for i=1:nbins
            for j=1:pbins
                pntheta_t(i,j,t)=pn(i)*p_theta_nt(i,j,t);
            end
        end
    end
    
    for t=1:timebins
        pntheta_t(:,:,t)=pntheta_t(:,:,t)/sum(sum(pntheta_t(:,:,t)));
    end
 
    pnthetat=zeros(nbins,pbins);
    for t=1:timebins
        pnthetat=pnthetat+pt(t)*pntheta_t(:,:,t);
    end
    pnthetat=pnthetat/sum(sum(pnthetat));
    
end