function I=VonMisesInformation(matrix,nt,aspect)

    [n m l]=size(matrix);
    [circular_mean circular_variance]=CircularAnalysis(matrix,nt);

    I0=zeros(m,1);
        for i=1:m
            if circular_mean(i)~=0
            phase_min=-pi;
            phase_max=pi;
            dtheta=36*pi/180;
            phase=phase_min;
            hn=0;
                while phase < phase_max
                    aux_nt=exp(cos(phase-circular_mean(i))/circular_variance(i))/(2*pi*besselj(1,circular_variance(i)));
                    aux_n=0;
                    for j=1:m
                        aux_n=aux_n+exp(cos(phase-circular_mean(j))/circular_variance(j))/(2*pi*besselj(1,circular_variance(j)));
                    end
                    aux_n=aux_n/m;
                    hn=hn+aux_nt*log(aux_nt/aux_n+eps)*dtheta;
                    phase=phase+dtheta;
                end     
                I0(i)=hn;
            end    
        end
        
            if aspect==3
                I=I0(1)*0.15+I0(2)*0.7+I0(3)*0.15;
            else if aspect==5
                I=I0(1)*0.075+I0(2)*0.075+I0(3)*0.35+I0(4)*0.35+I0(5)*0.75+I0(6)*0.075;
            else
                I=mean(I0);
            end


end