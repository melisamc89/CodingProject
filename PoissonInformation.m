% This function computes information by the suposition that the p(n) has an
% poisson distribution and that p(n|t) is algo a poisson distribution. 
% Where n0 is taken from data

function I=PoissonInformation(matrix,nt,aspect)

    [l m]=size(matrix);
    [r var]=MeanRate(matrix,nt);
   
    I0=zeros(m,1);
        for i=1:m
            nmin=0;
            n=nmin;
            nmax=100;
            hn=0;
            if var(i) ~= 0 && r(i) ~= 0
                while n < nmax
                    aux_nt=exp(-r(i))*(r(i)^n)/factorial(n);
                    aux_n=0;
                    for j=1:m
                        if var(j) ~= 0 && r(j) ~= 0
                        aux_n=aux_n+exp(-r(j))*r(j)^n/factorial(n);
                        end
                    end
                    aux_n=aux_n/m;
                    hn=hn+aux_nt*log(aux_nt/aux_n+eps);
                    n=n+1;
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
    