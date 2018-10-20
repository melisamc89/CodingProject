function [p_theta p_thetat pt]=VonMisesProbability_v2018(matrix_firing,matrix1,num)

    pmin=0;
    pmax=2*pi;
    [n m l]=size(matrix1);
    %Para pasar de la matrix 3D a una 2D
    matrix=zeros(n*l,m);
    
    for j=1:m
         aux1=matrix(1,j,:);
         for i=2:n
              aux1=cat(3,aux1,matrix1(i,j,:));
         end
        matrix(1:length(aux1),j)=aux1(1,1,:);
    end
    
    [ntrials ns]=size(matrix);
    
    [circular_mean,circular_variance]=CircularAnalysis(matrix);
   
    kappa=zeros(1,length(circular_variance));
    for i=1:length(circular_variance)
        kappa(i)=(-1/(2*log(1-circular_variance(i))));
    end
    circ_var_min=min(circular_variance);
    kmin=(-1/(2*log(1-circ_var_min)));

    
    ns=ns-1;
    p_theta=zeros(num,1);
    p_thetat=zeros(num,ns);
    pt=ones(1,ns);
    time_counter=ns;
    
    for i=1:ns
        x=find(matrix(:,i));
        x1=find(matrix_firing(:,i));
        if length(x)<=20 || circular_mean(i)== 0 || length(x1)<5
            kappa(i)=kmin;
            pt(i)=0;
            time_counter=time_counter-1;
        end
        n=1;
        o=pmin;
        while o<pmax && n<=num 
            p_thetat(n,i)=exp(cos(o-circular_mean(i))*kappa(i))/(2*pi*besseli(0,kappa(i)));
            o=n*(2*pi)/num;
            n=n+1;
        end
    end
    
    if time_counter
        pt=pt/time_counter;
    end
    suma=sum(p_thetat,1);
    for t=1:ns
        p_thetat(:,t)=p_thetat(:,t)/suma(t);
    end

    p_theta=p_thetat*pt';
    suma=sum(p_theta);
    p_theta=p_theta/suma;
    
end