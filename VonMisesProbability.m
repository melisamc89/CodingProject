function [p_theta p_thetat pt]=VonMisesProbability(matrixa,matrixb,phases)
    
    matrix=[matrixa,matrixb];
    [ntrials ns]=size(matrix);
    
    %Para pasar de la matrix 3D a una 2D
    [n m l]=size(matrixa);
    matrix1=zeros(n*l,m);
    matrix2=zeros(n*l,m);
    
    for j=1:m
         aux1=matrixa(1,j,:);
         aux2=matrixb(1,j,:);
         for i=2:n
              aux1=cat(3,aux1,matrixa(i,j,:));
              aux2=cat(3,aux2,matrixb(i,j,:));
         end
        matrix1(1:length(aux1),j)=aux1(1,1,:);
        matrix2(1:length(aux2),j)=aux2(1,1,:);
    end
    
    matrix=[matrix1,matrix2];
    [ntrials ns]=size(matrix);
    
    [circular_mean,circular_variance]=CircularAnalysis(matrix);
   
    for i=1:length(circular_variance)
        kappa(i)=(-1/(2*log(1-circular_variance(i))));
    end
    kmax=max(kappa);
    numero_puntos=floor(20*pi*sqrt(kmax));
    dtheta=2*pi/numero_puntos;
    
    p_thetat=zeros(phases,ns);
    pt=ones(1,ns);
    
    contador=0;
    for i=1:ns
        x=find(matrix(:,i));
                    contador=contador+1;

        if length(x)>3 && circular_variance(i)~=0
            pt(1,i)=1;
            for n=1:phases
                o=n*(2*pi)/phases;
                p_thetat(n,i)=0;
                while (o < (n+1)*2*pi/phases)
                    p_thetat(n,i)=p_thetat(n,i)+(exp(cos(o-circular_mean(i))*kappa(i))/(2*pi*besseli(0,kappa(i))))*dtheta;
                    o=o+dtheta;
                end
            end
            suma=sum(p_thetat(:,i));
            p_thetat(:,i)=p_thetat(:,i)/suma;
        end
    end
    
    if contador~=0
        pt=pt/ns;
    end

    p_theta=zeros(1,phases);
    for n=1:phases
        if contador~=0
        p_theta(n)=sum(p_thetat(n,:))/contador;
        end
    end
    
    suma=sum(p_theta);
    if suma ~=0
        p_theta=p_theta/suma;
    end
end