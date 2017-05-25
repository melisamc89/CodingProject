function [p_thetat p_theta pt]=ConditionalVonMisesProbability(n_matrix,matrixa,matrixb)

    nbins=100;
    size(matrixa);
    [ntrials ns]=size(n_matrix);
    r=mean(n_matrix,1);
    variance=var(n_matrix,1);
    auxvar=min(variance(find(variance)));
    
    r_mean=mean(r);
    variance1=var(reshape(n_matrix,ntrials*ns,1));
    nmin=r_mean-3*sqrt(variance1);
    nmax=r_mean+3*sqrt(variance1);
    dn=1.0*(nmax-nmin)/nbins;

    pmin=0;
    pmax=2*pi;
    num=40;
    matrix=[matrixa,matrixb];
    covariance=NandFaseCorrelation(n_matrix,matrix);
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
   
    kappa=zeros(1,length(circular_variance));
    for i=1:length(circular_variance)
        kappa(i)=(-1/(2*log(1-circular_variance(i))));
    end
    %kmax=max(kappa);
    %numero_puntos=floor(num*pi*sqrt(kmax));

    ns=ns-1;
    p_thetat=zeros(nbins,num,ns);
    pt=ones(1,ns);
    contador=ns;
    for i=1:ns
        x=find(matrix(:,i));
        if length(x)<=5
            pt(i)=0;
            contador=contador-1;
        end
        if length(x)>5  && r(i)~=0
        n=nmin;
        it=1;
        while n<nmax && it<nbins+1
            circular_mean(i)=circular_mean(i)+covariance(i)*(n-r(i))/variance(i);
            j=1;
            o=pmin;
            while o<pmax && j<num 
                p_thetat(it,j,i)=exp(cos(o-circular_mean(i))*kappa(i))/(2*pi*besseli(0,kappa(i)));
                o=j*(2*pi)/num;
                j=j+1;
            end
            n=n+dn;
            it=it+1;
        end
        end
    end
    pt=pt/contador;
    suma=sum(p_thetat,2);
    for t=1:ns
        for n=1:nbins
            p_thetat(n,:,t)=p_thetat(n,:,t)/suma(n,t);
        end
    end
    
    p_theta=zeros(nbins,num);
    for n=1:nbins
        for j=1:num
            for t=1:ns
                if pt(t)~=0
                p_theta(n,j)=p_theta(n,j)+p_thetat(n,j,t)*pt(t);
                end
            end
        end
    end
end