%Given data in form of matrix (where rows are trials and columns are
%binnings of time), this function calculetes the probabilitity of the
%events assuming Poisson dsitribution. It computes the mean in each colum
%and aproximates the distibution of events given the bin of time
%p(n|t)= Poisson distribucion. The rate of this distribution is taken as
%the mean of the data in each column.

function [pn pnt pt]=PoissonProbability(matrix1,matrix2)
    
    nmax=50;
    matrix=[matrix1,matrix2];
    [ntrials ns]=size(matrix);
    r=mean(matrix,1);
    
    pnt=zeros(nmax,ns);
    pt=zeros(1,ns);
    
    contador=0;
    for i=1:ns
        x=find(matrix(:,i));
        if length(x)>5 && r(i) ~= 0
            contador=contador+1;
            pt(1,i)=1;
            for n=0:nmax-1
                pnt(n+1,i)=exp(-r(i))*(r(i)^n)/factorial(n);
            end
            suma=sum(pnt(:,i));
            pnt(:,i)=pnt(:,i)/suma;
        end
    end

    if contador ~=0
        pt=pt/contador;
    end

    pn=zeros(1,nmax);
    for n=0:nmax-1
        if contador ~=0
            pn(n+1)=sum(pnt(n+1,:))/contador;
        end
    end
    suma=sum(pn);
    if suma~=0
        pn=pn/suma;
    end


end