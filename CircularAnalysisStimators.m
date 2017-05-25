function [theta variance]=CircularAnalysisStimators(matrixa,matrixb)
    
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
    [theta variance]=CircularAnalysis(matrix);
    
  
end