function new_matrix=RandomSpikesTimes_v2018(matrix)

    for i=1:size(matrix,1)
        aux=matrix(i,randperm(size(matrix,2)));
        new_matrix(i,:)=aux;
    end
    
end