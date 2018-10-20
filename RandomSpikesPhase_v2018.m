function new_pmatrix=RandomSpikesPhase_v2018(pmatrix)

    aux1=reshape(pmatrix,[size(pmatrix,1)*size(pmatrix,2)*size(pmatrix,3) 1]);
    phase_aux=aux1(find(aux));
    aux2=phase_aux(randperm(length(phase_aux)));
    aux3=aux1;
    aux3(find(aux1))=aux2;
    new_pmatrix=reshape(aux3,[size(pmatrix,1) size(pmatrix,2) size(pmatrix,3)]);
    
end