function dkl=KLDivergence_v2018(p,q)
    dkl=sum((p+eps)'.*log((p+eps)'./(q+eps)'));
end