function  h = ShannonEntropy(dat,M)
%%% function ht calculates the Shannpon entropy of a set of responses dat
%%% having M possible symbols
%%% INPUT:
%%% dat - set of all possible responses
%%% M - number of all possible symbols (for M=4, symbols: 0,1,2,3) 
%%% OUTPUT:
%%% h - Shannon entropy

L = size(dat,2);  %number of columns: L - length of the 'words'
n = size(dat,1);  %number of rows: n - number of 'words'
%there are M^L possible words

wi = 1 + dat*(M.^[0:L-1])'; %wi is a vector of all words in the set of responses
                        %dat (each possible word is a unique number)
count = histc(wi,[1:M^L+eps]); %count is the number of time each unique word  
                             %is present in the set of responses dat
p = (count'/n);               %p is the probability of each possible 'word' 
h = -sum(p.*log2(p+eps));     %h is the Shannon entropy
%Adding eps is used so that when p is 0 the log2 of 0+eps returns a value 
%(-52) which is zeroed latter by multiplying with p being 0. If eps is not
%used, it would return NaN whenever p is 0.
