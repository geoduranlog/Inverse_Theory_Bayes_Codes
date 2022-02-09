function [RandV] = RandHistValues(V,X)
%From a histogram, generate random values following such
%distribution(histogram). The values are chosen randomly but the frequency
%a given value appears remains constant.

%X: vector containing the values
%V:vector containing the frequencies of the different values (patterns).
%The total number of patterns is length(X)


% a simple but effective run-length decoding scheme
clear ix ; 
ix([1 cumsum(V(1:end))+1]) = 1 ;
ix = cumsum(ix(1:end-1)) ;
% retrieve the values
X2 = X(V>0) ;
V = X2(ix) ; % all length(X) values
RandV = V(randperm(numel(V))); % in random order

end