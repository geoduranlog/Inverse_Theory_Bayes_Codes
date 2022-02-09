%% -- Misfit ---

function mis = misfit(syndata,data,sigma)   

%Convert matrices into vectors
data= data(:);
syndata=syndata(:);

mis =-(norm( data - syndata ) ^ 2) / (2*sigma^2);

% %-Small test-
% A=rand(2); B=rand(2)
% A=A(:); B=B(:)
% (norm( A - B ) ^ 2)
% ( A - B)'*( A - B)


end