function B = Pframe(Q)
% Conditions on the input of this function
% 1. The first row of Q must be a 1xN unit vector 
% 2. The rest of the rows of Q are orthogonal to the first row of Q and
% they form an LxN matrix.
% 3. output 'B': vXN matrix whose rows contain the high - pass filter
% coefficients that define an s - dimensional PFramelet for L2(R^s).
N = size(Q,2);   
%% SVD on Q
[~, Sigma1, V] = svd(Q);
% Getting D2 
Sigma2 = sqrt(eye(N) - Sigma1'*Sigma1);
% Getting rid of near-zero elements
Sigma2 = Sigma2.*(abs(Sigma2) > 10^(-6));
D2 = Sigma2*V';
% matrix whose rows form a Parseval Frame for R^N
QD2 = [Q ; D2]; 
% getting rid of zero-rows (corresponding singular values = 1)
% look at how Sigma2 is defined
QD2(~any(QD2, 2), :) = [];
% matrix whose rows give the high-pass filter coefficients that define 
% a PFramelet for L2(R^s)
B = QD2 * diag(Q(1,:));
B = B.*(abs(B) > 10^(-6));
end
