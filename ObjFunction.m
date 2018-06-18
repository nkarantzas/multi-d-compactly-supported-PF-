 function f = ObjFunction(lambda,Q)
% lambda is an L - dimensional vector with first component equal to 1 
% and the rest unknown. Q is an LxN given matrix. We're multiplying
% each row of Q by the corresponding component of lambda. The goal is 
% to maximize the singular values of the resulting LxN matrix 
% under the constraint that the largest singular value 
% is less than or equal to  1.

% defining the objective function.
D = diag([1 lambda]);
M = D*Q;
f = -trace(M'*M);
end

