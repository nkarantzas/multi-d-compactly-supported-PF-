function Bmat = fminNEW(a,B1)
%% This function takes in an s-dimensional low-pass filter matrix 
% and #L pre-designed high-pass filters and returns a collection 
% of high-pass filter coefficients defining a Parseval Framelet 
% for L^2(R^s). The resulting matrix Bmat includes the pre-designed 
% filters but in a rescaled form that allows the existence of the PF. 

% 'a': A low-pass filter matrix with '>' coefficients adding up to 1.
% 'B1': An LxN matrix. Each row of B1 is an N - dim vector coming
% from a vectorized multi - dimensional high - pass filter array with
% numel = N. The elements of each row of B1 must add up to 0.

a = reshape(a,[1,numel(a)]); % turning a multi-d low-pass into a vector
c = sqrt(a);
if isempty(B1) == 1
    Q = c;
    Bmat = Pframe(Q);
    return
else
    D1 = B1*diag(1./c);
    Q = [c;D1];
end
nvars = size(B1,1);
LB = zeros(1,nvars); %specified lower bound
UB = ones(1,nvars); %specified upper bound
x0 = ones(1,nvars); %you may have to play around with the init - value
%x0 = (1/norm(Q(2:end,:),'fro'))*ones(1,nvars);
options = optimset('Largescale','on',...%'Display','iter',...
    'MaxIter',1000,'MaxFunEvals',5000,'TolX',1e-11,'TolFun',1e-11);
while max(svd(diag([1 x0])*Q)) > 1 + 10^(-6)
    x = fmincon(@(x)ObjFunction(x,Q),...
        x0,[],[],[],[],LB,UB,@(x)NonLinearCon(x,Q),options);
    x0 = x;
end
newQ = diag([1 x])*Q;
[U,S,V] = svd(newQ);
% If the constraint is violated from above even though an
% eps - machine - precision PF is obtained,
if S(1) < 1+10^(-6)
    S(1) = 1;
end
newQ = U*S*V';
Bmat = Pframe(newQ);
end