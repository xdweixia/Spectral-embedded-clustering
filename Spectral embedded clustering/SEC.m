function [F, W, b] = SEC(X, L, para)
% X: each column is a data point
% L: Laplacian matrix

[dim,n] = size(X);
Lc = eye(n) - 1/n*ones(n,n);

if dim < n
    Xprod = X*Lc*X';
    G = Xprod + 1/para.gamma*eye(dim);
    A = inv(G)*X*Lc;
    clear G Xprod;
    M = L + para.mu*para.gamma*(Lc - Lc*X'*A);
else
    Xprod = Lc*X'*X*Lc;
    G = Xprod + 1/para.gamma*eye(n);
    clear Xprod;
    Gi = inv(G);
    clear G;
    M = L + para.mu*(Lc * Gi);
end

clear Lc;
M = (M+M')/2;
[v,d] = eig(M);
[d, idx] = sort(diag(d));
if para.mu == 0
    F = v(:, idx(2:para.nc+1));
else
    F = v(:,idx(1:para.nc));
end
clear M v;

if nargout > 2
% calculate W and b for out-of-sample extension to clustering
    if dim < n
        W = A*F;
        b = mean(F',2) - mean(W'*X,2);
    else
        Lc = eye(n) - 1/n*ones(n,n);
        W = X*(Lc*Gi*F);
        b = mean(F',2) - mean(W'*X,2);
    end
end






