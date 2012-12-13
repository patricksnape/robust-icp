% Least squares normal estimation from point clouds using PCA
%
% H. Hoppe, T. DeRose, T. Duchamp, J. McDonald, and W. Stuetzle. 
% Surface reconstruction from unorganized points. 
% In Proceedings of ACM Siggraph, pages 71:78, 1992.
%
% p should be a matrix containing the horizontally concatenated column
% vectors with points. k is a scalar indicating how many neighbors the
% normal estimation is based upon.
%
% Jakob Wilm 2010

function n = lsqnormest(p, k)
m = size(p, 2);
n = zeros(3, m);

neighbors = transpose(knnsearch(transpose(p), transpose(p), 'k', k + 1));

for i = 1:m
    x = p(:, neighbors(2:end, i));
    p_bar = (1 / k) * sum(x, 2);
    
    % spd matrix P
    P = (x - repmat(p_bar, 1, k)) * transpose(x - repmat(p_bar, 1, k)); 
    
    [V, D] = eig(P);
    
    % choses the smallest eigenvalue
    [~, idx] = min(diag(D)); 
    
    % returns the corresponding eigenvector
    n(:, i) = V(:, idx);   
end