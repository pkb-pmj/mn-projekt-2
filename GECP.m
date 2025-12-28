function [X] = GECP(A,Y)
%GECP Gaussian Elimination with Complete Pivoting
%   Performs gaussian elimination on matrix [A Y] to get [I X]
arguments (Input)
    A % left matrix (square)
    Y % right vector/matrix (must have the same number of rows as A)
end

arguments (Output)
    X % output matrix (same dimensions as Y)
end

n = length(Y);
M = [A Y];
cols = 1:n;
for k = 1:n
    % find indices of max element
    [~, argmax] = max(M(k:n,k:n), [], "all", ComparisonMethod="abs");
    [i, j] = ind2sub(size(M(k:n,k:n)), argmax);
    i = i+k-1;
    j = j+k-1;
    % swap rows
    M([k, i],:) = M([i, k],:);
    % swap columns
    M(:,[k,j]) = M(:,[j,k]);
    cols([k,j]) = cols([j,k]);
    % zero elements in k-th column
    for i = k+1:n
        M(i,k:end) = M(i,k:end) - M(k,k:end) * M(i,k) / M(k,k);
    end
    % zero elements below explicitly to be sure they are exactly 0
    M(k+1:n,k) = 0;
end
% normalize pivot elements
for k = 1:n
    M(k,k:end) = M(k,k:end) / M(k,k);
end
% back substitution
for k = n:-1:1
    for i = 1:k-1
        M(i,k:end) = M(i,k:end) - M(k,k:end) * M(i,k); % / M(k,k) == 1
    end
    % zero elements above explicitly to be sure they are exactly 0
    M(1:k-1,k) = 0;
end
% extract X
X = M(:,n+1:end);
X(cols,:) = X(:,:);

end