function [v, H] = householder_vector(x)
    % HOUSEHOLDER_VECTOR computes the Householder vector v and matrix H
    % to reflect x onto the first basis vector.
    
    n = length(x);
    e1 = eye(n, 1);
    
    % Compute the Householder vector v
    v = x + sign(x(1)) * norm(x) * e1;
    v = v / norm(v); % Normalize for stability
    
    % Form the Householder matrix H = I - 2*v*v'
    H = eye(n) - 2 * (v * v');
end