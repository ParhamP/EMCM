function W = exp_weight(X)
    norm = X(:, 1) + 0.0001;
    numer = exp(bsxfun(@rdivide, -X, norm));
    denom = sum(numer, 2);
    W = bsxfun(@rdivide, numer, denom);
end