function Q = quadricov(X)
% Copyright: Anastasia Podosinnikova 2019
  X = X - repmat( mean(X,2), 1, size(X,2)); % zero-mean
  Q = quadricov_in(X);
  Q = (Q - diag(diag(Q)))' + Q;
end