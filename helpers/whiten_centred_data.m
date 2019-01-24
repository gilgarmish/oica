function [WX, W, Winv] = whiten_centred_data(X)
% Copyright: Anastasia Podosinnikova 2019
  if norm( mean(X,2) ) > 1e-8, error('data must be centred!'); end
  covx = cov(X');
  [V,D] = eig(covx);
  [d,inds] = sort( diag(D), 'descend' );
  V = V(:,inds);
  W = diag( d.^(-0.5) ) * V';
  Winv = V * diag( d.^(0.5) );
  WX = W*X;
end

