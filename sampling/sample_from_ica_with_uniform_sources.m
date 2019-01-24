function [X, Alpha] = sample_from_ica_with_uniform_sources(ds, N)
% ds : mixing matrix
% N  : number of observations, n = 1, 2, ..., N
%
% model: x[n] = D*alpha[n]
%        alpha[n] ~ uniform( 0.5,1 )

% Copyright: Anastasia Podosinnikova 2019

  if ~( nargin==2 )
    error('sample data: wrong number of inputs')
  end

  [p,k] = size(ds);
  X = zeros(p,N);
  Alpha = zeros(k,N);
  
  % sample in batches (faster)
  n = 1000;
  times = floor(N/n);
  rest = mod(N,n);
  
  for i=1:times
    inds = (i-1)*n + 1:i*n;
    [x, alpha] = sample_batch(n);
    X(:,inds) = x;
    Alpha(:,inds) = alpha;
  end
  
  inds = times*n + 1 : times*n + rest;
  [x, alpha] = sample_batch(rest);
  X(:,inds) = x;
  Alpha(:,inds) = alpha;
  
  function [x, alpha] = sample_batch(n)
    alpha = rand(k,n).*abs( randn(k,n) );
    x = sparse( ds*alpha );
%     alpha = 0.5 + 0.5*rand(k,n);
%     x = sparse( ds*alpha );
    % WARNING: IS DATA (ALWAYS) SPARSE?
  end
  
end


