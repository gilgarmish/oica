function ds_est = fourier_pca1(X, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X: p x N, p is dimension, N # samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright: Anastasia Podosinnikova 2019

  [p,~] = size(X);
  
  u1 = randn(p,1)/sqrt(p);
  u2 = randn(p,1)/sqrt(p);
  
  Q1 = genquadricov(X, u1);
  Q2 = genquadricov(X, u2);
  
  [W,S,~] = svd(Q1);
  [~,inds] = sort(diag(S),'descend');
  W = W(:,inds); % just in case 
  W = W(:,1:k);
  
  q1 = W'*Q1*W;
  q2 = W'*Q2*W;
  
  M = q1/q2;
  [V,~] = eig(M);
  
  C = W*V;
  
  % THEIR DEFINITION OF C IS NOT UNIQUE -> only up to sign, but that's not important
  for j = 1:k
    c = C(:,j);
    
    % computing the angle; there are always two solutions, we choose one
    a = real(c); b = imag(c);
    theta = atan( -( 2 * sum(a.*b) ) / sum( a.^2 - b.^2 ) )/2;
    while theta<0, theta = theta + pi; end
    while theta>2*pi, theta = theta - pi; end
    
    temp = real(exp(1i*theta)*c);
    C(:,j) = temp/norm(temp);
  end
  
  ds_est = zeros(p,k);
  for j = 1:k
    c = C(:,j);
    [v,s,~] = svd( reshape(c,p,p) );
    [~,inds] = sort(diag(s), 'descend');
    v = v(:,inds);
    ds_est(:,j) = v(:,1);
  end
  
end
