function ds = fourier_pca2(X, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X: p x N, p is dimension, N # samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright: Anastasia Podosinnikova 2019

  [p,~] = size(X);
  
  u = randn(p,1); u = u/norm(u); % or u = randn(p,1) / sqrt(p);
  
  Q = genquadricov(X, u);
  Q1 = real(Q); Q2 = imag(Q);
  
  [W,S,~] = svd(Q1);
  [~,inds] = sort(diag(S),'descend');
  W = W(:,inds); % just in case 
  W = W(:,1:k);
  
  q1 = W'*Q1*W;
  q2 = W'*Q2*W;
  
  M = q1/q2;
  [V,~] = eig(M);
  
  C = W*V;
  
  for j = 1:k
    c = C(:,j);
    
    a = real(c); b = imag(c);
    theta = atan( -( 2 * sum(a.*b) ) / sum( a.^2 - b.^2 ) )/2;
    while theta<0, theta = theta + pi; end
    while theta>2*pi, theta = theta - pi; end
    
    temp = real(exp(1i*theta)*c);
    C(:,j) = temp/norm(temp);
  end
  
  ds = zeros(p,k);
  for j = 1:k
    c = C(:,j);
    [v,s,~] = svd( reshape(c,p,p) );
    [~,inds] = sort(diag(s), 'descend');
    v = v(:,inds);
    ds(:,j) = v(:,1);
  end
  
end

