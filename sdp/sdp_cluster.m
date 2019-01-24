function [ds_est, Ds_est] = sdp_cluster(Hs, k, ctype)
% Copyright: Anastasia Podosinnikova 2019

  if ~( nargin==2 || nargin==3), error('Wrong number of inputs'); end
  if nargin==2, ctype = 'h'; end
  if nargin==3 && ~( strcmp( ctype,'h' ) || strcmp( ctype,'km' ) ), ctype = 'h'; end
  
  p = sqrt( size(Hs,1) );
  nclust = 3*k;
  
  [~, Fs] = extract_basis(Hs, k);
  Dss = zeros(p^2, nclust);
  for irep = 1:nclust
    u = randn(p,1); u = u/norm(u);
    G = u*u'; 
    D = majorize_minimize(G, Fs);
    Dss(:,irep) = D(:);
  end
  
  Ds_est = cluster_Dss(Dss, k, ctype);
  ds_est = approx_ds_from_Ds(Ds_est);
  
end