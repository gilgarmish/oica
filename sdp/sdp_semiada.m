function [ds_est, Ds_est] = sdp_semiada(Hs, k, ctype)
% Copyright: Anastasia Podosinnikova 2019

  if ~( nargin==2 || nargin==3), error('Wrong number of inputs'); end
  if nargin==2, ctype = 'h'; end
  if nargin==3 && ~( strcmp( ctype,'h' ) || strcmp( ctype,'km' ) ), ctype = 'h'; end
  
  [~, Ds_clust] = sdp_cluster(Hs, k, ctype);
  
  % keep only the atoms which are well separated
  G = abs(Ds_clust'*Ds_clust) - eye(k);
  [~, mind] = min(max(G)); % choose the most separated atom
  Ds_est(:,1) = Ds_clust(:,mind);
  ind = 2;
  for i = setdiff( 1:k, mind )
    if max( G(:,i) ) < .8
      Ds_est(:,ind) = Ds_clust(:,i);
      ind = ind + 1;
    end
  end
  
  [~, Fs] = extract_basis(Hs, k);
  Ds_est = adaptive_deflation(Fs, k, Ds_est);
  ds_est = approx_ds_from_Ds(Ds_est);
  
end