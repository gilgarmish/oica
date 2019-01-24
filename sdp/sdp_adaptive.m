function [ds_est, Ds_est] = sdp_adaptive(Hs, k)
% Copyright: Anastasia Podosinnikova 2019

  if nargin~=2, error('wrong input'); end
  
  [~, Fs] = extract_basis(Hs, k);
  Ds_est = adaptive_deflation(Fs, k);
  ds_est = approx_ds_from_Ds(Ds_est);
  
end