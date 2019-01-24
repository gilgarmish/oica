function Ds_est = adaptive_deflation(Fs, k, Ds_est)
% Copyright: Anastasia Podosinnikova 2019

  p = sqrt( size(Fs,1) );
  if ~ (nargin==2 || nargin==3), error('wrong number of inputs'); end
  if nargin==2, Ds_est = []; end
  if size(Ds_est,2) < k
    kloc = k - size(Ds_est,2);
    % updating Fs basis
    Fsloc = [Fs Ds_est];
    
    for i = 1:kloc
      [uuu,~,~] = svd(Fsloc);
      Fsloc = uuu(:,1:end-(kloc-i+1));
      Esloc = uuu(:,end-(kloc-i):end);
      G = reshape(Esloc(:,1),p,p);
      D = majorize_minimize(G, Fsloc);
      Ds_est = [Ds_est D(:)]; %#ok
      Fsloc = [Fsloc D(:)]; %#ok
    end
    
  end
end
