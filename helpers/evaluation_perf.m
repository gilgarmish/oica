function [perf, perfout] = evaluation_perf(ds_est, ds, perf_type)
% Copyright: Anastasia Podosinnikova 2019

  if ~( (nargin==2) || (nargin==3) )
    error('perf: wrong input')
  end

  if sum( sum( size(ds_est) == size(ds) ) ) ~= 2
    error('perf: wrong input')
  end

  if nargin==2
    perf_type = 2;
  end
  
  ds = normalize_2(ds);
  ds_est = normalize_2(ds_est);

  switch perf_type
    case 1 % l1
      [perf, perfout] = perf_l1(ds_est, ds);
    case 2 % l2
      [perf, perfout] = perf_l2(ds_est, ds);
    case 3 % fro
      [perf, perfout] = perf_fro(ds_est, ds);
    case 4 % cos
      [perf, perfout] = perf_cos(ds_est, ds);
    otherwise
      error('perf: specify type!')
  end

end


function [perf, perfout] = perf_l2(ds_est, ds)
  k = size(ds,2);
  
  F = zeros(k);
  for i=1:k
    for j=1:k 
      sij = sign( ds_est(:,i)'*ds(:,j) );
      F(i,j) = norm( sij*ds_est(:,i) - ds(:,j), 2);  
    end
  end
  
  [matching, cost] = HungarianBipartiteMatching(F);
  
  [perm, ~] = find(sparse(matching));
  [~, perfout] = update_ds_est(ds_est, ds, perm);
  
  perf = cost/k;
  perf = perf/sqrt(2);  % to scale in [0,1]
  
end

function [perf, perfout] = perf_fro(ds_est, ds)

  %perm = hungarian(-abs(ds_est'*ds) );
  
  matching = HungarianBipartiteMatching(-abs(ds_est'*ds));
  [perm, ~] = find(sparse(matching));

  [ds_est, perfout] = update_ds_est(ds_est, ds, perm);

  perf = norm(ds_est-ds,'fro').^2 / norm(ds,'fro').^2;

end


function [perf, perfout] = perf_l1(ds_est, ds)

  k = size(ds,2);
  
  F = zeros(k);
  for i = 1:k 
    for j = 1:k 
      sij = sign( ds_est(:,i)'*ds(:,j) );
      F(i,j) = norm( sij*ds_est(:,i)-ds(:,j), 1); 
    end
  end
  
  [matching, cost] = HungarianBipartiteMatching(F);
  
  [perm, ~] = find(sparse(matching));
  [~, perfout] = update_ds_est(ds_est, ds, perm);
  
  perf = cost/k;
  perf = perf/2; % to scale in [0,1]
end

function [perf, perfout] = perf_cos(ds_est, ds)

  k = size(ds,2);

  F = zeros(k,k);
  
  for i = 1:k
    for j = 1:k
      cosij = abs( ds_est(:,i)'*ds(:,j) ) / ( norm(ds_est(:,i)) * norm(ds(:,j)) );
      loc = acos( cosij );
%       if abs(imag(loc)) > 0
%         loc
%         %disp('keyboard!')
%         %keyboard
%       end
      F(i,j) = real(loc);
    end
  end
  
  if sum(sum(imag(F))) > 0
    1;
  end

  [matching, cost] = HungarianBipartiteMatching(F);
  
  [perm, ~] = find(sparse(matching));
  [~, perfout] = update_ds_est(ds_est, ds, perm);
  
  perf = 2*cost/pi;
  perf = perf/k;

end


% function ds = normalize_1(ds)
%   k = size(ds,2);
%   for i = 1:k
%     ds(:,i) = ds(:,i) / norm(ds(:,i),1);
%   end
% end

function ds = normalize_2(ds)
  k = size(ds,2);
  for i = 1:k
    ds(:,i) = ds(:,i) / norm(ds(:,i));
  end
end

function [ds_est, perfout] = update_ds_est(ds_est, ds, perm)
  k = size(ds_est, 2);
  ds_est = ds_est(:, perm);
  signs = zeros(k,1);
  for i = 1:k
    signi = sign( ds_est(:,i)' * ds(:,i) );
    signs(i) = signi;
    ds_est(:,i) = ds_est(:,i) * signi; 
  end
  perfout = struct('ds_est', ds_est, 'perm', perm, 'signs', signs);
end


