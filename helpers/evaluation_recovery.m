function [nrec, recov, perfout] = evaluation_recovery(ds_est, ds, th)
% Copyright: Anastasia Podosinnikova 2019

  if nargin ~= 3
    error('wrong input')
  end
  
  ds = normalize_2(ds);
  ds_est = normalize_2(ds_est);
  [perf, perfout] = perf_cos(ds_est,ds);
  k = size(ds, 2);
  recov = zeros(k,1); % number of atoms recovered
  for i = 1:k
    recov(i) = sum( perf < th ) >= i;
  end
  nrec = sum( perf < th );
end

function [test, perfout] = perf_cos(ds_est,ds)
  k = size(ds,2);
  F = zeros(k,k);
  for i = 1:k
    for j = 1:k
      cosij = abs( ds_est(:,i)'*ds(:,j) ) / ( norm(ds_est(:,i)) * norm(ds(:,j)) );
      loc = acos( cosij );
      F(i,j) = real(loc);
    end
  end
  [matching, cost] = HungarianBipartiteMatching(F);
  [perm, ~] = find(sparse(matching));
  [~, perfout] = update_ds_est(ds_est, ds, perm);
  test = zeros(1,k);
  for i = 1:k
    test(i) = acos( abs( ds_est(:,perm(i))'*ds(:,i) ) / ( norm(ds_est(:,perm(i))) * norm(ds(:,i)) ) );
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

function ds = normalize_2(ds)
  k = size(ds,2);
  for i = 1:k
    ds(:,i) = ds(:,i) / norm(ds(:,i));
  end
end
