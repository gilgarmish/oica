function reproduce_fig_4_patches_atoms( data )
% This script reproduces Figure 4 displaying the overcomplete atoms
%     estimated from CIFAR-10 (training batch 1)
% Dependency: CIFAR-10 (https://www.cs.toronto.edu/~kriz/cifar.html)

% data is the training set from batch 1 of CIFAR-10

% Copyright: Anastasia Podosinnikova 2019

  p = 7;
  X = get_patches(data, p); % converts to grayscale & returns 7x7 patches
  
  XC = X - repmat( mean(X,2), 1, size(X,2) );
  opts.('sdp') = 'semiada';
  opts.('t') = 0.1;
  k = 150;
  
  rand('stat',0)
  randn('stat',0)
  ds_est = overica( XC, k, opts );
  
  plot_atoms( ds_est, p, k, 10, 15 )

end

function X = get_patches(data, p)
  sr = 32; sc = 32; % image size is sr x sc
  [N, ~] = size(data);
  d = floor(p/2) + mod(p,2) - 1;
  step = (sr - 2*d - 1)*(sc - 2*d - 1);
  X = zeros( p*p, N * step );
  ind = 1;

  for n = 1:N
    disp(['n = ', num2str(n)])
    rr = data(n,:);
    gg = rgb2gray( reshape( rr, 32, 32, 3 ) );
    pp = nlfilter(gg, [p p], @(block) {block});

    for i = 1:32
      for j = 1:32
        patch = pp{i,j};
        X(:,ind) = patch(:);
        ind = ind + 1;
      end
    end
  end

end

function plot_atoms( ds, p, k, a, b )
  
   d1 = ds(:,1);
   for i = 2:k
     ds(:,i) = sign(d1'*ds(:,i))*ds(:,i);
   end
  ds = ds./(ones(size(ds,1),1)*max(abs(ds)));

  ff=figure; hold on
    screensize = get( groot, 'Screensize' );
    rr = min(screensize(3), screensize(4));
    set(ff,'Position', [0 0 rr*b/a rr])

    for i = 1:k
      subplot(a,b,i)
      rr = reshape( ds(:,i), [p p] );
      pcolor(rr)
      shading flat
      colormap('gray')
      axis off
    end

  hold off
end
