function [ds_est, Ds_est, times, Hs] = overica( X, k, opts )
% Overcomplete Independent Component Analysis via SDP (OverICA) introduced in:
%
%  Anastasia Podosinnikova, Amelia Perry, Alexander Wein, Francis Bach,
%  Alexandre d'Aspremont, David Sontag. Overcomplete Independent Component
%  Analysis via SDP. In Proceedings of the 22nd International Conference on
%  Artificial Intelligence and Statistics (AISTATS), 2019.
%
% Please cite this paper if you use this code for your research.
%
% Input:
%   X is the p-times-n data matrix with observations in columns
%   k is the desired latent dimension (k < p^2/4 for our guarantees to hold)
%   opts (optional):
%      * opts.sub can take either 'quad' (uses quadricovariance, i.e. the
%        fourth order cumulant) or 'gencov' (uses generalized covariance; default)
%      * opts.s is defines the number (s*k) of generalized covariances to be
%        used (default is s = 5, i.e. 5*k gencovs will be constructed)
%      * opts.t is the parameter for gencovs (we found that t = 0.1 works
%        well in practice; default is t = 0.05/sqrt(max(max(abs(cov(X')))));
%      * opts.sdp can take either 'clust' or 'ada' or 'semiada' (default) 
%        and defines the type of deflation procedure to be used
%        (clustering-based vs adaptive vs semi-adaptive)

% Copyright: Anastasia Podosinnikova 2019

  if ~( nargin==2 || nargin==3), error('Wrong number of inputs'); end
  if nargin==2, opts = set_default_opts; end
  if nargin==3, opts = check_opts(opts); end
  
  globtt = tic;
  if strcmp( opts.sub, 'quad' )
    disp('Computing quadricovariance')
    tt=tic; C = quadricov(X); toc(tt)
  end
  
  if strcmp( opts.sub, 'gencov' )
    disp('Computing generalized covariances')
    tt = tic;
    s = 5*k;
    t = 0.05/sqrt( max( max( abs( cov(X') ) ) ) );
    if isfield( opts, 's' ), s = opts.s*k; end
    if isfield( opts, 't' ), t = opts.t; end
    C = estimate_gencovs(X, s, t);
    toc(tt)
  end
  toc(globtt)
  cum_time = toc(globtt);
  
  
  disp('Computing SVD')
  [CU,~,~] = svds(C,k);
  Hs = CU(:, 1:k);
  toc(globtt)
  svd_time = toc(globtt) - cum_time;
  
  
  
  if strcmp( opts.sdp, 'clust' )
    disp('Deflation via clustering')
    if isfield( opts, 'ctype' )
      [ds_est, Ds_est] = sdp_cluster(Hs, k, opts.ctype);
    else
      [ds_est, Ds_est] = sdp_cluster(Hs, k);
    end
  end
  if strcmp( opts.sdp, 'ada' )
    disp('Adaptive deflation')
    [ds_est, Ds_est] = sdp_adaptive(Hs, k);
  end
  if strcmp( opts.sdp, 'semiada' )
    disp('Semiadaptive deflation')
    if isfield( opts, 'ctype' )
      [ds_est, Ds_est] = sdp_semiada(Hs, k, opts.ctype);
    else
      [ds_est, Ds_est] = sdp_semiada(Hs, k);
    end
  end
  toc(globtt)
  sdp_time = toc(globtt) - cum_time - svd_time;
  
  total_time = toc(globtt);
  
  times.('total_time') = total_time;
  times.('cum_time')   = cum_time;
  times.('svd_time')   = svd_time;
  times.('sdp_time')   = sdp_time;
  
end

function C = estimate_gencovs(X, s, t0)

  [p, n] = size(X);
  t = t0;
  C = zeros(p^2, s);
  G0 = gencov(X,zeros(p,1)); G0 = G0(:);
  for i = 1:s
    omega = randn(p,1);
    if t0 == -1, t = choose_t(X'*omega, n/20); end
    omega = t*omega;
    Gi = gencov(X,omega);
    C(:,i) = Gi(:) - G0;
  end
  
end

function opts = check_opts(opts)
  if ~isstruct(opts), opts = set_default_opts; end
  if ~isfield( opts, 'sub' ), opts.('sub') = 'gencov'; end
  if ~isfield( opts, 'sdp' ), opts.('sdp') = 'semiada'; end
  if isfield( opts, 'ctype' )
    ctype = opts.ctype;
    if ~( strcmp( ctype, 'h' ) || strcmp( ctype, 'km' ) )
      disp('The opts.ctype value is changed to h')
      opts.('ctype') = 'h'; 
    end
  end
  if isfield( opts, 'sub' )
    sub = opts.sub;
    if ~( strcmp( sub, 'quad' ) || strcmp( sub, 'gencov' ) )
      disp('The opts.sub value is changed to gencov')
      opts.('sub') = 'gencov'; 
    end
  end
  if isfield( opts, 'sdp' )
    sdp = opts.sdp;
    if ~( strcmp( sdp, 'ada' ) || strcmp( sdp, 'semiada' ) || strcmp( sdp, 'clust' ) )
      disp('The opts.sdp value is changed to semiada')
      opts.sdp = 'semiada';
    end
  end
end

function opts = set_default_opts
  opts.('sub') = 'gencov';
  opts.('sdp') = 'semiada';
end
