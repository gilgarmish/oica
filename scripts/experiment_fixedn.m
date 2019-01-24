function expopts = experiment_fixedn( p, ks, nrep, algs, n )
% p is the observed dimension
% ks is the array of latent dimensions k
% nrep the number of times the data is resampled (experiment repeated)
% algs a cell array with the algorithms names as strings
%      possible entries include 'foobi', 'fpca', 'overica', 'rica', 'rand'
% n is the sample size

% Copyright: Anastasia Podosinnikova 2019

  if nargin~=5, error('Wrong input'); end
  
  expopts.('p') = p;
  expopts.('ks') = ks;
  expopts.('nrep') = nrep;
  expopts.('n') = n;
  expopts.('algs') = algs;
  expname = get_expname( expopts );
  expopts.('expname') = expname;
  expdir = strcat( pwd, '/expres/', expname );
  expopts.('expdir') = expdir;
  
  disp(['Data will be saved to: ', expdir])
  if exist( expdir, 'dir' ) ~= 7
    mkdir(expdir)
  end
  save( strcat( expdir, '/expopts' ), 'expopts' )
  
  disp('Sampling data...')
  sample_data( expopts );
  
  disp('Computing Quadricovariances...')
  compute_quadricov_bases( expopts, algs );
  
  for ialg = 1:length(algs)
    alg = algs{ialg};
    disp(['Running ', alg])
    run_alg( expopts, alg )
  end
  
  make_time_plot( expopts )
  make_nerr_plots( expopts )

  
end

function run_alg( expopts, alg )

  ks = expopts.ks;
  expdir = expopts.expdir;
  nrep = expopts.nrep;
  
  for k = ks
    disp(['     k = ', num2str(k)])
    data = get_datak( expdir, k );
    
    algpath = get_algpath( expdir, alg, k );

    if exist( strcat( algpath, '.mat' ), 'file' ) == 2
      ll = load( algpath );
      expres = ll.expres;
    else
      expres = cell(nrep, 1);
      for irep = 1:nrep
        
        if strcmp( alg, 'oica' ) % oica-gencov-semiada
          X = data{irep}.('X');
          opts.('sub') = 'gencov';
          opts.('s') = 10;
          opts.('t') = 0.1;
          opts.('sdp') = 'semiada';
          tt=tic; [ds_est, ~, times] = overica( X, k, opts ); time=toc(tt);
          expres{irep}.('ds_est') = ds_est;
          expres{irep}.('time')   = time;
          expres{irep}.('times')  = times;
          expres{irep}.('opts')   = opts;
        end

        if strcmp( alg, 'fpca' )
          X = data{irep}.('X');
          tt = tic;
          ds_est = fourier_pca2(X, k);
          time = toc(tt);
          expres{irep}.('ds_est') = ds_est;
          expres{irep}.('time') = time;
        end

        if strcmp( alg, 'foobi' )
          Hs = data{irep}.('Hs_quad');
          time_cum = data{irep}.('time_cum');
          tt = tic;
          [ds_est, Ds_est] = foobi_core(Hs, k); 
          time_core = toc(tt);
          time = time_cum + time_core;
          expres{irep}.('ds_est') = ds_est;
          expres{irep}.('Ds_est') = Ds_est;
          expres{irep}.('time') = time;
          expres{irep}.('time_cum') = time_cum;
          expres{irep}.('time_core') = time_core;
        end

        if strcmp( alg, 'rand' )
          p = expopts.p;
          ds_est = sample_mixing_matrix(p,k);
          expres{irep}.('ds_est') = ds_est;
          expres{irep}.('time') = 0;
        end

      end
      save( algpath, 'expres' );

    end

    resk = compute_resk(data, expres, expopts, k); %#ok
    save( get_reskpath( expdir, alg, k ), 'resk' )

  end

end


function resk = compute_resk(data, expres, expopts, k)
  nrep = expopts.nrep;
  theta = acos(.99);
  resk.('times_loc') = zeros(nrep, 1);
  resk.('aerrs_loc') = zeros(nrep, 1);
  resk.('ferrs_loc') = zeros(nrep, 1);
  resk.('recveck')   = zeros(k,1);
  
  for irep = 1:nrep
    ds = data{irep}.ds;
    ds_est = expres{irep}.ds_est;
    
    time = expres{irep}.time;
    resk.times_loc(irep) = time;
    
    resk.aerrs_loc(irep) = a_error(ds_est, ds);
    resk.ferrs_loc(irep) = f_error(ds_est, ds);
    
    [~, recov] = evaluation_recovery(ds_est, ds, theta);
    resk.recveck = resk.recveck + recov;
  end
end


function sample_data( expopts )

  expdir = expopts.expdir;
  ks = expopts.ks;
  nrep = expopts.nrep;
  p = expopts.p;
  n = expopts.n;
  
  for k = ks
    datapath = get_datapath( expdir, k );
    if ~( exist( strcat( datapath,'.mat' ), 'file' ) == 2 )
      data = cell(nrep, 1);
      for irep = 1:nrep
        ds = sample_mixing_matrix(p,k);
        data{irep}.('ds') = ds;
        X = sample_from_ica_with_uniform_sources(ds, n);
        X = X - repmat( mean(X,2), 1, size(X,2) );
        data{irep}.('X') = X;
      end
      save( datapath, 'data' )
    end
  end
  
end


function expname = get_expname( expopts )
  p = expopts.p;
  n = expopts.n;
  expname = strcat( 'fixn', num2str(n), 'p', num2str(p) );
end
function data = get_datak( expdir, k )
  datak = get_datapath( expdir, k );
  ll = load( datak );
  data = ll.data;
end
function datapath = get_datapath( expdir, k )
  datapath = strcat( expdir, '/datak', num2str(k) );
end
function algpath = get_algpath( expdir, alg, k )
  algpath = strcat( expdir, '/', alg, 'k', num2str(k) );
end
function resk = get_resk( expdir, alg, k )
  reskpath = get_reskpath( expdir, alg, k );
  ll = load( reskpath );
  resk = ll.resk;
end
function reskpath = get_reskpath( expdir, alg, k )
  reskpath = strcat( expdir, '/resk', num2str(k), alg );
end


function compute_quadricov_bases( expopts, algs )
  
  expdir = expopts.expdir;
  nrep = expopts.nrep;
  ks = expopts.ks;
  
  isquad = 0;
  for ialg = 1:length(algs)
    alg = algs(ialg);
    if strcmp( alg, 'foobi' ), isquad = 1; end
    if strcmp( alg, 'oica-quad-clust' ), isquad = 1; end
    if strcmp( alg, 'oica-quad-ada' ), isquad = 1; end
    if strcmp( alg, 'oica-quad-semiada' ), isquad = 1; end
    if isquad==1, break; end
  end
  
  for k = ks
    disp(['    k = ',num2str(k)])
    data = get_datak( expdir, k );
    issave = 0;
    
    for irep = 1:nrep
      
      if ~isfield(data{irep}, 'Hs_quad') && isquad == 1
        tt = tic;
        X = data{irep}.X;
        C = quadricov(X);
        [CU,~,~] = svd(C);
        Hs = CU(:, 1:k);
        time_cum = toc(tt);
        data{irep}.('Hs_quad') = Hs;
        data{irep}.('time_cum') = time_cum;
        issave = 1;
      end
      
    end
    
    if issave == 1, save( get_datapath( expdir, k ), 'data' ); end
    
  end
  
end

function make_time_plot( expopts )

  algs = algs_exclude_rand( expopts.algs );
  nalg = length(algs);
  ks = expopts.ks;
  nrep = expopts.nrep;

  times_min = zeros(length(ks), nalg);
  times_max = zeros(length(ks), nalg);
  times_med = zeros(length(ks), nalg);
  
  for j = 1:length(ks)
    k = ks(j);
    times_loc = zeros( nrep, nalg );
    for i = 1:nalg
      resk = get_resk( expopts.expdir, algs{i}, k );
      times_loc(:, i) = log( resk.times_loc );
    end
    times_min( j, : ) = min(times_loc)';
    times_max( j, : ) = max(times_loc)';
    times_med( j, : ) = median(times_loc)';
  end
  
  plotopts.('yname') = strcat('Runtime ($\log$-Scale)');
  plotopts.('xname') = 'Latent Dimension ($k$)';
  plotopts.('legends') = get_legends( algs );
  plotopts.('algs' ) = algs;
  
  plotopts.('greenlines') = expopts.p;
  plotopts.('ylims') = [min(min(times_min)) max(max(times_max))];
  plotopts.('xlims') = [ks(1) ks(end)];
  plotopts.('xtickpos') = ks;
  plotopts.('xticks') = ks;
  
  [xs, ys, ys_L, ys_U] = plotdata2cell( ks, times_med, times_min, times_max );
  make_single_plot_cells(xs, ys, ys_L, ys_U, plotopts);
  

end

function make_nerr_plots( expopts )
  
  ks = expopts.ks;
  algs = expopts.algs;
  nalg = length(algs);
  nrep = expopts.nrep;
  
  aerr_min = zeros(length(ks), nalg);
  aerr_max = zeros(length(ks), nalg);
  aerr_med = zeros(length(ks), nalg);

  ferr_min = zeros(length(ks), nalg);
  ferr_max = zeros(length(ks), nalg);
  ferr_med = zeros(length(ks), nalg);

  
  for j = 1:length(ks)
    k = ks(j);

    aerrs_loc = zeros(nrep, nalg);
    ferrs_loc = zeros(nrep, nalg);

    for i = 1:length(algs)
      resn = get_resk( expopts.expdir, algs{i}, k );
      aerrs_loc(:,i) = resn.aerrs_loc;
      ferrs_loc(:,i) = resn.ferrs_loc;
    end

    aerr_min( j, : ) = min(aerrs_loc)';
    aerr_max( j, : ) = max(aerrs_loc)';
    aerr_med( j, : ) = median(aerrs_loc)';

    ferr_min( j, : ) = min(ferrs_loc)';
    ferr_max( j, : ) = max(ferrs_loc)';
    ferr_med( j, : ) = median(ferrs_loc)';

  end
  
  plotopts.('xname') = 'Latent Dimension ($k$)';
  plotopts.('legends') = get_legends( expopts.algs );
  plotopts.('algs' ) = algs;
  plotopts.('greenlines') = expopts.p;
  plotopts.('xlims') = [ks(1) ks(end)];
  plotopts.('xtickpos') = ks;
  plotopts.('xticks') = ks;
  
  plotopts.('yname') = 'A-Error';
  plotopts.('ylims') = [0 1];
  [xs, ys, ys_L, ys_U] = plotdata2cell( ks, aerr_med, aerr_min, aerr_max );
  make_single_plot_cells(xs, ys, ys_L, ys_U, plotopts);
  
  plotopts.('yname') = 'F-Error';
  plotopts.('ylims') = [ 0 max( max(max(ferr_max)), 1 ) ];
  [xs, ys, ys_L, ys_U] = plotdata2cell( ks, ferr_med, ferr_min, ferr_max );
  make_single_plot_cells(xs, ys, ys_L, ys_U, plotopts);
  
end

function [xs, ys, ys_L, ys_U] = plotdata2cell( ks, data_med, data_min, data_max )
  nalg = size(data_med,2);
  xs = cell(nalg,1);
  ys = cell(1,nalg);
  ys_L = cell(1,nalg);
  ys_U = cell(1,nalg);
  for i = 1:nalg
    xs{i} = ks;
    ys{i} = data_med(:,i);
    ys_L{i} = data_min(:,i);
    ys_U{i}  = data_max(:,i);
  end
end
