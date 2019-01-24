function experiment_fixedk( p, k, nrep, algs, ns, isshort )
% p is the observed dimension
% k is the latent dimension
% nrep the number of times the data is resampled (experiment repeated)
% algs a cell array with the algorithms names as strings
%      possible entries include 'foobi', 'fpca', 'overica', 'rica', 'rand'
% ns is the array of sample sizes

% Copyright: Anastasia Podosinnikova 2019

  if nargin~=6, error('Wrong input'); end
  
  expopts.('p') = p;
  expopts.('k') = k;
  expopts.('nrep') = nrep;
  expopts.('ns') = ns;
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
  
  disp('Computing quadricovariances...')
  compute_quadricov_bases( expopts, algs );
  
  for ialg = 1:length(algs)
    alg = algs{ialg};
    disp(['Running ', alg])
    run_alg( expopts, alg )
  end
  
  
  if isshort == 1
    islog = 0;
    make_nerr_plots( expopts, islog )
  else
    islog = 0;
    make_nerr_plots( expopts, islog )
    islog = 1;
    make_nerr_plots( expopts, islog )
    make_recovery_plots( expopts )
  end
  
end

function run_alg( expopts, alg )

  ns = expopts.ns;
  expdir = expopts.expdir;
  nrep = expopts.nrep;
  k = expopts.k;
  
  for n = ns
    disp(['     n = ', num2str(n)])
    data = get_datan( expdir, n );
    
    algpath = get_algpath( expdir, alg, n );
    if exist( strcat( algpath, '.mat' ), 'file' ) == 2
      ll = load( algpath );
      expres = ll.expres;
    else
      expres = cell(nrep, 1);
      for irep = 1:nrep

        if strcmp( alg, 'oica-quad-semiada' )
          Hs = data{irep}.('Hs_quad');
          time_cum = data{irep}.('time_cum');
          tt = tic;
          ds_est = sdp_semiada(Hs, k);
          time = toc(tt) + time_cum;
          expres{irep}.('ds_est') = ds_est;
          expres{irep}.('time') = time;
        end

        if strcmp( alg, 'oica' )
          X = data{irep}.('X');
          opts.('sub') = 'gencov';
          opts.('s') = 10;
          opts.('sdp') = 'semiada';
          tt = tic;
          [ds_est, ~, times] = overica( X, k, opts );
          time = toc(tt);
          expres{irep}.('ds_est') = ds_est;
          expres{irep}.('time') = time;
          expres{irep}.('times') = times;
          expres{irep}.('opts') = opts;
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
          ds_est = foobi_core(Hs, k); 
          time = toc(tt) + time_cum;
          expres{irep}.('ds_est') = ds_est;
          expres{irep}.('time') = time;
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
    
    resn = compute_resn(data, expres, expopts, k); %#ok
    save( get_resnpath( expdir, alg, n ), 'resn' )

  end

end


function resn = compute_resn(data, expres, expopts, k)
  nrep = expopts.nrep;
  theta = acos(.99);
  
  resn.('times_loc') = zeros(nrep, 1);
  resn.('aerrs_loc')  = zeros(nrep, 1);
  resn.('ferrs_loc')  = zeros(nrep, 1);
  resn.('recveck')   = zeros(k,1);
  
  for irep = 1:nrep
    ds = data{irep}.ds;
    ds_est = expres{irep}.ds_est;
    
    time = expres{irep}.time;
    resn.times_loc(irep) = time;
    
    resn.aerrs_loc(irep) = a_error(ds_est, ds);
    resn.ferrs_loc(irep) = f_error(ds_est, ds);
    
    [~, recov] = evaluation_recovery(ds_est, ds, theta);
    resn.recveck = resn.recveck + recov;
  end
end



function expname = get_expname( expopts )
  p = expopts.p;
  k = expopts.k;
  expname = strcat( 'fixk', num2str(k), 'p', num2str(p) );
end
function data = get_datan( expdir, n )
  datapath = get_datapath( expdir, n );
  ll = load( datapath );
  data = ll.data;
end
function datapath = get_datapath( expdir, n )
  datapath = strcat( expdir, '/datan', num2str(n) );
end
function algpath = get_algpath( expdir, alg, n )
  algpath = strcat( expdir, '/', alg, 'n', num2str(n) );
end
function resn = get_resn( expdir, alg, n )
  resnpath = get_resnpath( expdir, alg, n );
  ll = load( resnpath );
  resn = ll.resn;
end
function resnpath = get_resnpath( expdir, alg, n )
  resnpath = strcat( expdir, '/resn', num2str(n), alg );
end


function sample_data( expopts )

  expdir = expopts.expdir;
  ns = expopts.ns;
  nrep = expopts.nrep;
  p = expopts.p;
  k = expopts.k;
  
  for n = ns
    datanpath = get_datapath( expdir, n );
    
    if ~ ( exist( strcat( datanpath,'.mat' ), 'file' ) == 2 )
      data = cell(nrep, 1);
      for irep = 1:nrep
        ds = sample_mixing_matrix(p,k);
        X = sample_from_ica_with_uniform_sources(ds, n);
        data{irep}.('ds') = ds;
        data{irep}.('X') = X;
      end
      save( datanpath, 'data' )
    end
    
  end
  
end

function compute_quadricov_bases( expopts, algs )
  
  expdir = expopts.expdir;
  nrep = expopts.nrep;
  ns = expopts.ns;
  k = expopts.k;
  
  isquad = 0;
  for ialg = 1:length(algs)
    alg = algs(ialg);
    if strcmp( alg, 'foobi' ), isquad = 1; end
    if strcmp( alg, 'oica-quad-clust' ), isquad = 1; end
    if strcmp( alg, 'oica-quad-ada' ), isquad = 1; end
    if strcmp( alg, 'oica-quad-semiada' ), isquad = 1; end
    if isquad==1, break; end
  end
  
  for n = ns
    disp(['    n = ', num2str(n)])
    data = get_datan( expdir, n );
    
    for irep = 1:nrep
      if isquad == 1 && ~ isfield( data{irep}, 'Hs_quad' )
        tt = tic;
        X = data{irep}.X;
        C = quadricov(X);
        [CU,~,~] = svd(C);
        Hs = CU(:, 1:k);
        time_cum = toc(tt);
        data{irep}.('Hs_quad') = Hs;
        data{irep}.('time_cum') = time_cum;
      end
    end
    
    save( get_datapath( expdir, n ), 'data' )
  end
  
end



function make_nerr_plots( expopts, islog )

  ns = expopts.ns;
  algs = expopts.algs;
  nalg = length(algs);
  nrep = expopts.nrep;
  
  aerr_min = zeros(length(ns), nalg);
  aerr_max = zeros(length(ns), nalg);
  aerr_med = zeros(length(ns), nalg);

  ferr_min = zeros(length(ns), nalg);
  ferr_max = zeros(length(ns), nalg);
  ferr_med = zeros(length(ns), nalg);

  for j = 1:length(ns)
    n = ns(j);

    aerrs_loc = zeros(nrep, nalg);
    ferrs_loc = zeros(nrep, nalg);

    for i = 1:length(algs)
      resn = get_resn( expopts.expdir, algs{i}, n );
      aerrs_loc(:,i) = resn.aerrs_loc;
      ferrs_loc(:,i) = resn.ferrs_loc;
      if islog == 1
        aerrs_loc(:,i) = log( aerrs_loc(:,i) );
        ferrs_loc(:,i) = log( ferrs_loc(:,i) );
      end
    end

    aerr_min( j, : ) = min(aerrs_loc)';
    aerr_max( j, : ) = max(aerrs_loc)';
    aerr_med( j, : ) = median(aerrs_loc)';

    ferr_min( j, : ) = min(ferrs_loc)';
    ferr_max( j, : ) = max(ferrs_loc)';
    ferr_med( j, : ) = median(ferrs_loc)';

  end
  
  
  plotopts.('xname') = 'Sample Size ($n$; in thousands)';
  plotopts.('xlims') = [ns(1) ns(end)];
  plotopts.('xtickpos') = ns(1:2:end);
  plotopts.('xticks') = ns(1:2:end) / 1000;
  plotopts.('legends') = get_legends( algs );
  plotopts.('xlims') = [ns(1) ns(end)];
  plotopts.('algs') = algs;
  
  
  plotopts.('yname') = 'A-Error';
  plotopts.('ylims') = [0 1];
  if islog == 1
    plotopts.('yname') = 'A-Error ($\log$-linear)';
    plotopts.('ylims') = [ min(min(aerr_min)) max( max(max(aerr_max)), 1 ) ];
  end
  [xs, ys, ys_L, ys_U] = plotdata2cell( ns, aerr_med, aerr_min, aerr_max );
  make_single_plot_cells(xs, ys, ys_L, ys_U, plotopts);
  
  
  plotopts.('yname') = 'F-Error';
  plotopts.('ylims') = [ 0 max( max(max(ferr_max)), 1 ) ];
  if islog == 1
    plotopts.('yname') = 'F-Error ($\log$-linear)';
    plotopts.('ylims') = [ min(min(ferr_min)) max( max(max(ferr_max)), 1 ) ];
  end
  [xs, ys, ys_L, ys_U] = plotdata2cell( ns, ferr_med, ferr_min, ferr_max );
  make_single_plot_cells(xs, ys, ys_L, ys_U, plotopts);
   
  
end

function [xs, ys, ys_L, ys_U] = plotdata2cell( ns, data_med, data_min, data_max )
  nalg = size(data_med,2);
  xs = cell(nalg,1);
  ys = cell(1,nalg);
  ys_L = cell(1,nalg);
  ys_U = cell(1,nalg);
  for i = 1:nalg
    xs{i} = ns;
    ys{i} = data_med(:,i);
    ys_L{i} = data_min(:,i);
    ys_U{i}  = data_max(:,i);
  end
end




function make_recovery_plots( expopts )

  k = expopts.k;
  ns = expopts.ns;
  nrep = expopts.nrep;
  theta = acos(.99);
  algs = expopts.algs;
  expdir = expopts.expdir;
  
  rr = length(ns);
  lenx = length(ns);
  
  for i = 1:length(algs)
    alg = algs{i};
    
    recmat  = zeros( k, length(ns) );
    for j = 1:length(ns)
      n = ns(j);
      resn = get_resn( expdir, alg, n );
      recmat(1:k, j) = resn.recveck;
    end
    
    nrecmat = zeros(k, lenx);
    for j = 1:lenx
      nrecmat(:,j) = recmat(:,j) / nrep;
    end

    nrecmat = abs(nrecmat - 1);


    ff=figure; hold on
    
      % the size and position of the figure on the screen
      screensize = get( groot, 'Screensize' );
      position = [0 0 screensize(3) screensize(4)];
      set(ff, 'Position', position)

      set(gcf, 'Color', 'None')

      % font type for text and ticks
      fontname = 'Times New Roman';
      set(gca, 'Fontname', fontname)

      fontsize = 80;
      set(gca, 'FontSize', fontsize)
      
      angle = round( theta * 180 / pi );
      algname = get_algname( alg );
      titlestr = sprintf('%s   (Angle=%d)', algname, angle );
      title(titlestr, 'FontSize', fontsize, 'Interpreter', 'latex')

      k = expopts.k;
      ylabel('Perfect Recovery', 'Interpreter', 'latex')

      xlim([1 rr])
      ylim([1 k])

      xlabel('Sample Size ($n$; in thousands)', 'Interpreter', 'latex')
      xl = xlim(); 
      inds = xl(1):2:xl(2);
      set(gca,'XTick', inds )
      set(gca,'XTickLabel', ns(inds) / 1000 )

      
      pcolor(nrecmat)
      colormap('gray')
      shading flat

      if sum(sum(nrecmat==0))
        pcolor( ones( size(nrecmat) ) )
        shading flat
        colormap('gray')
        pcolor( nrecmat )
        shading flat
        colormap('gray')
      end
      
      %box on
      ax = gca;

      x1 = ax.XLim(1);
      x2 = ax.XLim(2);
      y1 = ax.YLim(1);
      y2 = ax.YLim(2);

      % instead of box on
      line([x1 x2], [y1 y1], 'Color', 'k')
      line([x1 x2], [y2 y2], 'Color', 'k')
      line([x1 x1], [y1 y2], 'Color', 'k')
      line([x2 x2], [y1 y2], 'Color', 'k')

    hold off
      
  end
end

