function experiment_asymptotic( p, ks, nrep, algs )
% p is the observed dimension
% ks is the array of latent dimensions k
% nrep the number of times the data is resampled (experiment repeated)
% algs a cell array with the algorithms names as strings
%      possible entries include 'foobi', 'fpca', 'overica', 'rica', 'rand'

% Copyright: Anastasia Podosinnikova 2019

  if nargin~=4, error('Wrong input'); end
  
  expopts.('p') = p;
  expopts.('ks') = ks;
  expopts.('nrep') = nrep;
  expopts.('algs') = algs;
  expname = get_expname( expopts );
  expdir = strcat( pwd, '/expres/', expname );
  expopts.('expdir') = expdir;
  
  disp(['Data will be saved to: ', expdir])
  if exist( expdir, 'dir' ) ~= 7
    mkdir(expdir)
  end
  save( strcat( expdir, '/expopts' ), 'expopts' )
  
  disp('Sampling data...')
  sample_data( expopts );
  
  disp('Constructing subspaces...')
  construct_subspaces( expopts );
  
  for ialg = 1:length(algs)
    alg = algs{ialg};
    disp(['Running: ', alg])
    run_alg( expopts, alg )
  end
  
  disp('Making plots...')
  %make_time_plot( expopts )
  make_err_plots( expopts )
  expopts.('ks') = ks( logical( ks <= p^2/4 ) );
  %make_err_plots( expopts )
  make_recovery_plots( expopts )
  
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
        Hs = data{irep}.('Hs');
        
        if strcmp( alg, 'oica-clust' )
          tt = tic;
          ds_est = sdp_cluster(Hs, k, 'h');
          time = toc(tt);
          expres{irep}.('ds_est') = ds_est;
          expres{irep}.('time') = time;
        end
        
        if strcmp( alg, 'oica-ada' )
          tt = tic;
          ds_est = sdp_adaptive(Hs, k);
          time = toc(tt);
          expres{irep}.('ds_est') = ds_est;
          expres{irep}.('time') = time;
        end
        
        if strcmp( alg, 'oica-semiada' )
          tt = tic;
          ds_est = sdp_semiada(Hs, k, 'h');
          time = toc(tt);
          expres{irep}.('ds_est') = ds_est;
          expres{irep}.('time') = time;
        end
        
        if strcmp( alg, 'foobi' )
          tt = tic;
          ds_est = foobi_core(Hs, k); 
          time = toc(tt);
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



function expname = get_expname( expopts )
  p = expopts.p;
  expname = strcat( 'asymp', num2str(p) );
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



function sample_data( expopts )

  expdir = expopts.expdir;
  ks = expopts.ks;
  
  for k = ks
    datapath = get_datapath( expdir, k );
    if exist( strcat( datapath,'.mat' ), 'file' ) == 2
      warning('Data exists. Skipping.');
    else
      nrep = expopts.nrep;
      p = expopts.p;
      data = cell(nrep, 1);
      for irep = 1:nrep
        ds = sample_mixing_matrix(p,k);
        data{irep} = struct( 'ds', ds );
      end
      save( datapath, 'data' )
    end
  end
  
end



function construct_subspaces( expopts )
  
  expdir = expopts.expdir;
  nrep = expopts.nrep;
  ks = expopts.ks;
  
  for k = ks
    data = get_datak( expdir, k );
    
    for irep = 1:nrep
      tt = tic;
      ds = data{irep}.ds;
      Ds = Ds_from_ds(ds);
      C = Ds*Ds';
      [CU,~,~] = svd(C);
      Hs = CU(:,1:k);
      time_cum = toc(tt);
      data{irep}.('Hs') = Hs;
      data{irep}.('time_cum') = time_cum;
    end
    
    save( get_datapath( expdir, k ), 'data' )
    
  end
  
end



function make_time_plot( expopts )

  algs = algs_exclude_rand( expopts.algs );
  nalg = length(algs);
  ks = expopts.ks;
  nrep = expopts.nrep;
  expdir = expopts.expdir;

  times_min = zeros(length(ks), nalg);
  times_max = zeros(length(ks), nalg);
  times_med = zeros(length(ks), nalg);
  
  for j = 1:length(ks)
    k = ks(j);
    times_loc = zeros( nrep, nalg );
    for i = 1:nalg
      resk = get_resk( expdir, algs{i}, k );
      times_loc(:, i) = log( resk.times_loc );
    end
    times_min( j, : ) = min(times_loc)';
    times_max( j, : ) = max(times_loc)';
    times_med( j, : ) = median(times_loc)';
  end
  
  
  
  plotopts.('yname') = 'Runtime ($\log$-linear)';
  plotopts.('ylims') = [min(min(times_min)) max(max(times_max))];
  plotopts.('xname') = 'Latent Dimension ($k$)';
  plotopts.('legends') = get_legends( algs );
  plotopts.('greenlines') = expopts.p;
  plotopts.('algs') = algs;
  plotopts.('xlims') = [ks(1) ks(end)];
  plotopts.('xtickpos') = ks;
  plotopts.('xticks') = ks;
  
  xs = cell(nalg,1);
  ys = cell(1,nalg);
  ys_L = cell(1,nalg);
  ys_U = cell(1,nalg);
  for i = 1:nalg
    xs{i} = ks;
    ys{i} = times_med(:,i);
    ys_L{i} = times_min(:,i);
    ys_U{i}  = times_max(:,i);
  end
  make_single_plot_cells(xs, ys, ys_L, ys_U, plotopts);

end


function make_err_plots( expopts )
  
  ks = expopts.ks;
  nrep = expopts.nrep;
  algs = expopts.algs;
  nalg = length(algs);

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
      resk = get_resk( expopts.expdir, algs{i}, k );
      aerrs_loc(:,i) = resk.aerrs_loc;
      ferrs_loc(:,i) = resk.ferrs_loc;
    end

    aerr_min( j, : ) = min(aerrs_loc)';
    aerr_max( j, : ) = max(aerrs_loc)';
    aerr_med( j, : ) = median(aerrs_loc)';

    ferr_min( j, : ) = min(ferrs_loc)';
    ferr_max( j, : ) = max(ferrs_loc)';
    ferr_med( j, : ) = median(ferrs_loc)';

  end
  
  
  plotopts.('legends') = get_legends( algs );
  plotopts.('algs' ) = algs;
  plotopts.('greenlines') = expopts.p;
  plotopts.('xname') = 'Latent Dimension ($k$)';
  plotopts.('xlims') = [ks(1) ks(end)];
  plotopts.('xtickpos') = ks;
  plotopts.('xticks') = ks;
  
  
  plotopts.('yname') = 'A-Error';
  plotopts.('ylims') = [ 0 max( max(max(aerr_max)), 1 ) ];
  
  xs = cell(nalg,1);
  ys = cell(1,nalg);
  ys_L = cell(1,nalg);
  ys_U = cell(1,nalg);
  for i = 1:nalg
    xs{i} = ks;
    ys{i} = aerr_med(:,i);
    ys_L{i} = aerr_min(:,i);
    ys_U{i}  = aerr_max(:,i);
  end
  make_single_plot_cells(xs, ys, ys_L, ys_U, plotopts);
  
  plotopts.('yname') = 'F-Error';
  plotopts.('ylims') = [ 0 max( max(max(ferr_max)), 1 ) ];
  
  xs = cell(nalg,1);
  ys = cell(1,nalg);
  ys_L = cell(1,nalg);
  ys_U = cell(1,nalg);
  for i = 1:nalg
    xs{i} = ks;
    ys{i} = ferr_med(:,i);
    ys_L{i} = ferr_min(:,i);
    ys_U{i}  = ferr_max(:,i);
  end
  make_single_plot_cells(xs, ys, ys_L, ys_U, plotopts);

  
end


function make_recovery_plots( expopts )

  ks = expopts.ks;
  nrep = expopts.nrep;
  algs = expopts.algs;
  p = expopts.p;
  theta = acos(.99);
  expdir = expopts.expdir;
  
  for ialg = 1:length(algs)
    alg = algs{ialg};

    recmat  = zeros( ks(end), length(ks) );
    for j = 1:length(ks)
      k = ks(j);
      resk = get_resk( expdir, alg, k );
      recmat(1:k, j) = resk.recveck / nrep;
    end

    nrecmat = zeros( size(recmat) );
    for i=1:length(ks)
      nrecmat(end:-1:1, i) = recmat(:,i);
    end

    ff=figure; hold on

      % the size and position of the figure on the screen
      screensize = get( groot, 'Screensize' );
      rr = min( screensize(3), screensize(4) );
      position = [0 0 rr rr];
      set(ff, 'Position', position)

      set(gcf, 'Color', 'None')

      % font type for text and ticks
      fontname = 'Times New Roman';
      set(gca, 'Fontname', fontname)

      % font size for title, axis, ticks
      fontsize = 80;
      set(gca, 'FontSize', fontsize)

      pcolor(recmat)
      shading flat

      colormap(flipud(gray(256)))
      caxis([0,1])
      %colorbar

      angle = round( theta * 180 / pi );

      algname = get_algname( alg );
      titlestr = sprintf('%s   (Angle=%d)', algname, angle );
      title(titlestr, 'FontSize', fontsize, 'Interpreter', 'latex')
      xlabel('Latent Dimension ($k$)', 'Interpreter', 'latex')
      ylabel('Perfect Recovery', 'Interpreter', 'latex')

      ylim([1 ks(end)])

      lenk = length(ks);
      inds = 1:1:lenk;
      xlim([1 lenk])
      set(gca,'XTick', inds);
      set(gca,'XTickLabel', ks(inds));

      pcolor(recmat)
      shading flat

      ax = gca;

      for i = 1:length(ks)-1
        k = ks(i);
        if ( (k <= p) && (ks(i+1) > p) )
          line( [i i], ax.YLim, 'Color', 'g')
        end
      end

      for i = 1:length(ks)-1
        k = ks(i);
        if ( (k <= round(p^2/4)) && (ks(i+1) > round(p^2/4)) )
          line( [i i], ax.YLim, 'Color', 'g')
        end
      end

      for i = 1:length(ks)-1
        k = ks(i);
        if ( (k <= p*(p-1)/2) && (ks(i+1) > p*(p-1)/2) )
          line( [i i], ax.YLim, 'Color', 'g' )
        end
      end

      for i = 1:length(ks)-1
        k = ks(i);
        if ( (k <= p*(p+1)/2) && (ks(i+1) > p*(p+1)/2) )
          line( [i i], ax.YLim, 'Color', 'g' )
        end
      end

      line([1 lenk], [ks(1) ks(end)], 'Color', 'r',  'Linewidth', 2)

      x1 = ax.XLim(1);
      x2 = ax.XLim(2);
      y1 = ax.YLim(1);
      y2 = ax.YLim(2);

      % because box on doesn't work correctly
      line([x1 x2], [y1 y1], 'Color', 'k')
      line([x1 x2], [y2 y2], 'Color', 'k')
      line([x1 x1], [y1 y2], 'Color', 'k')
      line([x2 x2], [y1 y2], 'Color', 'k')

    hold off
  end
end
