 function make_single_plot_cells(xs, ys, ys_L, ys_U, plotopts)
% Copyright: Anastasia Podosinnikova 2019

  if nargin ~= 5, error('wrong input'); end
  
  nalg = length(xs);
  if (length(ys) ~= nalg) || (length(ys_L) ~= nalg) || (length(ys_U) ~= nalg)
    error('wrong input')
  end
  
  if ~isfield( plotopts, 'algs' )
    error('specify algs')
  end
  if length(plotopts.algs) ~= nalg
    error('wrong input (algs vs nalg)')
  end
  
  % plots line specifications
  [LF, LC, LW, MF, MC, MS] = get_alg_specs( plotopts.algs );
  
  ff = figure('Visible','On'); hold on

    % the size and position of the figure on the screen
    screensize = get( groot, 'Screensize' );
    position = [0 0 screensize(3) screensize(4)];
    set(ff, 'Position', position)
    
    % font type for text and ticks
    fontname = 'Times New Roman';
    set(gca, 'Fontname', fontname)
    
    % font size for title, axis, ticks
    fontsize = 80;
    set(gca, 'FontSize', fontsize)
    
    % setting ticks, limits, etc
    if isfield(plotopts, 'title')
      title(plotopts.title, 'FontSize', fontsize, 'Interpreter', 'latex')
    end
    if isfield(plotopts, 'xname')
      xlabel(plotopts.xname, 'Interpreter', 'latex')
    end
    if isfield(plotopts, 'yname')
      ylabel(plotopts.yname, 'Interpreter', 'latex')
    end
    if isfield(plotopts, 'xlims')
      xlim(plotopts.xlims)
    end
    if isfield(plotopts, 'ylims')
      ylim(plotopts.ylims)
    end
    if isfield(plotopts, 'xtickpos')
      set(gca,'XTick', plotopts.xtickpos)
    end
    if isfield(plotopts, 'xticks')
      set(gca,'XTickLabel', plotopts.xticks); 
    end
    if isfield(plotopts, 'ytickpos')
      set(gca,'YTick', plotopts.ytickpos);
    end
    if isfield(plotopts, 'yticks')
      set(gca,'YTickLabel', plotopts.yticks); 
    end

    
    
    % if need to plot green vertical lines
    if isfield(plotopts, 'greenlines')
      p = plotopts.greenlines;
      ks = xs{1};
      SP = p;
      ax = gca;
      
      line([SP SP],ax.YLim,'Color',[0 1 0],'LineWidth',3)

      SP = p^2/4;
      if SP <= ks(end)
        line([SP SP],ax.YLim,'Color',[0 1 0],'LineWidth',3)
      end

      SP = p*(p-1)/2;
      if SP <= ks(end)
        line([SP SP],ax.YLim,'Color',[0 1 0], 'LineWidth',3)
      end

      SP = p*(p+1)/2;
      if SP <= ks(end)
        line([SP SP],ax.YLim,'Color',[0 1 0],'LineWidth',3)
      end
      
    end
    
    
    % set transparent background
    set(gcf, 'Color', 'None')
    
    % make actual plots
    for i = 1:nalg
      xsloc = xs{i};
      ysloc = ys{i};
      ys_Lloc = ys_L{i};
      ys_Uloc = ys_U{i};
      len_x = length(xsloc);
      
      % want to avoide intersecting markers
      mstep = 1/(nalg+1);
      mlen = len_x - 1;
      mxs = zeros(mlen,1);
      mys = zeros(mlen,1);
      for j = 1:mlen
        mxs(j) = xsloc(j) + (1-mstep*i)*(xsloc(j+1)-xsloc(j));
        mys(j) = ysloc(j) + (1-mstep*i)*(ysloc(j+1)-ysloc(j));
      end
      
      plot( [xsloc(1); mxs; xsloc(end)], ...
        [ysloc(1); mys; ysloc(end)], LF{i}, 'Color', LC{i}, 'LineWidth', LW)
      
      inds = 1:mlen;
      errorbar(mxs, mys, ysloc(inds) - ys_Lloc(inds), ys_Uloc(inds) - ysloc(inds), ...
        LF{i}, 'Color', LC{i}, 'LineWidth', LW)
      
      plot( mxs, mys, '.', 'LineWidth', LW,...
        'Marker', MF{i}, 'MarkerSize', MS{i},...
        'MarkerFaceColor', MC{i}, 'MarkerEdgeColor', MC{i});
      
    end
      

    % rotate x-tick-labels by 45 degrees
    if isfield(plotopts, 'XRotate') && plotopts.XRotate == 1
      xtickangle(45)
    end
    
    
    % to display legends correctly
    if isfield( plotopts, 'legends' )
      fakelines = zeros(nalg,1);
      for i = 1:nalg
        xsloc = xs{i};
        fakelines(i) = plot(xsloc, Inf*ones(1,length(xsloc)), ...
          LF{i}, 'Color', LC{i}, 'LineWidth', LW, ...
          'Marker', MF{i}, 'MarkerSize', MS{i}, ...
          'MarkerFaceColor', MC{i}, 'MarkerEdgeColor', MC{i});
      end
      [leg, legobj] = legend(fakelines, plotopts.legends);
      textobj = findobj(legobj, 'type', 'text');
      legendsize = 50;
      set(textobj, 'Interpreter', 'latex', 'fontsize', legendsize);
      set(leg, 'FontSize', legendsize );
      set(leg, 'Location', 'NorthEast')
      set(leg, 'Box', 'on')
    end

    box on
    
  hold off
  
end 


function [LF, LC, LW, MF, MC, MS] = get_alg_specs( algs )

  nalg = length(algs);
  c = get_colors;

  LW = 4; % line width
  LF = cell(nalg, 1); % line form
  LC = cell(nalg, 1); % line color
  MF = cell(nalg, 1); % marker form
  MC = cell(nalg, 1); % marker color
  MS = cell(nalg, 1); % marker size
  
  for i = 1:nalg
    alg = algs{i};
    
    if strcmp( alg, 'oica-gencov' )
      LF{i} = '-';
      LC{i} = c.dviolet;
      MF{i} = '*';
      MC{i} = c.red;
      MS{i} = 17;
    end
    if strcmp( alg, 'oica-gencov-clust' )
      LF{i} = '--';
      LC{i} = c.magenta;
      MF{i} = '^';
      MC{i} = c.lgreen;
      MS{i} = 17;
    end
    if strcmp( alg, 'oica-gencov-ada' )
      LF{i} = ':';
      LC{i} = c.yellow;
      MF{i} = 'o';
      MC{i} = c.lblue;
      MS{i} = 17;
    end
    if strcmp( alg, 'oica-gencov-semiada' )
      LF{i} = '-';
      LC{i} = c.red;
      MF{i} = '*';
      MC{i} = c.dviolet;
      MS{i} = 17;
    end
    if strcmp( alg, 'oica-semiada' )
      LF{i} = '-';
      LC{i} = c.red;
      MF{i} = '*';
      MC{i} = c.dviolet;
      MS{i} = 17;
    end
    if strcmp( alg, 'oica' )
      LF{i} = '-';
      LC{i} = c.red;
      MF{i} = '*';
      MC{i} = c.dviolet;
      MS{i} = 17;
    end
    if strcmp( alg, 'oica-quad-semiada-h' )
      LF{i} = '-';
      LC{i} = c.red;
      MF{i} = '*';
      MC{i} = c.dviolet;
      MS{i} = 17;
    end
    if strcmp( alg, 'oica-quad-semiada' )
      LF{i} = '-';
      LC{i} = c.lblue;
      MF{i} = '*';
      MC{i} = c.dviolet;
      MS{i} = 17;
    end
    if strcmp( alg, 'oica-quad-clust' )
      LF{i} = '--';
      LC{i} = c.cian;
      MF{i} = '^';
      MC{i} = c.lgreen;
      MS{i} = 17;
    end
    if strcmp( alg, 'oica-quad-ada' )
      LF{i} = ':';
      LC{i} = c.blue;
      MF{i} = 'o';
      MC{i} = c.lblue;
      MS{i} = 17;
    end
    if strcmp( alg, 'foobi' )
      LF{i} = '-';
      LC{i} = c.green;
      MF{i} = '^';
      MC{i} = c.orange;
      MS{i} = 17;
    end
    if strcmp( alg, 'rica' )
      LF{i} = '-';
      LC{i} = c.lorange;
      MF{i} = 'o';
      MC{i} = c.green;
      MS{i} = 17;
    end
    if strcmp( alg, 'fpca' )
      LF{i} = '-';
      LC{i} = c.blue;
      MF{i} = 'd';
      MC{i} = c.lpink;
      MS{i} = 17;
    end
    if strcmp( alg, 'rand' )
      LF{i} = ':';
      LC{i} = c.lorange;
      MF{i} = '*';
      MC{i} = c.lorange;
      MS{i} = 17;
    end
      
    
  end
  
  for i = 1:nalg
    if isempty(LF{i})
      error('No specs for one of the algs')
    end
  end
  
end


function c = get_colors

  c.('orange')  = [0.8 0.3 0.0];
  c.('lorange') = [0.9 0.5 0.0];
  c.('pink')    = [1.0 0.0 0.5];
  c.('lpink')   = [1.0 0.5 0.5];
  c.('blue')    = [0.0 0.2 0.6];
  c.('lblue')   = [0.5 0.5 1.0];
  c.('green')   = [0.1 0.4 0.1];
  c.('lgreen')  = [0.2 0.5 0.2];
  c.('dviolet') = [0.5 0.0 0.5];
  c.('cian')    = [0.0 1.0 1.0];
  c.('dcian')   = [0.0 0.5 0.5];
  c.('red')     = [1.0 0.0 0.0];
  c.('dgreen')  = [0.0 1.0 0.0];
  c.('magenta') = [1.0 0.0 1.0];
  c.('yellow')  = [0.9 0.9 0.1];

end



