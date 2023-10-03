 function [ p, sp, func_str ] = linfitxy( xdatain, ydatain, xerrin, yerrin, varargin)

% [P] = LINFITXY (XDATA, YDATA, 0 ,0) returns estimation of P=[A B] so that
% ydata is best estimated by A * XDATA + B , with uncertainty of unknown
% magnitude on YDATA.
%
% [P,SP] = LINFITXY (XDATA, YDATA, 0 , 0) returns also estimated
% uncertainty (i.e. one standard deviation) of fitted parameters A and B 
% in SP, with uncertainty of unknown magnitude on YDATA.
%
% [P,SP] = LINFITXY(XDATA, YDATA, XERR, YERR) returns estimation of P=[A B]
% with known uncertainty on XDATA and on YDATA (resp. XERR and YERR). It is
% either a scalar or vector of length XDATA (or YDATA), containing one
% standard deviation associated to each measurement.
%
% The uncertainty of the fitted parameters is computed using Monte Carlo
% simulation with NBLOOP=500 iterations and assuming errors are
% Gaussian and centered, and given with NSIGMA times standard deviation. By
% default, NSIGMA equals to 1.
%
% [P,SP] = LINFITXY(...,'NbLoop', 1500, 'NSigma',2.5) will specify the
% numbers of loops to compute and the number of standard deviation for
% error parameter estimates.
%
% By default a plot of the data and errorbars is displayed, together with
% the nonsimultaneous prediction bounds for the linear fit function.
% Confidence level is given considering NSIGMA * standard deviations.
%
% Finally the result of the fit is displayed is a self-explanatory format.
%
% Other options : 
% - 'Verbosity' : default 1. No output 0. Greater verbosity
% means more output informations.
% - 'Plotting' : default true. By default a plot of the data is done, with
%                a fitting line.
% - 'PlotHull' : default true. By default, when it's possible, the envelope
%                of all possible line fits is plotted
% - 'OnlyPlot' : default false. If wanted, it is possible to plot the data 
%                without performing a linear fit. 
% - 'FitResultOnGraph' : default true. If true, put the fit result as the 
%                title of the graph.
%
% NOTA: Least square computations are done using FMINSEARCH of the standard
% matlab toolbox. 
%
% EXAMPLE 1 (assume no error on x and unknown error on y)
%   xdata = (0:10);
%   ydata = 1.2 - 3.4*xdata + 0.5*randn(1,length(xdata));
%   fittedparam=linfitxy(xdata, ydata,0,0);
%
% EXAMPLE 1-bis (idem as 1 but to get error on fitted param)
%   xdata = (0:10);
%   ydata = 1.2 - 3.4*xdata + 0.5*randn(1,length(xdata));
%   [fittedparam,fittederrorparam]=linfitxy(xdata, ydata,0,0)
%
% EXAMPLE 2 (with error estimation)
%   xdata = (0:10);
%   ydata = 1.2 - 3.4*xdata + 2*randn(1,length(xdata));
%   linfitxy(xdata, ydata,0.2,3.,'Verbosity',2);
%
% EXAMPLE 3 (with variable error estimation)
%   xdata = (0:10);
%   ydata = 0 + 1*xdata + [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1] ...
%                         .*randn(1,length(xdata));
%   xerr  = 0 ;
%   yerr  = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1] ;
%   linfitxy(xdata, ydata,xerr,yerr,'Verbosity',2, 'NbLoop', 1000, 'NSigma',2);
% 
%SEE ALSO fminsearch, optimset
%
% COMPATIBILITY : Standard Matlab 2010A / no additional library needed
%                 or Octave 3.4.0 with optim 1.0.16 package
%
%
%   - -  Paris Diderot feb 2014                             - -
%        Contacts : tristan.beau@univ-paris-diderot.fr , 
%                   julien.browaeys@univ-paris-diderot.fr
%
%        File available at : 
%          http://www.mathworks.fr/matlabcentral/fileexchange/45711-linear-fit-with-both-errors-in-x-and-in-y
%

%---- Default values ----

nhull=100;
nbloop=500;
nsigma=1;
verbosity=1;
plotting=true;
plothull=true;
onlyplot=false;
fitresultongraph=true;

%---- Test if matlab or octave mode
matlabmode=1;
cur_ver=ver('MATLAB');
if isempty(cur_ver)
  matlabmode=0;
end

%---- Manage hold status
holdstatus=ishold;

%---- Checking Arguments of the function ----

if nargin < 4
  error('Not enough arguments. Please read help');
end

nOptArgin = length(varargin); %Number of arguments input

if nOptArgin == 0
  disp('')
elseif nOptArgin == 1
  error('Need either no error arguments or 2 error arguments');
elseif mod(nOptArgin,2) ~= 0
  error('Bad optional argument statement');
else
  for k=1:2:nOptArgin
    switch varargin{k}
      case 'NbLoop'
        nbloop = varargin{k+1};
        if nbloop<1
          error('NbLoop must be a positive integer');
        end
      case 'NSigma'
        nsigma = varargin{k+1};
        if nsigma<=0
          error('NSigma must be a positive number');
        end
      case 'Verbosity'
        verbosity = varargin{k+1};
        if verbosity < 0
          error('Verbosity must be positive');
        end
      case 'Plotting'
        plotting=varargin{k+1};
      case 'PlotHull'
        plothull=varargin{k+1};
      case 'OnlyPlot'
        onlyplot=varargin{k+1};
        plotting=true;
        plothull=false;
        verbosity=0;
      case 'FitResultOnGraph'
        fitresultongraph=varargin{k+1};
      otherwise
        error('Unknown %s argument',varargin{k});
    end
  end
end

% Impose line vector 
xdata=xdatain(:)';
ydata=ydatain(:)';
xerr=xerrin(:)';
yerr=yerrin(:)';

% some tests on inputs
if length(xdata)~=length(ydata)
  error('XDATA and YDATA must be of the same length');
end
if ( length(xdata)<=1 )
  error('DATA should contain more than one value');
end

if (length(yerr)==1) && (yerr == 0)
  if verbosity
    warning('no error given on yerr. Yerr is assumed unknown and deduced from data and model')
  end
  nbloop=0;
else
  if sum((xerr.^2 + yerr.^2) == 0) > 1
    error('More than one point with null uncertainty, no fit possible');
  end
end
  
if ( length(xerr)==1 )
  xerr = xerr * ones(1,length(xdata));
end
if ( length(yerr)==1 )
  yerr = yerr * ones(1,length(ydata));
end  
if length(xerr) ~= length(xdata)
  error('XDATA and XERR must be of the same length');
end
if length(yerr) ~= length(ydata)
  error('YDATA and YERR must be of the same length');
end

% one more test to check if errors are infinite : restriction of data
wherefinite = isfinite(xerr) & isfinite(yerr); % finds where xerr AND yerr are both finite
xerr=xerr(wherefinite) ; % using logical subscriting, removes datapoints where error is infinite
xdata=xdata(wherefinite);
yerr=yerr(wherefinite) ;
ydata=ydata(wherefinite);
if any(not(wherefinite)) % display a warning and some explanation
  warning('warn:infiniteerr','Some errors given were infinite. \n The corresponding datapoints are removed, \n since information content associated to them is null.');
end


% ---- Fitting if NO error estimation of the parameters is given ---- %
 
[PolyCoef,S] = polyfit(xdata,ydata,1);
p=PolyCoef;

if ( abs(p(2))/max(abs(ydata)) < 1e-3 ) 
  p(2)=0;
end % This will avoid the fminsearch bug when p(2) is negligeable

invR=inv(S.R);
covMatrix= (invR*invR')*S.normr^2/S.df;
sp=nsigma*sqrt(covMatrix([1,4]));

if nbloop==0
   if plotting
     if plothull
       warning('no hull available yet with no y error given')
     end
     plot(xdata,ydata,'x')
     hold on      
     xplot = min(xdata)+(max(xdata)-min(xdata)).*(0-0.2:1./(nhull-1):1.+0.2);
     if not(onlyplot)
       plot([min(xplot),max(xplot)],p(1)*[min(xplot),max(xplot)]+p(2),'r')
     end
     
% if error not input then calculate yerr from data and model

     yerr = S.normr/sqrt(length(xdata)-2) ;
     
     if not(onlyplot)
       disp(['Estimated error on Y is : ',num2str(yerr)])
     end
     
     if verbosity
       disp(['Linear fit function : Y = (', ...
         num2str(p(1)),' +/- ',num2str(sp(1)), ') * X + (', ...
         num2str(p(2)),' +/- ',num2str(sp(2)),')'])
       func_str = ['Linear fit function : Y = (', ...
         num2str(p(1)),' +/- ',num2str(sp(1)), ') * X + (', ...
         num2str(p(2)),' +/- ',num2str(sp(2)),')'];
     end
     
     for i=1:length(xdata)
         plot([xdata(i) xdata(i)],[ydata(i)-yerr ydata(i)+yerr],'-k')
     end   
     if not(holdstatus) % restore hold status
         hold off
     end
       return
   else
       return
   end
end

% ---- Fitting WHEN error estimation over XDATA and YDATA are given ---- %

if not(onlyplot)
% some inits
if matlabmode
  options = optimset('MaxFunEvals',10000); % to avoid optimset warnings (?)
end
pLoop = zeros(nbloop,2); % to allocate the full matrix

% use previous fit as initial parameters

if verbosity > 1 
  disp 'Computing...'
end

% test if there is a fixed point (no xerr neither yerr on it)
bFixedPoint=0;
if ( sum((xerr.^2 + yerr.^2) == 0) == 1 )
  bFixedPoint=1;
  % then fit should be over the slope only
  idxFixedPoint = find((xerr.^2 + yerr.^2) == 0);
  xFixed = xdata(idxFixedPoint);
  yFixed = ydata(idxFixedPoint);
  %Now remove fixed point
  xdataNoFixed = xdata;
  ydataNoFixed = ydata;
  xerrNoFixed =xerr;
  yerrNoFixed = yerr;
  
  xdataNoFixed(idxFixedPoint) = [];
  ydataNoFixed(idxFixedPoint) = [];
  xerrNoFixed(idxFixedPoint) = [];
  yerrNoFixed(idxFixedPoint) = [];
end

for i = 1 : nbloop   % i.e. do not enter in that loop if nbloop == 0
  
  if bFixedPoint == 1
      % Fit with y-y0 = a[1] (x-x0)
    x=xdataNoFixed+randn(1,length(xdataNoFixed)).*xerrNoFixed;
    y=ydataNoFixed+randn(1,length(ydataNoFixed)).*yerrNoFixed;
    
    chi2 = @(a) sum( ( (y-yFixed)-a(1).*(x-xFixed)).^2 ./ ( yerrNoFixed.^2 + a(1)^2*xerrNoFixed.^2 ) );
    
    if matlabmode
       pLoop(i,1) = fminsearch(chi2,[p(1)],options);
     else
       pLoop(i,1) = fminsearch(chi2,[p(1)]);
    end
    pLoop(i,2) = yFixed - pLoop(i,1)*xFixed ; 
     
  else
    % Fit with y = ax+b == a[1]* x + a[2]
    x=xdata+randn(1,length(xdata)).*xerr;
    y=ydata+randn(1,length(ydata)).*yerr;
    
    chi2 = @(a) sum( (y-a(1).*x-a(2)).^2 ./ ( yerr.^2 + a(1)^2*xerr.^2 ) );
    
     if matlabmode
       pLoop(i,:) = fminsearch(chi2,[p(1),p(2)],options);
     else
       pLoop(i,:) = fminsearch(chi2,[p(1),p(2)]);
     end
  
  end
  
  %%%%%% REM de JB
  % peut-etre faut-il se servir du fait que l'ordonnee est determinee par
  % la pente pour resoudre une equation a une inconnue avec fminbnd
  % aussi on peut penser a regarder les dernieres pages du poly de M2 qui
  % donne une equation de degre 3 ? resoudre pour tout savoir
  
  if verbosity > 1 && mod(i,100)==0
    if i>100
      for ix=1:22
        fprintf('\b')
      end
    end
    fprintf('  ... %04.1f %% completed',i*100./nbloop)
    if matlabmode==0
      fflush(stdout);
    end
  end
end

if verbosity > 1 
    disp ' '
end

% first version
p=mean(pLoop);
sp=std(pLoop)*nsigma;

% MC version
threshold=0.5*erfc(nsigma/sqrt(2));
for k=1:2
  sorted_pt_k=sort(pLoop(:,k)); 
  low_index_k=floor(nbloop*threshold);
  high_index_k=ceil(nbloop*(1-threshold));
  plow(k)=sorted_pt_k(1+low_index_k); 
  phigh(k)=sorted_pt_k(high_index_k);
end
p=(plow+phigh)/2;
sp=phigh-p;

end %onlyplot)

% Now plotting
if plotting
  
  plot(xdata,ydata,'kx')

  hold on
  for i=1:length(xdata)
    x0=xdata(i);
    x1=x0-xerr(i);
    x2=x0+xerr(i);
    y0=ydata(i);
    y1=y0-yerr(i);
    y2=y0+yerr(i);
    plot([x1 x2],[y0 y0],'-k')
    plot([x0 x0],[y1 y2],'-k')
  end
  
  if not(onlyplot)
  xplot = min(xdata)+(max(xdata)-min(xdata)).*(0-0.2:1./(nhull-1):1.+0.2);  
  plot([min(xplot),max(xplot)],p(1)*[min(xplot),max(xplot)]+p(2),'r')
  
  if plothull
    axipbi = pLoop(:,1)*xplot+pLoop(:,2)*ones(1,length(xplot));
    saxipbi = sort(axipbi,1);
    low_alpha=floor(nbloop*threshold);
    high_alpha=ceil(nbloop*(1-threshold));
    if (low_alpha == 0 || high_alpha==nbloop) 
      low_alpha = 1;
      high_alpha = nbloop;
      warning('warn:hullpbl','NbLoop too low to draw the correct Hull Plot. \n Hull Plot drawn with extremal computed values only.');
    end
    plot(xplot,saxipbi(low_alpha,:),'r--')
    plot(xplot,saxipbi(high_alpha,:),'r--')
  end
  
  end %onlyplot
  
  if not(holdstatus) % restore hold status 
    hold off
  end

  FitResultString=['Y = (', ...
         num2str(p(1)),' +/- ',num2str(sp(1)), ') * X + (', ...
         num2str(p(2)),' +/- ',num2str(sp(2)),')'] ;
  
  if verbosity
    disp(FitResultString);
  end

  if plotting && fitresultongraph
    title(FitResultString)
  end
end
