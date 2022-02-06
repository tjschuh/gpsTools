function gps2his(files,protype)
% GPS2HIS(files,protype)
%
% Given Precise Point Position time series of four different units, computes
% their pairwise distances, calculates a least-squares regression line and
% produces histograms of residuals (with quantile-quantile plots).
%
% INPUT:
% 
% files        cell with MAT-filename strings containing data structures
% protype      type of prd file ('ppp' or 'rtk')
%
% EXAMPLE:
%
% gps2his({'0001-05340.mat','0002-05340.mat','0003-05340.mat','0004-05340.mat'})
% gps2his({'0001-F089.mat','0002-F089.mat','0003-F089.mat','0004-F089.mat'})
%
% Originally written by tschuh-at-princeton.edu, 12/01/202
% Last modified by tschuh-at-princeton.edu, 02/03/2022
% Last modified by fjsimons-at-princeton.edu, 02/06/2022

% New output filename made from first, you'll save the full info
[~,fname,~] = fileparts(files{1});
fname=sprintf('000X-%s.mat',suf(fname,'-'));

switch protype
 case 'ppp'
  % light blue --> PPP
  defval('col',[0.400 0.6667 0.8431])
 case 'rtk'
  % lime green --> RTK
  defval('col',[0.466 0.6740 0.1880])
end

% keep rows where nsats > nthresh and pdop < pthres and pdop~=0
nthresh = 4; pthresh = 15;
% outlier removal by percentile
percs=[1 99];

if exist(fname)~=2
  % convert data to all be same time spans with no time gaps
  [d,tmax] = mat2mod(files);
  % compute pairwise Euclidean distances between receivers
  nk=nchoosek(1:length(d),2);

  for k=1:size(nk,1)
    i=nk(k,1); j=nk(k,2);
    
    % Remember the times etc were prematched
    dest{k} = sqrt([d(i).xyz(:,1)-d(j).xyz(:,1)].^2 + ...
		   [d(i).xyz(:,2)-d(j).xyz(:,2)].^2 + ...
		   [d(i).xyz(:,3)-d(j).xyz(:,3)].^2);
    % keep the original length
    dn(k)=length(dest{k});
    % find the good data condition
    cond=d(i).pdop<pthresh & d(i).pdop~=0 & d(i).nsats(:,1)>nthresh & ...
	 d(j).pdop<pthresh & d(j).pdop~=0 & d(j).nsats(:,1)>nthresh;
    % Calculate the residuals of the linear fit, applying condition
    thetimes=seconds(d(i).t(cond)-d(i).t(cond(1)));

    % May start an inner loop  - something recursive - "while"
    p{k}=polyfit(thetimes,dest{k}(cond),1);
    % Calculate residuals as observed minus prediction
    e{k}=dest{k}(cond)-polyval(p{k},thetimes);
    % remove outliers to get better results
    try
      [ee{k},er{k}]=rmoutliers(e{k},'gesd');
      em{k}='gesd';
    catch
      [ee{k},er{k}]=rmoutliers(e{k},'percentiles',percs);
      em{k}='percentiles';
    end
  end
  % Might enter into loop again... to redo fits after outlier removal

  % Save whatever you need, all still in standard units
  save(fname,'e','p','ee','em','percs','nthresh','pthresh','nk','tmax','dn')
else
  load(fname)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting in case you have exactly 4 files, hence 6 distances
% plot histograms with normal curves over top
% make histograms by using histcounts and bar
% make curves by using fitdist and pdf
f=figure(1); clf
[ah,ha]=krijetem(subnum(3,2));
f.Position=[250 500 1100 600];

% Convert to from SI in m to mm
ucon=1000;

for k=1:length(ah)
  % Unit conversion
  e{k}=e{k}*ucon;
  ee{k}=ee{k}*ucon;

  axes(ah(k))
  % nbins is computed for each dataset using Freedman-Diaconis' Rule
  % this is a personal choice over other methods including Sturges' Rule
  % experimented with other methods and this was my preference
  nbins=round((max(ee{k})-min(ee{k}))/(2*iqr(ee{k})*(length(ee{k}))^(-1/3)));

  % calculate goodness of fit (gof) compared to a normal distribution
  [~,~,stats]=chi2gof(ee{k},'NBins',nbins);

  % divide chi squared by degrees of freedom to reduce to 1 DoF
  % with 1 Dof, chi squared <= 4 signifies ~90% chance data are normal
  % will make red curve dotted if gof > 4
  gof=stats.chi2stat/stats.df;
  % Calculate histogram
  [yvals,edges]=histcounts(ee{k},nbins);
  % Calculate bin centers 
  barc{k}=0.5*(edges(1:end-1)+edges(2:end));
  % Plot the histogram
  b{k}=bar(barc{k},yvals,'BarWidth',1);
  % Cosmetics
  str1=sprintf('%i-%i, # of points = %i/%i/%i',nk(k,1),nk(k,2),...
			 length(ee{k}),length(e{k}),dn(k));
  str2=sprintf('%i-%i',nk(k,1),nk(k,2));
  [lain{k},xl(k),yl(k)]=cosmo1(ah(k),...
		 str2,...
		 'residuals [mm]','counts',ee{k},gof,b{k},nthresh,pthresh,col);
end

% finishing touches - you should keep minmax times from before
tt=supertit(ah([1 2]),sprintf('GNSS Pairwise Distance Hourly Residuals %s to %s',...
			       datestr(tmax(1)),datestr(tmax(2))));
movev(tt,0.3)
%set(tt,'FontName','Courier')
% Centering is off for non-fixed-width fonts
moveh(tt,0.07775)
delete(xl(1:4))
delete(yl([2 4 6]))

moveh(ah([1 3 5]),.05)

% It is smart enough to strip the extension but we do oit anyway
figdisp([],pref(fname,'.'),'',2,[],'epstopdf')

xver=0;
if xver==1
  % Not for now, we don't

  % plot qq plots for each dataset to measure how "normal" data is
  g=figure;
  g.Position = [250 500 1100 600];

  ah2(1) = subplot(3,2,1);
  qq{k} = qqplot(ee12);
  cosmo2('GPS Pair 1-2',qq12)

  % finishing touches
  tt=supertit(ah1([1 2]),sprintf('QQ Plots of Residuals vs Standard Normals (%s to %s)',...
				 datestr(d1.t(1)),datestr(d1.t(end))));
  movev(tt,0.3)

  %figdisp(sprintf('qqplot-%s',fname),[],'',2,[],'epstopdf')
end


% cosmetics for histogram and pdf plots
function [lain,xl,yl] = cosmo1(ax,titl,xlab,ylab,data,gof,b,nthresh,pthresh,col)
 
ax.XGrid = 'on';
ax.YGrid = 'off';
ax.GridColor = [0 0 0];
ax.TickLength = [0 0];

xl=xlabel(xlab);

% Establish a scale
stdd=std(data);
nstd=3;

% later versions could use XTICKS and XLABELS 
xel=round([-nstd nstd]*stdd,2);
xlim(xel)
ax.XTick=round([-nstd:nstd]*stdd,2);
% Trick to do it properly
bla=round([-nstd:nstd]*stdd,0);
for in=1:length(bla)
  blace{in}=sprintf('%.0f',bla(in));
end
ax.XTickLabel=blace;
yl=ylabel(ylab);
yel=[0 max(b.YData)+0.1*max(b.YData)];
ylim(yel)
longticks(ax,2)
ax.YTick=unique(0:50:max(yel));
if length(ax.YTick)>6
  ax.YTick=unique(0:100:max(yel));
end

% If you move the title, can't do more adjustments of axes
tto=title(titl);
tto.Position=tto.Position.*[0 1 1];
moveh(tto,stdd/30);

% Quote data stats
mrg1=0.95;
mrg2=0.995;
mrg3=0.95;
fs=9;
t(1)=text(mrg2*nstd*stdd,mrg3*ax.YLim(2),...
     sprintf('N = %3.0f\nstd = %3.0f\nmed = %3.0f\navg = %3.0f\ngof = %3.0f',...
	     length(data),stdd,median(data),mean(data),gof));

% This is calculating the percent of the data shown due to axis limitation
pct=(length(data(data<=xel(2) & data>=xel(1)))/length(data))*100;

t(2)=text(-mrg1*nstd*stdd,mrg3*ax.YLim(2),...
     sprintf('%05.2f%%\nmin = %.0f\nmax = %.0f\ndop < %.0f\nsat > %.0f',...
	     pct,min(data),max(data),pthresh,nthresh));

% Forget LISTFONTS
set(t(1),'HorizontalAlignment','right')
set(t,'VerticalAlignment','top','FontSize',fs);

% Overlayes
% light blue bars
b.FaceColor =col;
% plot vertical line at median
hold on
% later versions could use XLINE (flilr on ylim?)
plot([1 1]*median(data),ylim,'k--','LineWidth',1);

% plot pdf on top of histogram
pd = fitdist(data,'Normal');
xvals = [-nstd*stdd:nstd*stdd];
yvals = pdf(pd,xvals);
area = sum(b.YData)*diff(b.XData(1:2));
lain = plot(xvals,yvals*area,'r','LineWidth',2);
% if gof > 4 make red curve dotted b/c data is not normal
if gof > 4
  lain.LineStyle = '--';
end
hold off

% cosmetics for qq plots
function cosmo2(titl,qq)
grid on
longticks([],2)
title(titl)
ylabel('Residual Quantiles')
%xlim([-4 4])
%ylim([-100 100])
qq(1).Marker = '.';
qq(1).MarkerSize = 8;
qq(3).LineWidth = 2;
qq(3).LineStyle = '--';
