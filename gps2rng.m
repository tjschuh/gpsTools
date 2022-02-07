function varargout=gps2rng(files,meth,xyz,v)
% [sr,st]=GPS2RNG(files,meth,xyz,v)
%
% Given Precise Point Position time series of four different units, compute
% a new average ship position time series from the four units and then
% compute the slant range to a fixed, imaginary beacon on the seafloor
%
% INPUT:
%
% files        cell with MAT-filename strings containing data structures
%              of which the field 'xyz' will be used, in com coordinates
% meth         averaging method, 'ave'
% xyz          1x3 matrix with nominal coordinates of the target, in com
% v            sound speed [m]
%
% OUTPUT:
%
% sr           the distance between the timeseries of points and the target
% st           the travel time between the points and the target
%
% EXAMPLE:
%
% gps2rng({'0001-05340.mat','0002-05340.mat','0003-05340.mat','0004-05340.mat'})
% gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'})
% gps2rng({'Unit2.camp.mat'})
%
% Originally written by tschuh-at-princeton.edu, 11/24/2021
% Last modified by tschuh-at-princeton.edu 02/03/2022
% Last modified by tschuh-at-princeton.edu 02/07/2022

% inspired by file 05040 first xyz positions, but ~5 km below surface
% eventually will need to be a changing value until it reaches seafloor
% after that will not be known and will be deduced from timestamps sent back
defval('xyz',[1.979e6 -5.074e6 3.30385e6])

defval('meth','ave')
% constant sound speed profile for now [m/s]
defval('v',1500)

% set beacon location on seafloor
dogx=xyz(1);
dogy=xyz(2);
dogz=xyz(3);

% Get extension from first filename
[~,fname,~]=fileparts(files{1});
filex=suf(fname,'-');

% combine all datasets into 1 with no plotting
[d,tmax]=mat2mod(files);

% Tick marks
ntix=7;
ttix=[1 round([1:ntix-1]*length(d(1).t)/(ntix-1))];
tixl=datestr(d(1).t(ttix),'HH:MM:SS');

% Averaging method, or keeping them invidual, or taking only one,
switch meth
  case 'ave'
   dxyz=squeeze(nanmean(reshape(cat(1,d(:).xyz),size(d(1).xyz,1),length(d),3),2));
end

% calculate slant range between ship and beacon for each second in meters
sr=sqrt((dxyz(:,1)-dogx).^2+(dxyz(:,2)-dogy).^2+(dxyz(:,3)-dogz).^2);

% calculate slant time from slant range and v [s]
% currently not doing anything with this
st = sr./v;

% optional output
varns={sr,st};
varargout=varns(1:nargout);

% Make a plot if you don't want output
if nargout==0
  % plot the distance [km] vs time
  % would be nice to eventually add error bars
  figure(1); clf
  ah=subplot(2,1,1);
  plot(d(1).t,sr,'LineWidth',2)
  xlim(tmax)
  % open it up a smidgen
  xel=xlim; yel=ylim;
  ah(1).XLim=xel+[-1 1]*range(xel)/20;
  ah(1).YLim=yel+[-1 1]*range(yel)/20;
  ylabel('slant range [m]')
  % Here you'll see that the total time represented could be wider than any
  % single one, since you've done MAT2MOD
  tl(1)=title(sprintf('Raw Slant Range Measurements: %s-%s',...
		      datestr(tmax(1)),datestr(tmax(2))));
  grid on
  longticks([],2)
  
  xticks(d(1).t(ttix))
  xticklabels(tixl)

  % Arbitrarily absolute, ok for the moment
  tt(1)=text(tmax(1),15000,sprintf('%i%s NaN',...
				   round(sum(isnan(sr))/length(sr)*100),'%'),...
	     'VerticalAlignment','top');

  stran='';
  for index=1:length(d)
    stran=sprintf('%s%s\n',stran,pref(files{index},'-'));
  end
  stran=stran(1:end-1);
  tt(2)=text(tmax(2),15000,stran,'HorizontalAlignment','right',...
	     'VerticalAlignment','top');
  set(tt(:),'FontSize',8);

  % Give it an XLABEL with the day!
  dat1=datestr(d(1).t(1),'dd mmm yyyy');
  dat2=datestr(d(1).t(end),'dd mmm yyyy');
  if dat1==dat2
    ah(1).XLabel.String=sprintf('%s',dat1);
  else
    ah(1).XLabel.String=sprintf('%s to %s',dat1,dat2);
  end
  
  % Delete title
  delete(tl)
  
  % save figure as pdf
  figdisp([],fname,[],2,[],'epstopdf')
end

