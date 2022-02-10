function varargout=gps2rng(files,meth,xyz,v)
% [st,dxyz,sr,xyz,v]=GPS2RNG(files,meth,xyz,v)
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
% st           the travel time between the points and the target
% dxyz         the GNSS time series, whichever the summary of the files
% sr           the distance between the timeseries of points and the target
% xyz          1x3 matrix with nominal coordinates of the target, in com
% v            sound speed [m]
%
% EXAMPLE:
%
% gps2rng({'0001-05340.mat','0002-05340.mat','0003-05340.mat','0004-05340.mat'})
% gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'})
% gps2rng({'Unit2-camp.mat'})
%
% Originally written by tschuh-at-princeton.edu, 11/24/2021
% Last modified by tschuh-at-princeton.edu 02/10/2022
% Last modified by fjsimons-at-princeton.edu 02/10/2022

% need this file: IFILES/TOPOGRAPHY/EARTH/GEBCO/GEBCO2014/GEBCO_2014_1D.nc
% if you dont have it: xyz=[2e6 -4.5e6 3e6]

% how the possibly multiple receivers will be jointly considered
defval('meth','ave')
% constant sound speed profile for now [m/s]
defval('v',1500)

% Get extension from first filename
[~,fname,~]=fileparts(files{1});
filex=suf(fname,'-');
% If more than one input, list them all
if length(files)>1
  for index=1:length(files)
    [~,fname,~]=fileparts(files{index});
    nr(index)=indeks(pref(fname,'-'),'end');
  end
  fname=sprintf('Unit%s-%s',nr(1:end),suf(pref(files{1},'.'),'-'));
end

% combine all datasets into 1 large matrix with no plotting
[d,tmax]=mat2mod(files);

% Tick marks
ntix=7;
ttix=[1 round([1:ntix-1]*length(d(1).t)/(ntix-1))];
tixl=datestr(d(1).t(ttix),'HH:MM:SS');

if length(files)>1
  switch meth
   case 'ave'
    % More or less one line version of mat2com.m
    dxyz=squeeze(nanmean(reshape(cat(1,d(:).xyz),size(d(1).xyz,1),length(d),3),2));
  end
else
  % Save yourself the trouble
  dxyz=d.xyz;
end

% Default beacon location is empty
defval('xyz',[])
% If empty, use either gps2dep or gps2guess to get a location
if isempty(xyz)
    % Determine a starting point from the middle section
    imeth='down'; 
    switch imeth
     case 'down'
      % Put in water depth read off the GPSTRAJECT map
      depth=5225; 
      % Just about
      % lonlat=[291.3001853867659   31.4661752724632];
      % z=gebco(lonlat(1)-360,lonlat(2));
      % Handpicked by GINPUT on a plot of dxyz... for DOG1
      [dogx,dogy,dogz]=gps2dep([1.977967 -5.073198 3.3101016]*1e6,depth);
     case 'guess'
      % Likely water depth from PrincetonSeafloorGeodesy-SURVEY3.pdf ORIGIN
      defval('depth',gebco(68+42/60,-(31+27/60)));
      % Rather put in a better water depth, so all units return uniform
      % results, and read off the GPSTRAJECT map
      depth=5225;
      msex=[1 3]/4;
      midsex=[round(msex(1)*size(dxyz,1)):round(msex(2)*size(dxyz,1))];
      [dogx,dogy,dogz]=gps2guess(dxyz(midsex,:),depth);
    end
    xyz=[dogx dogy dogz];
% If nonempty, use it and define dogx, dogy, dogz as such
else
     dogx = xyz(1);
     dogy = xyz(2);
     dogz = xyz(3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the forward model
% calculate slant range between ship and beacon for each second in meters
sr=sqrt((dxyz(:,1)-dogx).^2+(dxyz(:,2)-dogy).^2+(dxyz(:,3)-dogz).^2);
% calculate slant time from slant range and v [s]
% currently not doing anything with this
st = sr./v;
% In the proper version, it's through complicated velocity, and refraction
% [sr,st]=raytracingofsomesort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% optional output
varns={st,dxyz,sr,xyz,v};
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
  % Be sure to mark the minimum
  
  % Here you'll see that the total time represented could be wider than any
  % single one, since you've done MAT2MOD
  tl(1)=title(sprintf('Raw Slant Range Measurements: %s-%s',...
		      datestr(tmax(1)),datestr(tmax(2))));
  grid on
  longticks([],2)
  
  xticks(d(1).t(ttix))
  xticklabels(tixl)

  % Arbitrarily absolute, ok for the moment
  yl=ylim;
  topy=yl(2)-[yl(2)-yl(1)]/20;
  tt(1)=text(tmax(1),topy,sprintf('%i%s NaN',...
				   round(sum(isnan(sr))/length(sr)*100),'%'),...
	     'VerticalAlignment','top');

  stran='';
  for index=1:length(d)
    stran=sprintf('%s%s\n',stran,pref(files{index},'-'));
  end
  stran=stran(1:end-1);
  tt(2)=text(tmax(2),topy,stran,'HorizontalAlignment','right',...
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
