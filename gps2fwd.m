function varargout = gps2fwd(d,tmax,xyz,v)
% [st,xyz,v] = GPS2FW(d,tmax,xyz,v)
%
%
%
% INPUT:
%
% d            single dataset with x,y,z,t ship locations
% tmax         last time in d
% xyz          1x3 matrix with nominal coordinates of the target, in com [m]
% v            sound speed [m/s]
%
% OUTPUT:
%
% st           the travel time between the points and the target
% xyz          1x3 matrix with nominal coordinates of the target, in com [m]
% v            sound speed [m/s]
%
% EXAMPLE:
%
% load Unit1234-camp.mat
% [st,xyz,v] = gps2fwd(d,tmax,[2e6 -4.5e6 3e6],1500);
%
% Originally written by tschuh-at-princeton.edu, 02/16/2022

% choose a beacon location at center of our real "cross" trajectory [m]
defval('xyz',[2e6 -4.5e6 3e6])

% constant sound speed profile for now [m/s]
defval('v',1500)

% This is the forward model
% calculate slant range between ship and beacon for each second in meters
sr = sqrt((d.xyz(:,1)-xyz(1)).^2 + (d.xyz(:,2)-xyz(2)).^2 + (d.xyz(:,3)-xyz(3)).^2);
% calculate slant time from slant range and v [s]
st = sr./v;
% In the proper version, it's through complicated velocity, and refraction
% [sr,st]=raytracingofsomesort

% optional output
varns={st,xyz,v};
varargout=varns(1:nargout);

% Make a plot if you don't want output
if nargout==0
  % plot the distance [km] vs time
  % would be nice to eventually add error bars
  figure(1); clf
  ah=plot(d.t,sr,'LineWidth',2);
  xlim(tmax)
  % open it up a smidgen
  xel=xlim; yel=ylim;
  xlim(xel+[-1 1]*range(xel)/20);
  ylim(yel+[-1 1]*range(yel)/20);
  ylabel('slant range [m]')
  % Be sure to mark the minimum
  
  % Here you'll see that the total time represented could be wider than any
  % single one, since you've done MAT2MOD
  tl(1)=title(sprintf('Raw Slant Range Measurements: %s-%s',...
		      datestr(tmax(1)),datestr(tmax(2))));
  grid on
  longticks([],2)

  % Tick marks
  ntix=7;
  ttix=[1 round([1:ntix-1]*length(d(1).t)/(ntix-1))];
  tixl=datestr(d(1).t(ttix),'HH:MM:SS');
  
  xticks(d.t(ttix))
  xticklabels(tixl)

  % Arbitrarily absolute, ok for the moment
  yl=ylim;
  topy=yl(2)-[yl(2)-yl(1)]/20;
  tt(1)=text(tmax(1),topy,sprintf('%i%s NaN',...
				   round(sum(isnan(sr))/length(sr)*100),'%'),...
	     'VerticalAlignment','top');

  stran='';
  % for index=1:length(d)
  %   stran=sprintf('%s%s\n',stran,pref(files{index},'-'));
  % end
  stran=stran(1:end-1);
  tt(2)=text(tmax(2),topy,stran,'HorizontalAlignment','right',...
	     'VerticalAlignment','top');
  set(tt(:),'FontSize',8);

  % Give it an XLABEL with the day!
  dat1=datestr(d.t(1),'dd mmm yyyy');
  dat2=datestr(d.t(end),'dd mmm yyyy');
  if dat1==dat2
    xlabel(sprintf('%s',dat1))
  else
    xlabel(sprintf('%s to %s',dat1,dat2))
  end
  
  % Delete title
  delete(tl)
  
  % save figure as pdf
  figdisp([],[],[],2,[],'epstopdf')
end
