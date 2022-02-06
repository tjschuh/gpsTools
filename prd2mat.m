function varargout=prd2mat(prdfile,protype,plt)
% d=PRD2MAT(prdfile,protype,plt)
%
% Turns a *.prd GNSS Precise Point Positioning solution file (created by
% PPP2PRD or RTK2PRD) into a structured *.mat file (and a plot)
%
% INPUT:
%
% prdfile      string with filename (e.g. '0002-05340.prd')
% protype      type of prd file ('ppp' or 'rtk')
% plt          0 for no plot, 1 for plot (default: 1)
%
% OUTPUT:
%
% d            actual data struct
% .mat file    output file saved as mat file to working directory
%
% EXAMPLE
%
% d=prd2mat('0002-05340.prd','ppp',1);
%
% REQUIRES:
%
% DEG2UTM.M
%
% TESTED ON:
%
% R2020a Update 4 (9.8.0.1417392)
% 9.0.0.341360 (R2016a) - without the timetabling... 
%
% Originally written by tschuh-at-princeton.edu, 10/06/2021
% Last modified by tschuh-at-princeton.edu, 02/03/2022
% Last modified by fjsimons-at-princeton.edu, 02/06/2022

% TO DO:
% extract trip section from fname somehow
% currently manually changing at end of code
%
% signify on plot when utmzone changes

% Default file name and type
defval('prdfile','0002-05340.prd')
defval('protype','ppp')

% New output filename made from input
[~,fname,~]=fileparts(prdfile);
outfile=sprintf('%s.mat',fname);

% Make a plot or not
defval('plt',1)

if protype=='ppp'
  % if outfile doesnt exist, make it and save it, otherwise load it
  if exist(outfile,'file')==0
    % explicitly define header from *.kin
    h = {'t','xyz','lat','lon','height','nsats','pdop'};
    % all possible satellite types
    sattypes = {'Total','GPS','GLONASS','Galileo','BDS-2','BDS-3','QZSS'};
    
    % original *.kin file columns were: 
    % Mjd, SoD, X, Y, Z, Latitude, Longitufe, Height, Nsats, PDOP
    dm=load(prdfile);

    % make datetime array from *.kin columns 1 and 2
    % convert Mjd to ymd
    ymd=datestr(dm(:,1)+678942);
    % convert SoD to hms
    hms=datestr(seconds(dm(:,2)),'HH:MM:SS');
    % convert tstr to datetime
    t=datetime([ymd repmat(' ',size(ymd,1),1) hms],...
		 'InputFormat','dd-MMM-yyyy HH:mm:ss');

    % convert lat & lon to utm easting & northing in meters
    warning off MATLAB:nargchk:deprecated
    % lat lon cols are 6 and 7
    [x,y,zone]=deg2utm(dm(:,6)',rem((dm(:,7))',180)-180);
    % if unique zone save only one
    if ~sum(sum(zone,1)/length(zone)-zone(1,:)); zone=zone(1,:); end
    warning on MATLAB:nargchk:deprecated
        
    % get rid of sat cols that are all zeros
    % cols 9:15 are sat info
    satcol=9:15;
    chuck=~sum(dm(:,satcol));
    hsat=sattypes(~chuck);
    sats=dm(:,satcol(~chuck))
    
    % make data structure explicitly
    d.(h{1}) = t;
    d.(h{2}) = dm(:,3:5);
    d.xyzunit = 'm';
    d.(h{3}) = dm(:,6);
    d.(h{4}) = dm(:,7);
    d.lonlatunit = 'deg';
    d.utmeasting = x;
    d.utmnorthing = y;
    d.utmunit = 'm';
    d.utmzone = zone;
    d.(h{5}) = dm(:,8);
    d.heightunit = 'm (rel to WGS84)';
    d.satlabels = hsat;
    d.(h{6}) = sats;
    d.(h{7}) = dm(:,end);
    % Save
    save(outfile,'d')
  else
    % Load
    disp(sprintf('Loading %s',outfile))
    load(outfile)
  end

elseif protype == 'rtk'
  % if outfile doesnt exist, make it and save it
  % otherwise load it
  if exist(outfile,'file') == 0
    
    % load data using readtable rather than load bc
    % columns 1 and 2 have / and : which dont work with load
    % need the 'MultipleDelimAsOne' to load columns in correct order
    dm = readtable(prdfile,'FileType','text','MultipleDelimsAsOne',true);

    % make datetime array from prdfile columns 1 and 2
    % convert column 1 from a cell array to datetime array
    col1=datetime(dm.Var1,'InputFormat','yyyy/MM/dd');
    % combine col1 and dm.Var2 into a datetime array
    dt = col1 + dm.Var2;

    % convert lat,lon,height to x,y,z in meters
    % specify reference ellipsoid as WGS84
    wgs84 = wgs84Ellipsoid('meter');
    % use built-in matlab function geodetic2ecef to do conversion
    [x,y,z] = geodetic2ecef(wgs84,dm.Var3,dm.Var4,dm.Var5);
    
    % convert lat & lon to utm easting & northing in meters
    warning off MATLAB:nargchk:deprecated
    % lat lon cols are 6 and 7
    [utmx,utmy,zone]=deg2utm(dm.Var3(:)',dm.Var4(:)');
    % if unique zone save only one
    if ~sum(sum(zone,1)/length(zone)-zone(1,:)); zone=zone(1,:); end
    warning on MATLAB:nargchk:deprecated
        

    % make data structure explicitly
    d.t = dt;
    d.xyz = [x y z];
    d.xyzunit = 'm';
    d.lat = dm.Var3;
    d.lon = dm.Var4;
    d.lonlatunit = 'deg';
    d.utmeasting = utmx;
    d.utmnorthing = utmy;
    d.utmunit = 'm';
    d.utmzone = zones;
    d.height = dm.Var5;
    d.heightunit = 'm (rel to WGS84)';
    d.satlabels = 'Total';
    d.nsats = dm.Var7;
    d.pdop = dm.Var15;
    % Save
    save(outfile,'d')
  else
    load(outfile)
  end

else
  error('Please select a valid protype')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting
if plt == 1

  % before anything else, need to fill timeskips with NaNs we don't do this
  % during the creation of the mat file because we want the mat file to be
  % consistent with the prd file we also want to be able to combine prd/mat
  % files and this method will still work here
  
  % Use (RE)TIME(TABLE) to fill in time skips/data gaps with NaNs
  drem={'xyzunit','lonlatunit','utmunit','heightunit','satlabels','utmzone'};
  d=retimes(d,drem,{'secondly','fillwithmissing'});
  
  % keep rows where nsats > nthresh and pdop < pthres and pdop~=0
  nthresh = 4; pthresh = 15;

  % find the good data condition
  cond=d.pdop<pthresh & d.pdop~=0 & d.nsats(:,1)>nthresh;

  % plotting interval
  int = 10;
  % Symbol size
  sz = 10;
  
  % Just open a figure and fiddle with positions outide
  figure(1)
  clf
  
  % plot utm coordinates
  % set the zero for the UTM coordinates based on the min and max of data
  x = d.utmeasting-(min(d.utmeasting)-.05*(max(d.utmeasting)-min(d.utmeasting)));
  y = d.utmnorthing-(min(d.utmnorthing)-.05*(max(d.utmnorthing)-min(d.utmnorthing)));
  tc = datetime(d.t,'Format','HH:mm:ss'); 

  % Only good data
  gx=x; gx(~cond)=NaN;
  gy=y; gy(~cond)=NaN;
  gh=d.height; gh(~cond)=NaN;
  % Only bad data
  bx=x; bx(cond)=NaN;
  by=y; by(cond)=NaN;
  bh=d.height; bh(cond)=NaN;

  % First panel - the ship track %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ah(1)=subplot(2,2,[1 3]);
  c = linspace(1,10,length(x(1:int:end)));
  scatter(gx(1:int:end)',gy(1:int:end)',sz,c,'filled')
  colormap(jet)
  colorbar('southoutside','Ticks',[1:3:10],'TickLabels',...
	   {datestr(tc(1),'HH:MM:SS'),datestr(tc(floor(end/3)),'HH:MM:SS'),...
	    datestr(tc(ceil(2*end/3)),'HH:MM:SS'),datestr(tc(end),'HH:MM:SS')})
  hold on
  % grey out "bad" data where nsats is too low or pdop is too high or 0
  scatter(bx(1:int:end)',by(1:int:end)',sz,[0.7 0.7 0.7],'filled')
  grid on
  longticks
  xlabel('easting [m]')
  ylabel('northing [m]')
  t(1)=title(sprintf('Ship Location (Every %dth Point)',int));

  % Second panel - the elevation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ah(2)=subplot(2,2,2);
  plot(d.t(1:int:end),gh(1:int:end),'color',[0.4660 0.6740 0.1880])
  hold on
  % grey out "bad" data where nsats is too low or pdop is too high or 0
  plot(d.t(1:int:end),bh(1:int:end),'color',[0.7 0.7 0.7])
  xlim([d.t(1) d.t(end)])
  xticklabels([])
  % remove outliers so plotting looks better
  htout = rmoutliers(d.height,'mean');
  outpct = (length(d.height)-length(htout))*100/length(d.height);
  ylim([min(htout,[],'all')-0.005*abs(min(htout,[],'all')) max(htout,[],'all')+0.005*abs(max(htout,[],'all'))])
  grid on
  longticks
  ylabel('height above WGS84 [m]')
  t(2)=title(sprintf('Ship Height (Every %dth Point)',int));
  
  % Third panel, nsats and pdop  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ah(3)=subplot(2,2,4);
  yyaxis left
  plot(d.t,d.nsats(:,1),'b','LineWidth',1)
  yticks([min(d.nsats(:,1))-1:max(d.nsats(:,1))+1])
  ylim([min(d.nsats(:,1))-0.5 max(d.nsats(:,1))+0.5])
  ylabel('number of observed satellites')
  yyaxis right
  plot(d.t,d.pdop,'r','LineWidth',1)
  ylim([min(d.pdop)-0.25 max(d.pdop)+0.25])
  xlim([d.t(1) d.t(end)])
  ylabel('position dilution of precision')
  % can only turn grid on for left axis
  grid on
  longticks
  t(3)=title('Total Number of Satellites and PDOP');
  
  % Get rid of the titles
  delete(t)
  
  keyboard
  
  

  tt=supertit(ah([1 2]),sprintf('Ship Data from %s to %s',datestr(d.t(1)),datestr(d.t(end))));
  movev(tt,0.3)

  a = annotation('textbox',[0.23 0.1 0 0],'String',['Unit 1: camp'],'FitBoxToText','on');
  a.FontSize = 12;
  
  figdisp([],fname,[],2)

end

% optional output
varns={d};
varargout=varns(1:nargout);
