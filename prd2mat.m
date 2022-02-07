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

switch protype
 case 'ppp'
  % light blue --> PPP
  defval('col',[0.400 0.6667 0.8431])
 case 'rtk'
  % lime green --> RTK
  defval('col',[0.466 0.6740 0.1880])
end

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

  % Symbol size
  sz=10;

  % plotting interval
  pint=20;
  deltat=seconds(d.t(pint+1)-d.t(1));
  
  % Average speed
  vave=nanmean(sqrt([diff(d.utmeasting).^2+diff(d.utmnorthing).^2])./...
	       seconds(diff(d.t)))*3600/1000;

  % subsample first, then look for condition
  dkp={'t','xyz','lat','lon','utmeasting','utmnorthing','height','nsats','pdop'};
  for index=1:length(dkp)
    % Don't forget that some are multicolumn, and keep the last one
    ssam=unique([1:pint:length(d.(dkp{index})) length(d.(dkp{index}))]);
    d.(dkp{index})=d.(dkp{index})(ssam,:);
  end
  
  % find the good data condition
  cond=d.pdop<pthresh & d.pdop~=0 & d.nsats(:,1)>nthresh;

  % Just open a figure and fiddle with positions outide
  figure(1)
  clf
  % plot utm coordinates referenced to the minima
  x=d.utmeasting -min(d.utmeasting);
  y=d.utmnorthing-min(d.utmnorthing);
  % and refer the times to the beginning of the sequence
  t=d.t-d.t(1);
  % Tick marks
  ntix=3;
  ttix=[1 round([1:ntix-1]*length(d(1).t)/(ntix-1))];
  tixl=datestr(d.t(ttix),'HH:MM:SS');

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
  % First the good data
  sg=scatter(gx',gy',sz,seconds(t),'filled'); hold on
  % Then the bad data
  sb=scatter(bx',by',sz,grey,'filled'); hold off
  xlabel('relative easting [m]')
  ylabel('relative northing [m]')
  tl(1)=title(sprintf('Ship Location (Every %dth Point)',pint));
  grid on; box on
  axis equal tight
  % open it up a smidgen
  xel=xlim; yel=ylim;
  ah(1).XLim=xel+[-1 1]*range(xel)/20;
  ah(1).YLim=yel+[-1 1]*range(yel)/20;
  % Some minor annotations
  tx(1)=text(0,0,sprintf('dop < %.0f\nsat > %.0f\n%s = %3.1f km/h\n%st = %i s',...
			 pthresh,nthresh,'v',vave,'\Delta',deltat),...
	     'VerticalAlignment','bottom');


  longticks(ah)
  % Then the color scale
  colormap(jet)
  cb=colorbar('southoutside','Ticks',seconds(t(ttix)),'TickLabels',tixl);
  %shrink(cb,1,1.5)
  longticks(cb)
  % Give the color bar an xlabel with the day!
  dat1=datestr(d.t(1),'dd mmm yyyy');
  dat2=datestr(d.t(end),'dd mmm yyyy');
  if dat1==dat2
    cb.XLabel.String=sprintf('%s',dat1);
  else
    cb.XLabel.String=sprintf('%s to %s',dat1,dat2);
  end

  % Second panel - the elevation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ah(2)=subplot(2,2,2);
  % First the good data
  plot(d.t,gh,'color',col); hold on
  % Then the bad data
  plot(d.t,bh,'color',grey); hold off
  xlim([d.t(1) d.t(end)])
  xticks(d.t(ttix))
  xticklabels(tixl)
  
  % remove outliers so plotting looks better
  % htout = rmoutliers(d.height,'mean');
  % How much was removed, in percent
  % outpct = (length(d.height)-length(htout))/length(d.height)*100;
  % ylim([min(htout,[],'all')-0.005*abs(min(htout,[],'all')) max(htout,[],'all')+0.005*abs(max(htout,[],'all'))])
  grid on
  longticks
  ylabel('height above WGS84 [m]')
  tl(2)=title(sprintf('Ship Height (Every %dth Point)',pint));
  
  % Third panel, nsats and pdop  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ah(3)=subplot(2,2,4);
  yyaxis left
  ps=plot(d.t,d.nsats(:,1),'b');
  yticks([min(d.nsats(:,1))-1:max(d.nsats(:,1))+1])
  ylim([min(d.nsats(:,1))-0.5 max(d.nsats(:,1))+0.5])
  % ylabel('number of observed satellites')
  ylabel('number of satellites')
  yyaxis right
  pp=plot(d.t,d.pdop,'r');
  ylim([min(d.pdop)-0.25 max(d.pdop)+0.25])
  xlim([d.t(1) d.t(end)])
  xticks(d.t(ttix))
  xticklabels(tixl)
  %ylabel('position dilution of precision')
  ylabel('dilution of precision')
  % can only turn grid on for left axis
  grid on
  longticks
  tl(3)=title('Total Number of Satellites and PDOP');

  tl(4)=supertit(ah([1 2]),sprintf('Ship Data from %s to %s',d.t(1),d.t(end)));
  movev(tl(4),0.3)

  % Get rid of the ugly titles
  delete(tl)
  % Print to file
  figdisp([],fname,[],2)
end

% optional output
varns={d};
varargout=varns(1:nargout);
