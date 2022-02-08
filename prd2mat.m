function varargout=prd2mat(prdfile,protype,plt)
% d=PRD2MAT(prdfile,protype,plt)
%
% Makes/saves/loads a *.prd GNSS Precise Point Positioning solution file
% (created by PPP2PRD or RTK2PRD) into a structured *.mat file (and plots),
% for subsequent joint handling by, e.g., MAT2MOD.
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
% EXAMPLE:
%
% d=prd2mat('0002-05340.prd','ppp',1);
%
% REQUIRES:
%
% DEG2UTM.M
%
% SEE ALSO:
%
% MAT2MOD
%
% TESTED ON:
%
% R2020a Update 4 (9.8.0.1417392)
% 9.0.0.341360 (R2016a) - without the timetabling... 
%
% Originally written by tschuh-at-princeton.edu, 10/06/2021
% Last modified by tschuh-at-princeton.edu, 02/03/2022
% Last modified by fjsimons-at-princeton.edu, 02/08/2022

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
  % Make a plot. Inside GPS2PLT we fill timeskips with NaNs. We don't do
  % this here during the creation of the mat file because we want the mat
  % file to be maximally consistent with the prd file, and consider
  % alternatives.
  % Name override
  ovr=mfilename;
  % Color schame
  switch protype
   case 'ppp'
    % light blue --> PPP
    defval('col',[0.400 0.6667 0.8431])
   case 'rtk'
    % lime green --> RTK
    defval('col',[0.466 0.6740 0.1880])
  end
  gps2plt(d,col,ovr,fname);
end

% optional output
varns={d};
varargout=varns(1:nargout);
