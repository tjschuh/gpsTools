function varargout=kinpro(prdfile,plt)
% d=PRD2MAT(prdfile,protype,plt)
%
% Makes/saves/loads a *.prd GNSS Precise Point Positioning solution file
% (created by PPP2PRD or RTK2PRD) into a structured *.mat file (and plots),
% for subsequent joint handling by, e.g., MAT2MOD.
%
% INPUT:
%
% prdfile      string with filename (e.g. '0002-05340.prd')
% plt          0 for no plot, 1 for plot (default: 1)
%
% OUTPUT:
%
% d            actual data struct
% .mat file    output file saved as mat file to working directory
%
% EXAMPLE:
%
% d=prd2mat('0002-05340.prd',1);
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
% Originally written by tschuh-at-princeton.edu, 03/04/2022

% TO DO:
% extract trip section from fname somehow
% currently manually changing at end of code
%
% signify on plot when utmzone changes

% Default file name
%defval('prdfile','0002-05340.prd')

% New output filename made from input
[~,fname,~]=fileparts(prdfile);
outfile=sprintf('%s.mat',fname);

% Make a plot or not
defval('plt',1)

% if outfile doesnt exist, make it and save it, otherwise load it
if exist(outfile,'file')==0
    % station information from RINEX file
    % approx position [m]
    x0 = 1288235.7406;
    y0 = -4694422.9216;
    z0 = 4107355.8820;
    % local timezone
    local = 'America/New_York';
    
    % explicitly define header from *.kin
    h = {'t','enu','nsats','pdop'};
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
               'InputFormat','dd-MMM-yyyy HH:mm:ss','TimeZone',local);
    % convert EST local time to UTC
    t.TimeZone = 'Z';

    % convert x,y,z to east-west,north-south,up-down (ENU)
    % using wgs84 and station reference location
    wgs84 = wgs84Ellipsoid('meter');
    [lat0,lon0,h0] = ecef2geodetic(wgs84,x0,y0,z0);
    [E,N,U] = ecef2enu(dm(:,3),dm(:,4),dm(:,5),lat0,lon0,h0,wgs84);
    % subtract mean from data b/c thats what seismologists do
    E = E-mean(E);
    N = N-mean(N);
    U = U-mean(U);

    % get rid of sat cols that are all zeros
    % cols 9:15 are sat info
    satcol=9:15;
    chuck=~sum(dm(:,satcol));
    hsat=sattypes(~chuck);
    sats=dm(:,satcol(~chuck));

    % make data structure explicitly
    scale{1,1} = 1e2; scale{1,2} = 'cm';
    d.(h{1}) = t;
    d.timezone = 'UTC';
    d.(h{2}) = [E N U]*scale{1,1};
    d.ENUunit = scale{1,2};
    d.(h{3}) = sats;
    d.(h{4}) = dm(:,end);
    % Save
    save(outfile,'d')
else
    % Load
    disp(sprintf('Loading %s',outfile))
    load(outfile)
end

% plotting
if plt == 1
    f=figure;
    f.Position=[70 120 560 685];
    ah(1)=subplot(4,1,1);
    plot(d.t,d.enu(:,1),'k')
    ylabel('east [cm]')
    cosmoenu(d.enu(:,1),d.t)
    
    ah(2)=subplot(4,1,2);
    plot(d.t,d.enu(:,2),'k')
    ylabel('north [cm]')
    cosmoenu(d.enu(:,2),d.t)
    
    ah(3)=subplot(4,1,3);
    plot(d.t,d.enu(:,3),'k')
    ylabel('up [cm]')
    cosmoenu(d.enu(:,3),d.t)
    
    ah(4)=subplot(4,1,4);
    yyaxis left
    ps=plot(d.t,d.nsats(:,1),'b');
    %ntix=8;
    %ttix=[1 round([1:ntix-1]*length(d(1).t)/(ntix-1))];
    %tixl=datestr(d.t(ttix),'HH:MM:SS');
    %yticks([min(d.nsats(:,1))-1:max(d.nsats(:,1))+1])
    ylim([min(d.nsats(:,1))-0.5 max(d.nsats(:,1))+0.5])
    ylabel('number of satellites')
    yyaxis right
    pp=plot(d.t,d.pdop,'r');
    %ylim([min(d.pdop)-0.25 max(d.pdop)+0.25])
    xlim([d.t(1) d.t(end)])
    %xticks(d.t(ttix))
    %xticklabels(tixl)
    xlabel('time [UTC]')
    ylabel('dilution of precision')
    % can only turn grid on for left axis
    grid on
    longticks([],2)
    %tl(4)=title('Total Number of Satellites and PDOP');

    figdisp(fname,[],[],2,[],'epstopdf')
end
    
% optional output
varns={d};
varargout=varns(1:nargout);

function cosmoenu(data,time)
xlim([time(1) time(end)])
slack=1.1;
ylim(slack*[min(data) max(data)])
longticks([],2)
grid on
xticklabels([])