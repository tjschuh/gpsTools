function varargout=kinpro(prdfile,xyzpos,plt,demean)
% d=KINPRO(prdfile,xyzpos,plt,demean)
%
% Makes/saves/loads a *.prd GNSS Precise Point Positioning solution file
% (created by PPP2PRD or RTK2PRD) into a structured *.mat file (and plots),
% for subsequent joint handling by, e.g., MAT2MOD.
%
% INPUT:
%
% prdfile      string with filename (e.g. 'kin_28-29_pton.prd')
% xyzpos       approximate position xyz of station copied from RINEX header
% plt          0 for no plot, 1 for plot (default: 1)
% demean       0 for no demeaning, 1 for demeaning ENU data (default: 0)
%
% OUTPUT:
%
% d            actual data struct
% .mat file    output file saved as mat file to working directory
% .txt file    output ascii file saved to working directory
%
% EXAMPLE:
%
% d=kinpro('kin_28-29_pton.prd',[1288235.7406 -4694422.9216 4107355.8820],1,0);
% kinpro('kin_2018230-2018231_raul.prd',[-5566082.7400 -201282.8434 -3097636.2460],1,0)
%
% TESTED ON:
%
% R2020a Update 4 (9.8.0.1417392)
%
% Originally written by tschuh-at-princeton.edu, 03/04/2022
% Last modified by tschuh-at-princeton.edu, 03/29/2022

% Default file name
defval('prdfile','kin_28-29_pton.prd')

% New output filename made from input
[~,fname,~]=fileparts(prdfile);
outfile1=sprintf('%s.mat',fname);
outfile2=sprintf('%s.txt',fname);

% Make a plot or not, 1 --> yes plot
defval('plt',1)

% Demean data or not, 0 --> dont demean
defval('demean',0)

% if outfile doesnt exist, make it and save it, otherwise load it
if exist(outfile1,'file')==0
    % station information from RINEX file
    % approx position [m]
    
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
               'InputFormat','dd-MMM-yyyy HH:mm:ss');
    % convert local time to UTC
    % dont need to convert to UTC b/c RINEX files are in GPS time aka UTC
    %t.TimeZone = 'Z';

    % convert x,y,z to east-west,north-south,up-down (ENU)
    % using wgs84 and station reference location
    wgs84 = wgs84Ellipsoid('meter');
    [lat0,lon0,h0] = ecef2geodetic(wgs84,xyzpos(1),xyzpos(2),xyzpos(3));
    [E,N,U] = ecef2enu(dm(:,3),dm(:,4),dm(:,5),lat0,lon0,h0,wgs84);

    % get rid of sat cols that are all zeros
    % cols 9:15 are sat info
    satcol=9:15;
    chuck=~sum(dm(:,satcol));
    hsat=sattypes(~chuck);
    sats=dm(:,satcol(~chuck));

    % make data structure explicitly
    scale{1,1} = 1e0; scale{1,2} = 'm';
    d.(h{1}) = t;
    d.timezone = 'UTC';
    d.(h{2}) = [E N U]*scale{1,1};
    d.ENUunit = scale{1,2};
    d.(h{3}) = sats;
    d.(h{4}) = dm(:,end);
    % save as mat file
    save(outfile1,'d')
    % save as ascii file
    % may want to save variable names and units
    T = table(t,E,N,U);
    writetable(T,outfile2,'WriteVariableNames',0,'Delimiter',' ')
else
    % Load
    disp(sprintf('Loading %s',outfile1))
    load(outfile1)
end

% plotting
if plt == 1
    % if user wants, subtract mean from data
    if demean == 1
        d.enu(:,1) = d.enu(:,1)-mean(d.enu(:,1));
        d.enu(:,2) = d.enu(:,2)-mean(d.enu(:,2));
        d.enu(:,3) = d.enu(:,3)-mean(d.enu(:,3));
    end
        
    f=figure;
    f.Position=[70 120 560 685];

    ah(1)=subplot(4,1,1);
    plot(d.t,d.enu(:,1),'k')
    ylabel('east [m]')
    cosmoenu(d.enu(:,1),d.t)
    
    ah(2)=subplot(4,1,2);
    plot(d.t,d.enu(:,2),'k')
    ylabel('north [m]')
    cosmoenu(d.enu(:,2),d.t)
    
    ah(3)=subplot(4,1,3);
    plot(d.t,d.enu(:,3),'k')
    ylabel('up [m]')
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

    tt=supertit(ah(1),sprintf('Station Approx Position: [%.4f %.4f %.4f] m',xyzpos(1),xyzpos(2),xyzpos(3)));
movev(tt,0.1)
    
    figdisp(fname,[],[],2,[],'epstopdf')
end
    
% optional output
varns={d};
varargout=varns(1:nargout);

function cosmoenu(data,time)
xlim([time(1) time(end)])
slack=1.1;
multi = 0.01;
ylim([min(data)-multi*abs(min(data)) max(data)+multi*abs(max(data))])
longticks([],2)
grid on
xticklabels([])