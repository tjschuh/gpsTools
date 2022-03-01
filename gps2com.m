function varargout=gps2com(files,plt)
% dd=GPS2COM(files,plt)
%
% Combine a series of GNSS time series into one dataset via averaging
% (eventually will be center of mass com) and then plot the results in
% similar fashion as prd2mat.m
%
% INPUT:
%
% files         cell with MAT-filename strings containing data structures
% plt           0 for no plot, 1 for plot (default: 1)
%
% OUTPUT:
%
% dd             actual data struct
%
% EXAMPLE:
%
% dd = gps2com({'0001-05340.mat','0002-05340.mat','0003-05340.mat','0004-05340.mat'},1)
%
% Originally written by tschuh-at-princeton.edu, 11/23/2021
% Last modified by tschuh-at-princeton.edu, 03/01/2022

[d,~] = mat2mod(files);

[~,fname,~] = fileparts(files{1});
fname = sprintf('000Z-%s.mat',suf(fname,'-'));


d e f i n i t e l y n o t r e a d y y e t

% Unlike GPS2RNG which computes a distance wrt to a point, this is just a
% time series. So its product is like those from PRD2MAT and could just
% use GPS2PLT to render them.


% combine all data into 1 dataset
% simply by taking the average of the data
% maybe will eventually use a different scheme
% ignoring d4 for now because bad data
dd.t = d(1).t;
alldx = [d(1).xyz(:,1) d(2).xyz(:,1) d(3).xyz(:,1) d(4).xyz(:,1)];
dd.x = nanmean(alldx,2);
alldy = [d(1).xyz(:,2) d(2).xyz(:,2) d(3).xyz(:,2) d(4).xyz(:,2)];
dd.y = nanmean(alldy,2);
alldz = [d(1).xyz(:,3) d(2).xyz(:,3) d(3).xyz(:,3) d(4).xyz(:,3)];
dd.z = nanmean(alldz,2);
alldlat = [d(1).lat d(2).lat d(3).lat d(4).lat];
dd.lat = nanmean(alldlat,2);
alldlon = [d(1).lon d(2).lon d(3).lon d(4).lon];
dd.lon = nanmean(alldlon,2);
alldutme = [d(1).utmeasting d(2).utmeasting d(3).utmeasting d(4).utmeasting];
dd.utme = nanmean(alldutme,2);
alldutmn = [d(1).utmnorthing d(2).utmnorthing d(3).utmnorthing d(4).utmnorthing];
dd.utmn = nanmean(alldutmn,2);
dd.utmz = d(1).utmzone;
alldht = [d(1).height d(2).height d(3).height d(4).height];
dd.ht = nanmean(alldht,2);
alldnsats = [d(1).nsats(:,1) d(2).nsats(:,1) d(3).nsats(:,1) d(4).nsats(:,1)];
dd.nsats = nanmean(alldnsats,2);
alldpdop = [d(1).pdop d(2).pdop d(3).pdop d(4).pdop];
dd.pdop = nanmean(alldpdop,2);

% plotting
defval('plt',1)
if plt == 1
    f=figure;
    % position = [left bottom width height]
    f.Position = [500 250 850 550];

    % find rows where nsats <= 4
    nthresh = 4;
    n = dd.nsats;

    % also should find rows where pdop is >= 10 or = 0
    % this doesnt always coincide with low nsats
    pthresh = 15;
    p = dd.pdop;

    % plotting interval
    int = 5;

    % plot utm coordinates
    % set the zero for the UTM coordinates based on the min and max of data
    x = dd.utme-(min(dd.utme)-.05*(max(dd.utme)-min(dd.utme)));
    y = dd.utmn-(min(dd.utmn)-.05*(max(dd.utmn)-min(dd.utmn)));
    tc = datetime(dd.t,'Format','HH:mm:ss'); 

    % find good (g) and bad (b) data
    % [gx bx] = x
    gx = x; bx = x;
    gy = y; by = y;
    gx(p>=pthresh | p==0 | n<=nthresh) = NaN;
    bx(p<pthresh & n>nthresh) = NaN;
    gy(p>=pthresh | p==0 | n<=nthresh) = NaN;
    by(p<pthresh & n>nthresh) = NaN;

    ah(1)=subplot(2,2,[1 3]);
    c = linspace(1,10,length(x(1:int:end)));
    sz = 10;
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
    xlabel('Easting [m]')
    %xticklabels({'0','5','10','15'})
    ylabel('Northing [m]')
    %yticklabels({'0','2','4','6','8','10','12','14','16','18','20'})
    title(sprintf('Ship Location (Every %dth point)',int))

    % plot heights relative to WGS84
    gh = dd.ht; bh = dd.ht;
    gh(p>=pthresh | p==0 | n<=nthresh) = NaN;
    bh(p<pthresh & n>nthresh) = NaN;

    ah(2)=subplot(2,2,2);
    plot(dd.t(1:int:end),gh(1:int:end),'color',[0.4660 0.6740 0.1880])
    hold on
    % grey out "bad" data where nsats is too low or pdop is too high or 0
    plot(dd.t(1:int:end),bh(1:int:end),'color',[0.7 0.7 0.7])
    xlim([dd.t(1) dd.t(end)])
    xticklabels([])
    % remove outliers so plotting looks better
    htout = rmoutliers(dd.ht,'mean');
    outpct = (length(dd.ht)-length(htout))*100/length(dd.ht);
    ylim([min(htout,[],'all')-0.005*abs(min(htout,[],'all')) max(htout,[],'all')+0.005*abs(max(htout,[],'all'))])
    a=annotation('textbox',[0.78 0.57 0 0],'String',[sprintf('%05.2f%% Outliers',outpct)],'FitBoxToText','on');
    a.FontSize = 8;
    grid on
    longticks
    ylabel('Height relative to WGS84 [m]')
    title(sprintf('Ship Height (Every %dth Point)',int))

    % plot nsats and pdop on same plot
    ah(3)=subplot(2,2,4);
    yyaxis left
    plot(dd.t,dd.nsats(:,1),'b','LineWidth',1)
    yticks([min(dd.nsats(:,1))-1:max(dd.nsats(:,1))+1])
    ylim([min(dd.nsats(:,1))-0.5 max(dd.nsats(:,1))+0.5])
    ylabel('Number of Observed Satellites') 
    yyaxis right
    plot(dd.t,dd.pdop,'r','LineWidth',1)
    ylim([min(dd.pdop)-0.25 max(dd.pdop)+0.25])
    xlim([dd.t(1) dd.t(end)])
    ylabel('Position Dilution Of Precision')
    % can only turn grid on for left axis
    grid on
    longticks
    title('Total Number of Satellites and PDOP')

    tt=supertit(ah([1 2]),sprintf('1 Hour of Averaged Ship Data Starting from %s',datestr(dd.t(1))));
    movev(tt,0.3)

    a = annotation('textbox',[0.23 0.1 0 0],'String',['camp'],'FitBoxToText','on');
    a.FontSize = 12;

    figdisp(fname,[],'',2,[],'epstopdf')

    close
end

% optional output
varns={dd};
varargout=varns(1:nargout);
