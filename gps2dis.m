function gps2dis(files,protype)
% GPS2DIS(files,protype)
%
% Given Precise Point Position time series of four different units,
% computes their pairwise distances, and plots them
%
% INPUT:
% 
% files        cell with MAT-filename strings containing data structures
% protype      type of prd file ('ppp' or 'rtk')
%
% EXAMPLE:
%
% gps2dis({'0001-05340.mat','0002-05340.mat','0003-05340.mat','0004-05340.mat'})
%
% Originally written by tschuh-at-princeton.edu, 11/12/2021
% Last modified by tschuh-at-princeton.edu, 06/06/2022

% to-do:
% fix annotations/text boxes outside axes
% add rtk vs ppp differences (protype not used currently)

% use mat2mod to convert data to all be same time spans with no time gaps
[d,tmax] = mat2mod(files);
[~,fname,~] = fileparts(files{1});

% keep rows where nsats > nthresh and pdop < pthresh
nthresh = 4; pthresh = 15;

% plotting intervals
pint1 = 30;
pint2 = 5;
pint3 = 15;

figure(1)
clf
cols = {'r','g','b','k'};

for i=1:length(d)
    % plot heights of all 4 GPS receivers
    % find the good data condition
    cond1=d(i).pdop<pthresh & d(i).pdop~=0 & d(i).nsats(:,1)>nthresh;
    % [g b] = h
    % only good data
    gh=d(i).height; gh(~cond1)=NaN;
    % only bad data
    bh=d(i).height; bh(cond1)=NaN;
    % keep all the data for the end
    allht(:,i) = gh;

    % also compute velocity components vx, vy, vz in cm/s
    vx = 100*diff(d(i).xyz(:,1))./seconds(diff(d(i).t));
    vy = 100*diff(d(i).xyz(:,2))./seconds(diff(d(i).t));
    vz = 100*diff(d(i).xyz(:,3))./seconds(diff(d(i).t));
    allvx(:,i) = vx;
    allvy(:,i) = vy;
    
    % actually plot heights of all 4 recs
    ah(1)=subplot(5,2,[1 3]);
    plot(d(i).t(1:pint1:end),gh(1:pint1:end),cols{i})
    hold on
    plot(d(i).t(1:pint1:end),bh(1:pint1:end),'color',[0.7 0.7 0.7])
    if i == length(d)
        % compute average spatial and time velocity in knots (1 knot = 1852 m/hr)
        vxavg = nanmean(allvx,2);
        vyavg = nanmean(allvy,2);
        vxymag = sqrt(vxavg.^2 + vyavg.^2);
        vavg = nanmean(vxymag,1);
        vavg = vavg*3600*1e-2/1852;

        % cosmetics for subplot 1
        xlim([d(i).t(1) d(i).t(end)])
        xticklabels([])
        ylabel('Height relative to WGS84 [m]')
        t1=title(sprintf('Ship Height (Every %dth Point)',pint1));
        grid on
        longticks
        % to set best ylim, remove outliers from alldht
        % then find the global min
        try
            allhtout = rmoutliers(allht,'gesd');
        catch
            allhtout = rmoutliers(allht,'mean');
        end
        outpct = (length(allht)-length(allhtout))*100/length(allht);
        ylim([min(allhtout,[],'all')-0.005*abs(min(allhtout,[],'all')) ...
              max(allhtout,[],'all')+0.005*abs(max(allhtout,[],'all'))])
        text(d(i).t(floor(0.03*length(d(i).t))),ah(1).YLim(1)+0.0025*abs(ah(1).YLim(1)),...
             sprintf('%05.2f%% Outliers',outpct),'FontSize',9)
        text(d(i).t(floor(0.475*length(d(i).t))),ah(1).YLim(1)+0.0025*abs(ah(1).YLim(1)),...
             sprintf('Nsats > %d & PDOP < %d',nthresh,pthresh),'FontSize',9)
        text(d(i).t(floor(0.65*length(d(i).t))),ah(1).YLim(2)-0.0025*abs(ah(1).YLim(2)),...
             sprintf('v = %.2f knots',vavg),'FontSize',9)
    end

    % plot GPS receiver acceleration components
    % compute acceleration components ax, ay, az in cm/s^2
    ax = diff(vx)./(2*seconds(diff(d(i).t(1:end-1))));
    ay = diff(vy)./(2*seconds(diff(d(i).t(1:end-1))));
    az = diff(vz)./(2*seconds(diff(d(i).t(1:end-1))));

    cond2=d(i).pdop(1:end-2)<pthresh & d(i).pdop(1:end-2)~=0 & d(i).nsats(1:end-2,1)>nthresh;
    
    % good data
    gax=ax; gax(~cond2)=NaN;
    gay=ay; gay(~cond2)=NaN;
    gaz=az; gaz(~cond2)=NaN;
    % only bad data
    bax=ax; bax(cond2)=NaN;
    bay=ay; bay(cond2)=NaN;
    baz=az; baz(cond2)=NaN;

    allax(:,i) = gax;
    allay(:,i) = gay;
    allaz(:,i) = gaz;
    
    % plot ax
    ah(3)=subplot(5,2,[5 6]);
    plot(d(i).t(1:pint3:end-2),gax(1:pint3:end),cols{i})
    hold on
    plot(d(i).t(1:pint3:end-2),bax(1:pint3:end),'color',[0.7 0.7 0.7])
    % cosmetics
    if i == length(d)
        grid on
        longticks([],3)
        xlim([d(i).t(1) d(i).t(end-2)])
        axcorr = corr(allax,'rows','complete');
        axavg = nanmean(nanmean(allax,2),1);
        try
            axout = rmoutliers(allax,'gesd');
        catch
            axout = rmoutliers(allax,'mean');
        end
        outpct = (length(allax)-length(axout))*100/length(allax);
        ylim([-max(axout,[],'all')-0.005*abs(max(axout,[],'all')) max(axout,[],'all')+0.005*abs(max(axout,[],'all'))])
        a=annotation('textbox',[0.77 0.61 0 0],'String',[sprintf('%05.2f%% Outliers',outpct)],'FitBoxToText','on');
        a.FontSize = 8;
        %boxtex('ur',ah(3),sprintf('%05.2f%% Outliers',outpct),10)
        b=annotation('textbox',[0.13 0.625 0 0],'String',[sprintf('%.2f, %.2f, %.2f,\n%.2f, %.2f, %.2f',axcorr(1,2),axcorr(1,3),axcorr(1,4),axcorr(2,3),axcorr(2,4),axcorr(3,4))],'FitBoxToText','on');
        b.FontSize = 8;
        text(d(i).t(10),0.8*max(axout,[],'all'),sprintf('mean = %f cm/s^2',axavg),'FontSize',8)
        ylabel('a_x [cm/s^2]')
        t3=title(sprintf('Ship Acceleration Components (Every %dth Point)',pint3));
        xticklabels([])
    end

    % plot ay
    ah(4)=subplot(5,2,[7 8]);
    plot(d(i).t(1:pint3:end-2),gay(1:pint3:end),cols{i})
    hold on
    plot(d(i).t(1:pint3:end-2),bay(1:pint3:end),'color',[0.7 0.7 0.7])
    % cosmetics
    if i == length(d)
        grid on
        longticks([],3)
        xlim([d(i).t(1) d(i).t(end-2)])
        aycorr = corr(allay,'rows','complete');
        ayavg = nanmean(nanmean(allay,2),1);
        try
            ayout = rmoutliers(allay,'gesd');
        catch
            ayout = rmoutliers(allay,'mean');
        end
        outpct = (length(allay)-length(ayout))*100/length(allay);
        ylim([-max(ayout,[],'all')-0.005*abs(max(ayout,[],'all')) max(ayout,[],'all')+0.005*abs(max(ayout,[],'all'))])
        a=annotation('textbox',[0.77 0.4375 0 0],'String',[sprintf('%05.2f%% Outliers',outpct)],'FitBoxToText','on');
        a.FontSize = 8;
        b=annotation('textbox',[0.13 0.45 0 0],'String',[sprintf('%.2f, %.2f, %.2f,\n%.2f, %.2f, %.2f',aycorr(1,2),aycorr(1,3),aycorr(1,4),aycorr(2,3),aycorr(2,4),aycorr(3,4))],'FitBoxToText','on');
        b.FontSize = 8;
        c=annotation('textbox',[0.4 0.45 0 0],'String',[sprintf('GPS 1 - red, GPS 2 - green,\nGPS 3 - blue, GPS 4 - black')],'FitBoxToText','on');
        c.FontSize = 8;
        text(d(i).t(10),0.8*max(ayout,[],'all'),sprintf('mean = %f cm/s^2',ayavg),'FontSize',8)
        ylabel('a_y [cm/s^2]')
        xticklabels([])        
    end

    % plot az
    ah(5)=subplot(5,2,[9 10]);
    plot(d(i).t(1:pint3:end-2),gaz(1:pint3:end),cols{i})
    hold on
    plot(d(i).t(1:pint3:end-2),baz(1:pint3:end),'color',[0.7 0.7 0.7])
    % cosmetics
    if i == length(d)
        grid on
        longticks([],3)
        xlim([d(i).t(1) d(i).t(end-2)])
        azcorr = corr(allaz,'rows','complete');
        azavg = nanmean(nanmean(allaz,2),1);
        try
            azout = rmoutliers(allaz,'gesd');
        catch
            azout = rmoutliers(allaz,'mean');
        end
        outpct = (length(allaz)-length(azout))*100/length(allaz);
        ylim([-max(azout,[],'all')-0.005*abs(max(azout,[],'all')) max(azout,[],'all')+0.005*abs(max(azout,[],'all'))])
        a=annotation('textbox',[0.77 0.265 0 0],'String',[sprintf('%05.2f%% Outliers',outpct)],'FitBoxToText','on');
        a.FontSize = 8;
        b=annotation('textbox',[0.13 0.2775 0 0],'String',[sprintf('%.2f, %.2f, %.2f,\n%.2f, %.2f, %.2f',azcorr(1,2),azcorr(1,3),azcorr(1,4),azcorr(2,3),azcorr(2,4),azcorr(3,4))],'FitBoxToText','on');
        b.FontSize = 8;
        c=annotation('textbox',[0.44 0.2775 0 0],'String',[sprintf('X12, X13, X14,\nX23, X24, X34')],'FitBoxToText','on');
        c.FontSize = 8;
        text(d(i).t(10),0.8*max(azout,[],'all'),sprintf('mean = %f cm/s^2',azavg),'FontSize',8)
        ylabel('a_z [cm/s^2]')
    end        
end

% plot distances between all GPS receivers with separate loop
% # of distances = all combos of # of GPS recs taken 2 at a time
pairs = nchoosek(1:length(d),2);
for i = 1:length(pairs)
    % compute Euclidean distance between GPS pair
    dist = sqrt((d(pairs(i,1)).xyz(:,1)-d(pairs(i,2)).xyz(:,1)).^2 + (d(pairs(i,1)).xyz(:,2)-d(pairs(i,2)).xyz(:,2)).^2 + (d(pairs(i,1)).xyz(:,3)-d(pairs(i,2)).xyz(:,3)).^2);
    % subtract off mean of data and offset each pair by constant
    normdist = dist - nanmean(dist) + i;
    % remove outliers
    dist(any(isnan(dist),2),:)=[];
    % fit line to data
    p = polyfit([1:length(dist)]',dist,1);
    % compute some statistics
    a = 1000*p(1); b = p(2); rmsd = rms(dist); stdd = std(1000*dist);
    x = (a/1000).*[1:length(dist)]' + b; e = x - dist; erms = 1000*rms(e);

    cond3=d(pairs(i,1)).pdop<pthresh & d(pairs(i,1)).pdop~=0 & d(pairs(i,1)).nsats(:,1)>nthresh ...
          & d(pairs(i,2)).pdop<pthresh & d(pairs(i,2)).pdop~=0 & d(pairs(i,2)).nsats(:,1)>nthresh;

    % only good data
    gd=normdist; gd(~cond3)=NaN;
    % only bad data
    bd=normdist; bd(cond3)=NaN;

    % actually plot GPS pair distance w.r.t time
    ah(2)=subplot(5,2,[2 4]);
    plot(d(1).t(1:pint2:end),gd(1:pint2:end))
    hold on
    plot(d(1).t(1:pint2:end),bd(1:pint2:end),'color',[0.7 0.7 0.7])
    text(d(1).t(10),i+0.4,sprintf('%f, %05.3f, %05.3f, %.0f, %.0f',a,b,rmsd,stdd,erms),'FontSize',9)
    ah(2).YTickLabel{i} = sprintf('%g-%g',pairs(i,1),pairs(i,2));
    % cosmetics
    if i == length(pairs)
        xlim([d(1).t(1) d(1).t(end)])
        xticklabels([])
        ylim([0.25 length(pairs)+0.75])
        ylabel('GPS Pair')
        t2=title(sprintf('Distances between GPS Receivers\n(Every %dth Point)',pint2));
        grid on
        longticks
        text(d(1).t(10),0.5,sprintf('a [mm/s], b [m], rms(x) [m], std [mm], rms(e) [mm]'),'FontSize',7.5)
    end
end

% finishing touches
tt=supertit(ah([1 2]),sprintf('Ship Data from %s to %s',datestr(d(1).t(1)),datestr(d(1).t(end))));
movev(tt,0.3)

%a = annotation('textbox',[0.465 0.085 0 0],'String',['camp'],'FitBoxToText','on');
%a.FontSize = 12;

figdisp(sprintf('gps2dis-%s',fname),[],'',2,[],'epstopdf')
