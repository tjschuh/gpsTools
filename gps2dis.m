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
% Last modified by fjsimons-at-princeton.edu, 05/17/2022

%NOT DONE

% need some serious edits to this:
% fix annotations/text boxes
% change input structure
% currently statistics (correlation coeff, ployfit, rms, std)
% are computed using all data including greyed out parts

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
    
    % plot heights of all 4 units all on 1 plot
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
        title(sprintf('Ship Height (Every %dth Point)',pint1))
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
        text(d(i).t(floor(0.05*length(d(i).t))),ah(1).YLim(1)+0.005*abs(ah(1).YLim(1)),...
             sprintf('%05.2f%% Outliers',outpct))
        text(d(i).t(floor(0.45*length(d(i).t))),ah(1).YLim(1)+0.005*abs(ah(1).YLim(1)),...
             sprintf('Nsats > %d & PDOP < %d',nthresh,pthresh))
        text(d(i).t(floor(0.6*length(d(i).t))),ah(1).YLim(2)-0.005*abs(ah(1).YLim(2)),...
             sprintf('v = %.2f knots',vavg))
    end

    %PUT DISTANCES CODE HERE
    
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
        b=annotation('textbox',[0.13 0.625 0 0],'String',[sprintf('%.2f, %.2f, %.2f,\n%.2f, %.2f, %.2f',axcorr(1,2),axcorr(1,3),axcorr(1,4),axcorr(2,3),axcorr(2,4),axcorr(3,4))],'FitBoxToText','on');
        b.FontSize = 8;
        text(d(i).t(10),0.8*max(axout,[],'all'),sprintf('mean = %f cm/s^2',axavg),'FontSize',8)
        ylabel('a_x [cm/s^2]')
        title(sprintf('Ship Acceleration Components (Every %dth Point)',pint3))
        xticklabels([])
    end

    % plot ay
    ah(4)=subplot(5,2,[7 8]);
    plot(d(i).t(1:pint3:end-2),gay(1:pint3:end),cols{i})
    hold on
    plot(d(i).t(1:pint3:end-2),bay(1:pint3:end),'color',[0.7 0.7 0.7])
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
keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot distances between all 4 GPS units
% compute all 6 distances between 4 GPS receivers
dist12 = sqrt((d1.xyz(:,1)-d2.xyz(:,1)).^2 + (d1.xyz(:,2)-d2.xyz(:,2)).^2 + (d1.xyz(:,3)-d2.xyz(:,3)).^2);
normdist12 = dist12 - nanmean(dist12) + 1;
dist12 = rmNaNrows(dist12); p = polyfit([1:length(dist12)]',dist12,1);
a12 = 1000*p(1); b12 = p(2); rms12 = rms(dist12); std12 = std(1000*dist12);
x12 = (a12/1000).*[1:length(dist12)]' + b12; e12 = x12 - dist12; erms12 = 1000*rms(e12);
dist13 = sqrt((d1.xyz(:,1)-d3.xyz(:,1)).^2 + (d1.xyz(:,2)-d3.xyz(:,2)).^2 + (d1.xyz(:,3)-d3.xyz(:,3)).^2);
normdist13 = dist13 - nanmean(dist13) + 2;
dist13 = rmNaNrows(dist13); p = polyfit([1:length(dist13)]',dist13,1);
a13 = 1000*p(1); b13 = p(2); rms13 = rms(dist13); std13 = std(1000*dist13);
x13 = (a13/1000).*[1:length(dist13)]' + b13; e13 = x13 - dist13; erms13 = 1000*rms(e13);
dist14 = sqrt((d1.xyz(:,1)-d4.xyz(:,1)).^2 + (d1.xyz(:,2)-d4.xyz(:,2)).^2 + (d1.xyz(:,3)-d4.xyz(:,3)).^2);
normdist14 = dist14 - nanmean(dist14) + 3;
dist14 = rmNaNrows(dist14); p = polyfit([1:length(dist14)]',dist14,1);
a14 = 1000*p(1); b14 = p(2); rms14 = rms(dist14); std14 = std(1000*dist14);
x14 = (a14/1000).*[1:length(dist14)]' + b14; e14 = x14 - dist14; erms14 = 1000*rms(e14);
dist23 = sqrt((d2.xyz(:,1)-d3.xyz(:,1)).^2 + (d2.xyz(:,2)-d3.xyz(:,2)).^2 + (d2.xyz(:,3)-d3.xyz(:,3)).^2);
normdist23 = dist23 - nanmean(dist23) + 4;
dist23 = rmNaNrows(dist23); p = polyfit([1:length(dist23)]',dist23,1);
a23 = 1000*p(1); b23 = p(2); rms23 = rms(dist23); std23 = std(1000*dist23);
x23 = (a23/1000).*[1:length(dist23)]' + b23; e23 = x23 - dist23; erms23 = 1000*rms(e23);
dist24 = sqrt((d2.xyz(:,1)-d4.xyz(:,1)).^2 + (d2.xyz(:,2)-d4.xyz(:,2)).^2 + (d2.xyz(:,3)-d4.xyz(:,3)).^2);
normdist24 = dist24 - nanmean(dist24) + 5;
dist24 = rmNaNrows(dist24); p = polyfit([1:length(dist24)]',dist24,1);
a24 = 1000*p(1); b24 = p(2); rms24 = rms(dist24); std24 = std(1000*dist24);
x24 = (a24/1000).*[1:length(dist24)]' + b24; e24 = x24 - dist24; erms24 = 1000*rms(e24);
dist34 = sqrt((d3.xyz(:,1)-d4.xyz(:,1)).^2 + (d3.xyz(:,2)-d4.xyz(:,2)).^2 + (d3.xyz(:,3)-d4.xyz(:,3)).^2);
normdist34 = dist34 - nanmean(dist34) + 6;
dist34 = rmNaNrows(dist34); p = polyfit([1:length(dist34)]',dist34,1);
a34 = 1000*p(1); b34 = p(2); rms34 = rms(dist34); std34 = std(1000*dist34);
x34 = (a34/1000).*[1:length(dist34)]' + b34; e34 = x34 - dist34; erms34 = 1000*rms(e34);

% find good (g) and bad (b) data
% [g b] = h
d12 = normdist12; d13 = normdist13; d14 = normdist14;
d23 = normdist23; d24 = normdist24; d34 = normdist34;
good12 = d12; bad12 = d12; good13 = d13; bad13 = d13; good14 = d14; bad14 = d14;
good23 = d23; bad23 = d23; good24 = d24; bad24 = d24; good34 = d34; bad34 = d34;
good12(p1>=pthresh | p1==0 | n1<=nthresh | p2>=pthresh | p2==0 | n2<=nthresh) = NaN;
bad12(p1<pthresh & n1>nthresh & p2<pthresh & n2>nthresh) = NaN;
good13(p1>=pthresh | p1==0 | n1<=nthresh | p3>=pthresh | p3==0 | n3<=nthresh) = NaN;
bad13(p1<pthresh & n1>nthresh & p3<pthresh & n3>nthresh) = NaN;
good14(p1>=pthresh | p1==0 | n1<=nthresh | p4>=pthresh | p4==0 | n4<=nthresh) = NaN;
bad14(p1<pthresh & n1>nthresh & p4<pthresh & n4>nthresh) = NaN;
good23(p2>=pthresh | p2==0 | n2<=nthresh | p3>=pthresh | p3==0 | n3<=nthresh) = NaN;
bad23(p2<pthresh & n2>nthresh & p3<pthresh & n3>nthresh) = NaN;
good24(p2>=pthresh | p2==0 | n2<=nthresh | p4>=pthresh | p4==0 | n4<=nthresh) = NaN;
bad24(p2<pthresh & n2>nthresh & p4<pthresh & n4>nthresh) = NaN;
good34(p3>=pthresh | p3==0 | n3<=nthresh | p4>=pthresh | p4==0 | n4<=nthresh) = NaN;
bad34(p3<pthresh & n3>nthresh & p4<pthresh & n4>nthresh) = NaN;

ah(2)=subplot(5,2,[2 4]);
plot(d1.t(1:pint2:end),good12(1:pint2:end))
hold on
plot(d1.t(1:pint2:end),good13(1:pint2:end))
plot(d1.t(1:pint2:end),good14(1:pint2:end))
plot(d1.t(1:pint2:end),good23(1:pint2:end))
plot(d1.t(1:pint2:end),good24(1:pint2:end))
plot(d1.t(1:pint2:end),good34(1:pint2:end))
xlim([d1.t(1) d1.t(end)])
xticklabels([])
ylim([0.25 6.75])
yticklabels({'1-2','1-3','1-4','2-3','2-4','3-4'})
ylabel('GPS Pair')
gpst=title(sprintf('Distances between GPS Receivers\n(Every %dth Point)',pint2));
grid on
longticks
% grey out bad data
plot(d1.t(1:pint2:end),bad12(1:pint2:end),'color',[0.7 0.7 0.7])
plot(d1.t(1:pint2:end),bad13(1:pint2:end),'color',[0.7 0.7 0.7])
plot(d1.t(1:pint2:end),bad14(1:pint2:end),'color',[0.7 0.7 0.7])
plot(d1.t(1:pint2:end),bad23(1:pint2:end),'color',[0.7 0.7 0.7])
plot(d1.t(1:pint2:end),bad24(1:pint2:end),'color',[0.7 0.7 0.7])
plot(d1.t(1:pint2:end),bad34(1:pint2:end),'color',[0.7 0.7 0.7])
text(d1.t(10),6.4,sprintf('%f, %05.3f, %05.3f, %.0f, %.0f',a34,b34,rms34,std34,erms34),'FontSize',9)
text(d1.t(10),5.4,sprintf('%f, %05.3f, %05.3f, %.0f, %.0f',a24,b24,rms24,std24,erms24),'FontSize',9)
text(d1.t(10),4.4,sprintf('%f, %05.3f, %05.3f, %.0f, %.0f',a23,b23,rms23,std23,erms23),'FontSize',9)
text(d1.t(10),3.4,sprintf('%f, %05.3f, %05.3f, %.0f, %.0f',a14,b14,rms14,std14,erms14),'FontSize',9)
text(d1.t(10),2.4,sprintf('%f, %05.3f, %05.3f, %.0f, %.0f',a13,b13,rms13,std13,erms13),'FontSize',9)
text(d1.t(10),1.4,sprintf('%f, %05.3f, %05.3f, %.0f, %.0f',a12,b12,rms12,std12,erms12),'FontSize',9)
text(d1.t(10),0.5,sprintf('a [mm/s], b [m], rms(x) [m], std [mm], rms(e) [mm]'),'FontSize',7.5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% finishing touches
tt=supertit(ah([1 2]),sprintf('Ship Data from %s to %s',datestr(d1.t(1)),datestr(d1.t(end))));
movev(tt,0.3)

a = annotation('textbox',[0.465 0.085 0 0],'String',['camp'],'FitBoxToText','on');
a.FontSize = 12;
keyboard
figdisp(sprintf('all4plt-%s',fname),[],'',2,[],'epstopdf')

close
