function varargout = gps2syn(st,d,vg,xyzg,expnum)
% GPS2SYN(st,d,vg,xyzg,expnum)
%
% This kinda works, plots st and hst, diff(st-hst) and bestfit line,
% (eventually DOP), ship trajectory with C-DOG location, and histogram of
% residuals for a given distance (x,y,z) from true location. Records if
% residuals for those perturbations are "normal" with either solid or empty
% circles on running x,y,z 2D plots. Does this for x number of trials (expnum)
%
% INPUT:
%
% st            slant times/arrival times of the signals from the
%               "source" to the receiver whose position we want to find
%               here. These are from the DOGS, or, synthetic, from GPS2RNG
% d
% vg            sound speed guess [m/s]
% xyzg          initial beacon location guess [x y z] [m]
%
% EXAMPLE:
%
% load Unit1234-camp.mat
% [st,xyz,v] = gps2fwd(d,tmax,[1.977967 -5.073198 3.3101016]*1e6,1500);
% gps2syn(st,d,v,xyz,expnum);
%
% Originally written by tschuh-at-princeton.edu, 02/23/2022
% Last modified by tschuh-at-princeton.edu, 03/01/2022

% need to get rid of rows with NaNs
% need to combine d.t and d.xyz into 1 matrix to remove all NaN rows
%d.t(any(isnan(d.t),2),:)=[];
% LEFT OFF HERE
keyboard
rowsx=find(isnan(d.x));
rowsy=find(isnan(d.y));
rowsz=find(isnan(d.z));
rows=unique([rowsx;rowsy;rowsz]);
d.t(rows,:)=[];
d.x(rows,:)=[];
d.y(rows,:)=[];
d.z(rows,:)=[];
st(rows,:)=[];
keyboard
% constant sound speed profile for now [m/s]
defval('vg',1500)

% guess the correct location [m]
defval('xyzg',[1.977967 -5.073198 3.3101016]*1e6)

v0 = vg;
xyz0 = xyzg;

% number of trials to run
defval('expnum',100)

sz=20;

for i=1:expnum
    % generate a random decimal value between -10 cm and 10 cm
    xyzn = [(randi(201)-101)./1000 (randi(201)-101)./1000 (randi(201)-101)./1000];
    sol = xyz0 + xyzn;
    hsr = sqrt((d.x-sol(1)).^2 + (d.y-sol(2)).^2 + (d.z-sol(3)).^2);
    hst = hsr./v0;
    %rmse = norm(st - hst);
    differ = st-hst;

    f=figure;
    f.Position = [250 500 900 500];

    ah(1)=subplot(3,4,[1 2]);
    p(1)=plot(d.t,st,'LineWidth',1.5);
    hold on
    p(2)=plot(d.t,hst,'LineWidth',1.5);
    hold off
    cosmot(d.t)
    ylabel('slant range time [s]')

    ah(2)=subplot(3,4,[3 4]);
    %trajectory of ship w/ C-DOG
    skp=1000;
    zex=d.utme(1:skp:end);
    zwi=d.utmn(1:skp:end);
    refx=min(d.utme);
    refy=min(d.utmn);
    sclx=1000;
    scly=1000;
    plot((d.utme-refx)/sclx,(d.utmn-refy)/scly,'LineWidth',1,'Color','k');
    hold on
    scatter((zex-refx)/sclx,(zwi-refy)/scly,5,...
                   'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')
    wgs84 = wgs84Ellipsoid('meter');
    [lat,lon,~] = ecef2geodetic(wgs84,xyz0(1),xyz0(2),xyz0(3));
    [dogx,dogy,~] = deg2utm(lat,mod(lon,360));
    box on
    grid on
    longticks([],2)
    xl(2)=xlabel('easting (km)');
    yl(2)=ylabel('northing (km)');
    openup(ah(2),5,10);
    openup(ah(2),6,10);
    scatter((dogx-refx)/sclx,(dogy-refy)/scly,10,...
               'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r')
    hold off

    ah(3)=subplot(3,4,[5 6]);
    %differ(any(isnan(differ),2),:) = [];
    time = d.t;
    ind = isnan(differ);
    differ(ind) = [];
    time(ind) = [];
    pf = polyfit([1:length(differ)]',differ,1);
    a = pf(1)*1e6; b = pf(2);
    bfline = (a/1e6).*[1:length(differ)]' + b;
    res = bfline - differ;
    p(3)=plot(time,bfline*1e6,'LineWidth',1.5);
    hold on
    p(4)=plot(time,differ*1e6,'LineWidth',1.5);
    hold off
    cosmot(time)
    ylabel('slant range time [\mus]')

    ah(4)=subplot(3,4,[7 8 11 12]);
    data = res*1e6;
    nbins=round((max(data)-min(data))/(2*iqr(data)*(length(data))^(-1/3)));
    % calculate goodness of fit (gof) compared to a normal distribution
    [~,~,stats]=chi2gof(data,'NBins',nbins);
    % divide chi squared by degrees of freedom to reduce to 1 DoF
    % with 1 Dof, chi squared <= 4 signifies ~90% chance data are normal
    % will make red curve dotted if gof > 100
    gof=stats.chi2stat/stats.df;
    % Calculate histogram
    [yvals,edges]=histcounts(data,nbins);
    % Calculate bin centers 
    barc=0.5*(edges(1:end-1)+edges(2:end));
    % Plot the histogram
    b=bar(barc,yvals,'BarWidth',1);
    b.FaceColor = [0.400 0.6667 0.8431];
    longticks([],2)
    stdd=std(data);
    nstd=3;
    xel=round([-nstd nstd]*stdd,2);
    xlim(xel)
    ax.XTick=round([-nstd:nstd]*stdd,2);
    % Trick to do it properly
    bla=round([-nstd:nstd]*stdd,0);
    for in=1:length(bla)
        blace{in}=sprintf('%.0f',bla(in));
    end
    ax.XTickLabel=blace;
    yel=[0 max(b.YData)+0.1*max(b.YData)];
    ylim(yel)
    ax.YTick=unique(0:50:max(yel));
    if length(ax.YTick)>6
        ax.YTick=unique(0:100:max(yel));
    end
    grid on
    xlabel('residuals [\mus]')
    t=text(-0.85*nstd*stdd,0.75*ah(4).YLim(2),...
           sprintf('N = %3.0f\nstd = %3.0f\nmed = %3.0f\navg = %3.0f\ngof = %3.0f',...
                        length(data),stdd,median(data),mean(data),gof));
    hold on
    pd = fitdist(data,'Normal');
    xvals = b.XData; 
    yvals = pdf(pd,xvals);
    area = sum(b.YData)*diff(b.XData(1:2));
    lain = plot(xvals,yvals*area,'r','LineWidth',2);
    thresh = 1000;
    if gof > thresh
        lain.LineStyle = '--';
    end
    hold off
    
    ah(5)=subplot(3,4,[9 10]);
    %dop
    %for i=1:length(d.x)
        %A = [(d.x(i,1)-sol(1)) (d.y(i,1)-sol(2)) (d.z(i,1)-sol(3))]./hsr(i);
        %A = [A -ones([length(A) 1])];
        %A = [A -1];
        %Q = inv(A'*A);
        %GDOP = sqrt(trace(Q));
    %end

    tt=supertit(ah([1 2]),sprintf('Distance from Truth = [%g cm %g cm %g cm]',100*xyzn(1),100*xyzn(2),100*xyzn(3)));

    %figdisp(sprintf('synthetic-%d',i),[],[],2,[],'epstopdf')
    %close(f)

    if i == 1
        g=figure;
        g.Position = [250 500 300 600];
    end
    
    ahh(4)=subplot(3,1,1);
    if gof > thresh
        scatter(100*xyzn(1),100*xyzn(2),sz)
    else
        scatter(100*xyzn(1),100*xyzn(2),sz,'filled')
    end
    hold on
    cosmoxyz()
    xlabel('perturbations in x [cm]')
    ylabel('perturbations in y [cm]')

    ahh(5)=subplot(3,1,2);
    if gof > thresh
        scatter(100*xyzn(1),100*xyzn(3),sz)
    else
        scatter(100*xyzn(1),100*xyzn(3),sz,'filled')
    end
    hold on
    cosmoxyz()
    xlabel('perturbations in x [cm]')
    ylabel('perturbations in z [cm]')

    ahh(6)=subplot(3,1,3);
    if gof > thresh
        scatter(100*xyzn(2),100*xyzn(3),sz)
    else
        scatter(100*xyzn(2),100*xyzn(3),sz,'filled')
    end
    hold on
    cosmoxyz()
    xlabel('perturbations in y [cm]')
    ylabel('perturbations in z [cm]')

end

%figdisp(sprintf('synthetic'),[],[],2,[],'epstopdf')
%close

function cosmot(t)
xlim([t(1) t(end)])
grid on
longticks([],2)

function cosmoxyz()
grid on
longticks
xlim([-10 10])
ylim([-10 10])