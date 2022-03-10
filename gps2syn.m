function varargout = gps2syn(d,tmax,xyz,v,expnum,exptype)
% GPS2SYN(d,tmax,xyz,v,expnum,exptype)
%
% This kinda works, plots st and hst, diff(st-hst) and bestfit line,
% (eventually DOP), ship trajectory with C-DOG location, and histogram of
% residuals for a given distance (x,y,z) from true location. Records if
% residuals for those perturbations are "normal" with either solid or empty
% circles on running x,y,z 2D plots. Does this for x number of trials (expnum)
%
% INPUT:
%
% d            single dataset with x,y,z,t ship locations
% tmax         last time in d
% xyz          beacon location [x y z] [m]
% v            sound speed [m/s]
% expnum       number of experiments to run
% exptype      type of number generation (random or system)
%
% EXAMPLE:
%
% load Unit1234-camp.mat
% gps2syn(d,tmax,[],[],[],[]);
%
% Originally written by tschuh-at-princeton.edu, 02/23/2022
% Last modified by tschuh-at-princeton.edu, 03/09/2022

% TO-DO
% play with velocity
% play with system method

% C-DOG location [x,y,z] [m]
defval('xyz',[1.977967 -5.073198 3.3101016]*1e6)

% constant sound speed profile for now [m/s]
defval('v',1500)

% call gps2fwd to generate perfect st data
[st,xyz0,v0] = gps2fwd(d,tmax,xyz,v);

% find the rows that dont contain NaNs
rowsx=find(~isnan(d.x));
rowsy=find(~isnan(d.y));
rowsz=find(~isnan(d.z));
rows=unique([rowsx;rowsy;rowsz]);

% number of trials to run
defval('expnum',1)

% experiment type
defval('exptype','system')

% marker size
sz=20;
% counter for keeping track of points we care about
counter = 1;
% unit multipliers across time and space
tmulti{1,1} = 1e6; tmulti{1,2} = '\mus';
xmulti{1,1} = 1e2; xmulti{1,2} = 'cm';
clf

% allocate array for number of experiments
xyzn = zeros(expnum,3);

for i=1:expnum
    if exptype == 'random'
        % generate random xyz perturbations between -10 xmulti{1,2} and 10 xmulti{1,2}
        xyzn(i,:) = [(randi(201)-101) (randi(201)-101) (randi(201)-101)]./(xmulti{1,1}*10);
    elseif exptype == 'system'
        % need to figure out how to generate multiple grid-based numbers
        xyzn(i,:) = [-3.3 -8.9 4.8]./xmulti{1,1};
    end
    % add xyzn to xyz0 to get perturbed C-DOG location
    xyzg = xyz0 + xyzn(i,:);
    % generate a random velocity perturbation
    vn = 0;
    % add vn to v0 to get perturbed sound speed
    vg = v0 + vn;
    % forward model: creating perturbed data
    hsr = sqrt((d.x-xyzg(1)).^2 + (d.y-xyzg(2)).^2 + (d.z-xyzg(3)).^2);
    hst = hsr./vg;
    % error between slant times and predicted slant times
    differ = st-hst;
    % calculate relative time difference
    rel = differ./st;
    % calculate rms(st-hst)
    rmse = norm(differ(rows));

    figure(1)
    clf
    lwidth = 1.5;

    % plot st and hst
    ah(1)=subplot(2,4,[1 2]);
    p(1)=plot(d.t,st,'-','LineWidth',lwidth);
    hold on
    p(2)=plot(d.t,hst,'--','LineWidth',lwidth);
    hold off
    datetick('x','HH')
    xticklabels([])
    cosmot(d.t)
    ylabel('slant range time [s]')
    title('st and hst')
    legend({'st','hst'})

    % plot absolute time differences
    ah(2)=subplot(2,4,[5 6]);
    times = seconds(d.t(rows)-d.t(rows(1)));
    pf = polyfit(times,differ(rows),1);
    bfline = polyval(pf,times);
    p(3)=plot(d.t(rows),bfline*tmulti{1,1},'LineWidth',lwidth);
    hold on
    p(4)=plot(d.t,differ*tmulti{1,1},'LineWidth',lwidth);
    %text(d.t(1000),0.9*ah(2).YLim(2),sprintf('avg = %3.0f',mean(differ(rows)*tmulti{1,1})));
    hold off
    datetick('x','HH')
    cosmot(d.t)
    title('Difference between st and hst')
    ylabel(sprintf('slant range time [%s]',tmulti{1,2}))
    xlabel('time [h]')

    % plot histogram of relative time differences
    ah(3)=subplot(2,4,[3 4]);
    thresh1=1000;
    nstd1 = 1;
    [b,gof1]=cosmoh(rel(rows)*tmulti{1,1},ah(3),thresh1,tmulti{1,2},nstd1);
    b.FaceColor = [0.400 0.6667 0.8431];
    title('Relative Time Differences')
    
    % plot histogram of absolute time differences
    ah(4)=subplot(2,4,[7 8]);
    thresh2=1000;
    nstd2 = 2;
    [c,gof2]=cosmoh(differ(rows)*tmulti{1,1},ah(4),thresh2,tmulti{1,2},nstd2);
    c.FaceColor = [0.8500 0.3250 0.0980];
    title('Absolute Time Differences')
    
    %dop
    %for i=1:length(d.x)
        %A = [(d.x(i,1)-xyzg(1)) (d.y(i,1)-xyzg(2)) (d.z(i,1)-xyzg(3))]./hsr(i);
        %A = [A -ones([length(A) 1])];
        %A = [A -1];
        %Q = inv(A'*A);
        %GDOP = sqrt(trace(Q));
    %end

    tt=supertit(ah([1 3]),sprintf('Distance from Truth = [%g %s %g %s %g %s] = |%3.3f %s|',xmulti{1,1}*xyzn(i,1),xmulti{1,2},xmulti{1,1}*xyzn(i,2),xmulti{1,2},xmulti{1,1}*xyzn(i,3),xmulti{1,2},xmulti{1,1}*sqrt(xyzn(i,1)^2+xyzn(i,2)^2+xyzn(i,3)^2),xmulti{1,2}));
    movev(tt,0.4);
    
    figdisp(sprintf('experiment-%d',i),[],[],2,[],'epstopdf')

    g=figure(2);
    g.Visible = 'off';
    ahh(1)=subplot(2,2,1);
    if gof2 > thresh2
        scatter(xmulti{1,1}*xyzn(i,1),xmulti{1,1}*xyzn(i,2),sz)
    else
        scatter(xmulti{1,1}*xyzn(i,1),xmulti{1,1}*xyzn(i,2),sz,'filled')
    end
    hold on
    cosmoxyz()
    xlabel(sprintf('perturbations in x [%s]',xmulti{1,2}))
    ylabel(sprintf('perturbations in y [%s]',xmulti{1,2}))

    ahh(2)=subplot(2,2,2);
    if gof2 > thresh2
        scatter(xmulti{1,1}*xyzn(i,1),xmulti{1,1}*xyzn(i,3),sz)
    else
        scatter(xmulti{1,1}*xyzn(i,1),xmulti{1,1}*xyzn(i,3),sz,'filled')
    end
    hold on
    cosmoxyz()
    xlabel(sprintf('perturbations in x [%s]',xmulti{1,2}))
    ylabel(sprintf('perturbations in z [%s]',xmulti{1,2}))

    ahh(3)=subplot(2,2,3);
    if gof2 > thresh2
        scatter(xmulti{1,1}*xyzn(i,2),xmulti{1,1}*xyzn(i,3),sz)
    else
        scatter(xmulti{1,1}*xyzn(i,2),xmulti{1,1}*xyzn(i,3),sz,'filled')
        % save perturbations to new matrix for later use
        ellip(counter,:) = xyzn(i,:);
        counter = counter + 1;
    end
    hold on
    cosmoxyz()
    xlabel(sprintf('perturbations in y [%s]',xmulti{1,2}))
    ylabel(sprintf('perturbations in z [%s]',xmulti{1,2}))
    g.Visible = 'off';

    ahh(4)=subplot(2,2,4);
    %trajectory of ship w/ C-DOG
    if i == expnum
        g.Visible = 'off';
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
        xl(2)=xlabel('easting [km]');
        yl(2)=ylabel('northing [km]');
        openup(ahh(4),5,10);
        openup(ahh(4),6,10);
        scatter((dogx-refx)/sclx,(dogy-refy)/scly,10,...
                'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r')
        [lat,lon,~] = ecef2geodetic(wgs84,xyzg(1),xyzg(2),xyzg(3));
        [fakex,fakey,~] = deg2utm(lat,mod(lon,360));
        %scatter((fakex-refx)/sclx,(fakey-refy)/scly,10,...
        %        'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','r')
        hold off
        g.Visible = 'off';
    end
    g.Visible = 'off';    

    ttt=supertit(ahh([1 2]),sprintf('Results of Experiment'));
    movev(ttt,0.3)
end

g.Visible = 'on';
figdisp(sprintf('trials'),[],[],2,[],'epstopdf')
%close

% if ellip exists, plot it
if exist('ellip','var') == 1
    figure(3)
    % maybe plot with contours?
    scatter3(ellip(:,1),ellip(:,2),ellip(:,3),sz,'filled')
    
    figdisp(sprintf('ellipsoid'),[],[],2,[],'epstopdf')
end

% optional output
varns={ellip};
varargout=varns(1:nargout);

function cosmot(t)
xlim([t(1) t(end)])
grid on
longticks([],2)

function cosmoxyz()
grid on
longticks
xlim([-10 10])
ylim([-10 10])

function [b,gof]=cosmoh(data,ax,thresh,multi,nstd)
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
longticks([],2)
stdd=std(data);
%nstd=3;
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
grid on
xlabel(sprintf('residuals [%s]',multi))
t=text(-0.95*nstd*stdd,0.75*ax.YLim(2),...
       sprintf('N = %3.0f\nstd = %3.0f\nmed = %3.0f\navg = %3.0f\ngof = %3.0f',...
                    length(data),stdd,median(data),mean(data),gof));
hold on
pd = fitdist(data,'Normal');
xvals = b.XData; 
yvals = pdf(pd,xvals);
area = sum(b.YData)*diff(b.XData(1:2));
lain = plot(xvals,yvals*area,'r','LineWidth',2);
if gof > thresh
    lain.LineStyle = '--';
end
hold off