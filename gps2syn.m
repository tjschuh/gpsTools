function varargout = gps2syn(st,d,vg,xyzg,expnum)
% GPS2SYN(st,d,vg,xyzg,expnum)
%
% This kinda works, plots st and hst, diff(st-hst) and bestfit line
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
% [st,xyz,v] = gps2fwd(d,tmax,[2e6 -4.5e6 3e6],1500);
% gps2syn(st,d,v,xyz,expnum);
%
% Originally written by tschuh-at-princeton.edu, 02/23/2022

% need to get rid of rows with NaNs
% need to combine d.t and d.xyz into 1 matrix to remove all NaN rows
%d.t(any(isnan(d.t),2),:)=[];
d.xyz(any(isnan(d.xyz),2),:)=NaN;
st(any(isnan(st),2),:)=NaN;

% constant sound speed profile for now [m/s]
defval('vg',1500)

% guess the correct location
defval('xyzg',[2e6 -4.5e6 3e6])

v0 = vg;
xyz0 = xyzg;

% number of trials to run
defval('expnum',100)

sz=20;

for i=1:expnum
    % generate a random decimal value between -10 cm and 10 cm
    xyzn = [(randi(201)-101)./1000 (randi(201)-101)./1000 (randi(201)-101)./1000];
    sol = xyz0 + xyzn;
    hsr = sqrt((d.xyz(:,1)-sol(1)).^2 + (d.xyz(:,2)-sol(2)).^2 + (d.xyz(:,3)-sol(3)).^2);
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
    %trajectory of ship

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

    ah(4)=subplot(3,4,[7 8]);
    nbins=round((max(res*1e6)-min(res*1e6))/(2*iqr(res*1e6)*(length(res*1e6))^(-1/3)));
    % calculate goodness of fit (gof) compared to a normal distribution
    [~,~,stats]=chi2gof(res,'NBins',nbins);
    % divide chi squared by degrees of freedom to reduce to 1 DoF
    % with 1 Dof, chi squared <= 4 signifies ~90% chance data are normal
    % will make red curve dotted if gof > 4
    gof=stats.chi2stat/stats.df;
    % Calculate histogram
    [yvals,edges]=histcounts(res*1e6,nbins);
    % Calculate bin centers 
    barc=0.5*(edges(1:end-1)+edges(2:end));
    % Plot the histogram
    b=bar(barc,yvals,'BarWidth',1);
    b.FaceColor = [0.400 0.6667 0.8431];
    longticks
    grid on
    xlabel('residuals [\mus]')
    hold on
    pd = fitdist(res*1e6,'Normal');
    xvals = b.XData; 
    yvals = pdf(pd,xvals);
    area = sum(b.YData)*diff(b.XData(1:2));
    lain = plot(xvals,yvals*area,'r','LineWidth',2);
    thresh = 100;
    if gof > thresh
        lain.LineStyle = '--';
    end
    hold off
    
    ah(5)=subplot(3,4,[9 10]);
    %pdop
    %A = [(d.xyz(:,1)-sol(1)) (d.xyz(:,2)-sol(2)) (d.xyz(:,3)-sol(3))]./hsr(:);
    %A = [A -ones([length(A) 1])];
    %keyboard

    tt=supertit(ah([1 2]),sprintf('Distance from Truth = [%g cm %g cm %g cm]',100*xyzn(1),100*xyzn(2),100*xyzn(3)));

    figdisp(sprintf('synthetic-%d',i),[],[],2,[],'epstopdf')
    close(f)

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

figdisp(sprintf('synthetic'),[],[],2,[],'epstopdf')
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