function varargout = gps2syn(st,d,vg,xyzg,expnum)
% GPS2SYN(st,d,vg,xyzg,expnum)
%
%
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
defval('expnum',1000)

xyzmat = zeros(expnum,4);

figure;
sz=20;

for i=1:50
    % generate a random decimal value between -10 cm and 10 cm
    xyzn = [(randi(201)-101)./1000 (randi(201)-101)./1000 (randi(201)-101)./1000];
    sol = xyz0 + xyzn;
    hsr = sqrt((d.xyz(:,1)-sol(1)).^2 + (d.xyz(:,2)-sol(2)).^2 + (d.xyz(:,3)-sol(3)).^2);
    hst = hsr./v0;
    %rmse = norm(st - hst);
    res = st-hst;

    ah(1)=subplot(3,3,[1 2]);
    p(1)=plot(d.t,st,'LineWidth',1.5);
    hold on
    p(2)=plot(d.t,hst,'LineWidth',1.5);
    hold off
    cosmot(d.t)
    ylabel('slant range time [s]')

    ah(2)=subplot(3,3,[4 5]);
    p(3)=plot(d.t,res*1e6,'LineWidth',1.5);
    cosmot(d.t)
    ylabel('slant range time [\mus]')

    ah(3)=subplot(3,3,[7 8]);
    nbins=round((max(res*1e6)-min(res*1e6))/(2*iqr(res*1e6)*(length(res*1e6))^(-1/3)));
    % calculate goodness of fit (gof) compared to a normal distribution
    %[~,~,stats]=chi2gof(res,'NBins',nbins);
    % divide chi squared by degrees of freedom to reduce to 1 DoF
    % with 1 Dof, chi squared <= 4 signifies ~90% chance data are normal
    % will make red curve dotted if gof > 4
    %gof=stats.chi2stat/stats.df;
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
    ylabel('counts')

    ah(4)=subplot(3,3,3);
    scatter(100*xyzn(1),100*xyzn(2),sz,'filled')
    hold on
    cosmoxyz()
    xlabel('perturbations in x [cm]')
    ylabel('perturbations in y [cm]')

    ah(5)=subplot(3,3,6);
    scatter(100*xyzn(1),100*xyzn(3),sz,'filled')
    hold on
    cosmoxyz()
    xlabel('perturbations in x [cm]')
    ylabel('perturbations in z [cm]')

    ah(6)=subplot(3,3,9);
    scatter(100*xyzn(2),100*xyzn(3),sz,'filled')
    hold on
    cosmoxyz()
    xlabel('perturbations in y [cm]')
    ylabel('perturbations in z [cm]')

    %delete(ah(1:3))
end

%keyboard

function cosmot(t)
xlim([t(1) t(end)])
grid on
longticks([],2)

function cosmoxyz()
grid on
longticks
xlim([-10 10])
ylim([-10 10])