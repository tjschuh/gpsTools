function varargout = gps2syn(d,tmax,xyz,xyzn,v,vn)
% xyzpdw = GPS2SYN(d,tmax,xyz,xyzn,v,vn)
%
% This kinda works, plots st and hst, diff(st-hst) and bestfit line,
% (eventually DOP), ship trajectory with C-DOG location, and histogram of
% residuals for a given distance (x,y,z) from true location. Records if
% residuals for those perturbations are "normal" with either solid or empty
% circles on running x,y,z 2D plots.
%
% INPUT:
%
% d            single dataset with x,y,z,t ship locations
% tmax         first and last time in d [2x1 datetime array]
% xyz          beacon location [x y z] [m]
% xyzn         perturbations to beacon location [x y z] [m]
% v            sound speed [m/s]
% vn           perturbation in sound speed [m/s]
%
% OUTPUT:
%
% xyzdwp       array containing xyzn, Durbin-Watson test result and p-value, 
%              and standard deviation of absolute time differences
% 
% EXAMPLE:
%
% load Unit1234-camp.mat
% xyzn=mesh(6,2);
% for i=1:length(xyzn)
% xyzdwp(i,:)=gps2syn(d,tmax,[],xyzn(i,:),[],[]);
% end
%
% Originally written by tschuh-at-princeton.edu, 02/23/2022
% Last modified by tschuh-at-princeton.edu, 04/11/2022

% C-DOG location [x,y,z] [m]
[x,y,z] = gps2dep([1.977967 -5.073198 3.3101016]*1e6,5225);
defval('xyz',[x y z])

% constant sound speed profile for now [m/s]
defval('v',1500)

% call gps2fwd to generate perfect st data
[st,xyz0,v0] = gps2fwd(d,tmax,xyz,v);

% find the rows that dont contain NaNs
rowsx=find(~isnan(d.x));
rowsy=find(~isnan(d.y));
rowsz=find(~isnan(d.z));
rows=unique([rowsx;rowsy;rowsz]);

defval('xyzn',[10 10 10])
defval('vn',0)

% unit multipliers across time and space
tmulti{1,1} = 1e6; tmulti{1,2} = '\mus';
xmulti{1,1} = 1e-3; xmulti{1,2} = 'mm';

% add xyzn to xyz0 to get perturbed C-DOG location
xyzg = xyz0 + xyzn*xmulti{1,1};
% add vn to v0 to get perturbed sound speed
vg = v0 + vn;

% add uncertainty to ship locations (+/- 0.02 m)
% draw a number from a normal distribution from -2 to 2
d.x0 = d.x; d.y0 = d.y; d.z0 = d.z;
for i = 1:length(d.x0)
    d.x(i) = d.x0(i) + (normrnd(0,2/3)*1e-2);
    d.y(i) = d.y0(i) + (normrnd(0,2/3)*1e-2);
    d.z(i) = d.z0(i) + (normrnd(0,2/3)*1e-2);
end

% forward model: creating perturbed data
%hsr = sqrt((d.x0-xyzg(1)).^2 + (d.y0-xyzg(2)).^2 + (d.z0-xyzg(3)).^2);
hsr = sqrt((d.x-xyzg(1)).^2 + (d.y-xyzg(2)).^2 + (d.z-xyzg(3)).^2);
hst = hsr./vg;
% error between slant times and predicted slant times
differ = st-hst;
% calculate relative time difference
rel = differ./st;
% calculate rms(st-hst)
rmse = norm(differ(rows));

% Durbin-Watson test to see if residuals (differ) are uncorrelated
% 0 <= dw <= 4, best result is dw = 2
% 1.5 <= dw <= 2.5 --> residuals are uncorrelated
% also want a p-value > 0.05
dwdata = differ;
dwdata(any(isnan(dwdata),2),:)=[];
% use dwtest with design matrix full of ones
[pval,dw] = dwtest(dwdata,ones(length(dwdata),1));
% set threshold
pthresh = 0.05;

% plotting
figure(1)
clf
lwidth = 1.5;

% plot st and hst
ah(1)=subplot(2,4,[1 2]);
p(1)=plot(d.t,hst,'-','color',[0.8500 0.3250 0.0980],'LineWidth',lwidth);
hold on
p(2)=plot(d.t,st,'b--','LineWidth',lwidth);
hold off
datetick('x','HH')
xticklabels([])
cosmot(d.t)
ylabel('slant range time [s]')
title('hst and st')
legend({'hst','st'})

% plot absolute time differences
ah(2)=subplot(2,4,[5 6]);
times = seconds(d.t(rows)-d.t(rows(1)));
pf = polyfit(times,differ(rows),1);
bfline = polyval(pf,times);
p(3)=plot(d.t,differ*tmulti{1,1},'color',[0.8500 0.3250 0.0980],'LineWidth',lwidth);
hold on
p(4)=plot(d.t(rows),bfline*tmulti{1,1},'b','LineWidth',lwidth);
hold off
datetick('x','HH')
cosmot(d.t)
ylim(1.1*[-max(abs(differ*tmulti{1,1})) max(abs(differ*tmulti{1,1}))])
title('Difference between st and hst')
ylabel(sprintf('slant range time [%s]',tmulti{1,2}))
xlabel('time [h]')
text(ah(2).XLim(1)+0.01*abs(ah(2).XLim(2)-ah(2).XLim(1)),0.9*ah(2).YLim(1),...
     sprintf('dw = %.3f, p = %.3g',dw,pval));
if pval >= pthresh
    text(ah(2).XLim(2)-0.25*abs(ah(2).XLim(2)-ah(2).XLim(1)),0.9*ah(2).YLim(1),...
         sprintf('ACCEPTED'))
else
    text(ah(2).XLim(2)-0.25*abs(ah(2).XLim(2)-ah(2).XLim(1)),0.9*ah(2).YLim(1),...
         sprintf('REJECTED'))
end

% plot histogram of relative time differences
ah(3)=subplot(2,4,[3 4]);
try
    thresh1=500;
    [b,lain1,gof1,stdr]=cosmoh(rel(rows)*tmulti{1,1},ah(3),thresh1,tmulti{1,2});
    b.FaceColor = [0.400 0.6667 0.8431];
    lain1.Color = 'blue';
    title('Relative Time Differences')
catch
end

% plot histogram of absolute time differences
ah(4)=subplot(2,4,[7 8]);
try
    thresh2=500;
    [c,lain2,gof2,stda]=cosmoh(differ(rows)*tmulti{1,1},ah(4),thresh2,tmulti{1,2});
    c.FaceColor = [0.8500 0.3250 0.0980];
    lain2.Color = 'red';
    title('Absolute Time Differences')
catch
end

%DOP
%for i=1:length(d.x)
%A = [(d.x(i,1)-xyzg(1)) (d.y(i,1)-xyzg(2)) (d.z(i,1)-xyzg(3))]./hsr(i);
%A = [A -ones([length(A) 1])];
%Q = inv(A'*A);
%GDOP = sqrt(trace(Q));
%end

tt=supertit(ah([1 3]),sprintf('Distance from Truth = [%g %s %g %s %g %s] = |%3.3f %s|, True Sound Speed = %g m/s\nSound Speed Error = %g m/s, GPS Perturbations = +/-[2 cm 2 cm 2 cm]',xyzn(1),xmulti{1,2},xyzn(2),xmulti{1,2},xyzn(3),xmulti{1,2},sqrt(xyzn(1)^2+xyzn(2)^2+xyzn(3)^2),xmulti{1,2},v0,vg-v0));
tt.FontSize = 12;
movev(tt,0.325);

figdisp(sprintf('experiment_%g_%g_%g',xyzn(1),xyzn(2),xyzn(3)),[],[],2,[],'epstopdf')

    % %trajectory of ship w/ C-DOG
    %     skp=1000;
    %     zex=d.utme(1:skp:end);
    %     zwi=d.utmn(1:skp:end);
    %     refx=min(d.utme);
    %     refy=min(d.utmn);
    %     sclx=1000;
    %     scly=1000;
    %     plot((d.utme-refx)/sclx,(d.utmn-refy)/scly,'LineWidth',1,'Color','k');
    %     hold on
    %     scatter((zex-refx)/sclx,(zwi-refy)/scly,5,...
    %             'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')
    %     wgs84 = wgs84Ellipsoid('meter');
    %     [lat,lon,~] = ecef2geodetic(wgs84,xyz0(1),xyz0(2),xyz0(3));
    %     [dogx,dogy,~] = deg2utm(lat,mod(lon,360));
    %     box on
    %     grid on
    %     longticks([],2)
    %     xl(2)=xlabel('easting [km]');
    %     yl(2)=ylabel('northing [km]');
    %     openup(ahh(4),5,10);
    %     openup(ahh(4),6,10);
    %     scatter((dogx-refx)/sclx,(dogy-refy)/scly,10,...
    %             'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r')
    %     [lat,lon,~] = ecef2geodetic(wgs84,xyzg(1),xyzg(2),xyzg(3));
    %     [fakex,fakey,~] = deg2utm(lat,mod(lon,360));
    %     %scatter((fakex-refx)/sclx,(fakey-refy)/scly,10,...
    %     %        'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','r')
    %     hold off
    %     g.Visible = 'off';

    % ttt=supertit(ahh([1 2]),sprintf('Results of Experiment'));
    % movev(ttt,0.3)


    %g.Visible = 'on';
    %figdisp(sprintf('trials'),[],[],2,[],'epstopdf')
    %close

% if ellip exists, plot it
% if exist('ellip','var') == 1
%     figure(3)
%     % maybe plot with contours?
%     scatter3(ellip(:,1)*xmulti{1,1},ellip(:,2)*xmulti{1,1},ellip(:,3)*xmulti{1,1},sz,'filled')
%     xlim([-10 10])
%     ylim([-10 10])
%     zlim([-10 10])
%     xlabel(sprintf('perturbations in x [%s]',xmulti{1,2}))
%     ylabel(sprintf('perturbations in y [%s]',xmulti{1,2}))
%     zlabel(sprintf('perturbations in z [%s]',xmulti{1,2}))
    
%     figdisp(sprintf('ellipsoid'),[],[],2,[],'epstopdf')
% end

% save some results
xyzdwp = [xyzn(1) xyzn(2) xyzn(3) dw pval stda];

% optional output
varns={xyzdwp};
varargout=varns(1:nargout);
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

function [b,lain,gof,stdd]=cosmoh(data,ax,thresh,multi)
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
nstd=3;
%xel=round(mean(data)+[-nstd nstd]*stdd,2);
%xlim(xel)
%ax.XTick=round([-nstd:nstd]*stdd,2);
% Trick to do it properly
%bla=round([-nstd:nstd]*stdd,0);
%for in=1:length(bla)
%    blace{in}=sprintf('%.0f',bla(in));
%end
%ax.XTickLabel=blace;
yel=[0 max(b.YData)+0.1*max(b.YData)];
ylim(yel)
grid on
xlabel(sprintf('residuals [%s]',multi))
t=text(ax.XLim(1)+0.05*abs(ax.XLim(2)-ax.XLim(1)),0.75*ax.YLim(2),...
       sprintf('N = %3.0f\nstd = %3.0f\nmed = %3.0f\navg = %3.0f\ngof = %3.0f',...
                    length(data),stdd,median(data),mean(data),gof));
hold on
pd = fitdist(data,'Normal');
xvals = b.XData; 
yvals = pdf(pd,xvals);
area = sum(b.YData)*diff(b.XData(1:2));
lain = plot(xvals,yvals*area,'LineWidth',2);
if gof > thresh | abs(mean(data)) > 2
    lain.LineStyle = '--';
end
hold off