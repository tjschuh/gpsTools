function sagegage(d,tmax,xyz,xyzn,v,vn)
% SAGEGAGE(d,tmax,xyz,xyzn,v,vn)
%
% Produces two time series of distances to a moving ship position (input,
% from data), one from a nominal subsurface point (the "truth") and another
% from a perturbation to that location, to assess error structure. The
% locations of the ship trajectory are corrupted by noise as part of this.
% 
% INPUT:
%
% d            single dataset with x,y,z,t ship locations
% tmax         first and last time in d [2x1 datetime array]
% xyz          subsurface beacon location [x y z] [m]
% xyzn         perturbations to beacon location [x y z] [m]
% v            sound speed [m/s]
% vn           perturbation in sound speed [m/s]
%
% load Unit1234-camp.mat
% 
% sagegage(d,tmax,[],[1 0.5 0.5],[],[])
% 
% Last modified by tschuh-at-princeton.edu,  05/18/2022
% Last modified by fjsimons-at-alum.mit.edu,  05/18/2022

% default plotting is off
defval('plt',0)
    
% constant sound speed profile for now [m/s]
defval('v',1500)
% velocity noise
defval('vn',0)

% unit multipliers across time and space
tmulti{1,1} = 1e6; tmulti{1,2} = '\mu';
xmulti{1,1} = 1e-3; xmulti{1,2} = 'mm';
permulti{1,1} = 1e-2; permulti{1,2} = 'cm'; 

% First subsurface position, the "C-DOG" location [x,y,z] [m]
oceandep = 5225;
% Compared to just about where we dropped it
[x,y,z] = gps2dep([1.977967 -5.073198 3.3101016]*1e6,oceandep);
defval('xyz',[x y z])

% Second subsfurface position, a perturbation to xyz to explore around C-DOG
defval('xyzn',[10 10 10])

% call GPS2FWD to generate perfect slant times data to the C-DOG
[st,xyz0,v0] = gps2fwd(d,tmax,xyz,v);

% find the rows that don't contain NaNs
rowsx=find(~isnan(d.x));
rowsy=find(~isnan(d.y));
rowsz=find(~isnan(d.z));
% Rows where none of three entries are NaN
rows=unique([rowsx;rowsy;rowsz]);
% So that needs to be the same as 
% rows=find(~isnan(st));

% add xyzn to xyz0 to get perturbed C-DOG location
xyzg = xyz0 + xyzn*xmulti{1,1};
% add vn to v0 to get perturbed sound speed
vg = v0 + vn;

% add uncertainty to ship locations
% |std| = sqrt(xstd^2 + ystd^2 + zstd^2)
% assuming xstd = ystd, zstd = 2*xstd, |std| = 2 (from previous research)
xstd = 0.8165; ystd = 0.8165; zstd = 1.6330;
perturbs = [xstd ystd zstd];
d.x0 = d.x; d.y0 = d.y; d.z0 = d.z;
d.x = d.x0 + xstd.*randn(length(d.x0),1).*permulti{1,1};
d.y = d.y0 + ystd.*randn(length(d.y0),1).*permulti{1,1};
d.z = d.z0 + zstd.*randn(length(d.z0),1).*permulti{1,1};
% forward model: creating perturbed data
%hsr = sqrt((d.x0-xyzg(1)).^2 + (d.y0-xyzg(2)).^2 + (d.z0-xyzg(3)).^2);
hsr = sqrt((d.x-xyzg(1)).^2 + (d.y-xyzg(2)).^2 + (d.z-xyzg(3)).^2);
hst = hsr./vg;
% This needs to be the same as GPS2FWD

% difference between the slant times from C-DOG to ship 
% and the slant times from perturbed C-DOG to perturbed ship
differ = st-hst;

% calculate relative time difference (unitless quantity)
rel = differ./st;
% calculate rms(st-hst)
rmse = norm(differ(rows));
% Same as  norm(differ(~isnan(differ))) 

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

stda=std(differ(rows)*tmulti{1,1});

% plotting
figure(1)
clf
lwidth = 0.5;

% plot st and hst
ah(1)=subplot(3,3,[1 3]);
p(1)=plot(d.t,hst,'-','color',[0.8500 0.3250 0.0980],'LineWidth',3*lwidth);
hold on
p(2)=plot(d.t,st,'k--','LineWidth',2*lwidth);
hold off
% Must set ticks every n hours
datetick('x','HH')
%xticklabels([])
nolabels(gca,1)
cosmot(d.t)
ylim([3.25 9])
yl(1)=ylabel('slant range time [s]');
%title(sprintf('hst and st, ocean depth = %g m',oceandep))
t(1)=title(sprintf('hst and st'));
axes(ah(1))
bx=text(737958,8.5,sprintf('%s = [%g %g %g] mm','\Delta',xyzn));
set(bx,'FontSize',6)

leg(1)=legend({'pred','obs'});
movev(leg(1),.035)
moveh(leg(1),.025)


% plot absolute time differences
ah(2)=subplot(3,3,[4 6]);
p(3)=plot(d.t,differ*tmulti{1,1},'color',[0.8500 0.3250 0.0980],'LineWidth',lwidth);
hold on
p(4)=plot(d.t,movmean(differ*tmulti{1,1},3600,'omitnan'),'y','LineWidth',2*lwidth);
hold off
leg(2)=legend({'pred-obs','1-hr average'},'Location','NorthWest');
movev(leg(2),.085)
moveh(leg(2),-.04)

% Must set ticks every n hours
datetick('x','HH')

cosmot(d.t)
if nansum(differ) ~= 0
    ylim(1.1*[-max(abs(differ*tmulti{1,1})) max(abs(differ*tmulti{1,1}))])
end
t(2)=title('Difference between st and hst');%
%yl(2)=ylabel(sprintf('slant range time [%ss]',tmulti{1,2}));
yl(2)=ylabel(sprintf('obs-pred (x1e-6)'));
set(ah(2),'YaxisLocation','right')

set(yl(2),'Interpreter','TeX')
xl(1)=xlabel('time [h]');
if pval <= 10^6*eps
  pval=0;
end
tx(1)=text(ah(2).XLim(1)+0.01*abs(ah(2).XLim(2)-ah(2).XLim(1)),0.9*ah(2).YLim(1),...
     sprintf('dw = %.3f, p = %.3g',dw,pval));
decx={'ACCEPTED','REJECTED'};
tx(2)=text(ah(2).XLim(2)-0.015*abs(ah(2).XLim(2)-ah(2).XLim(1)),0.9*ah(2).YLim(1),...
     decx{2-[pval>=pthresh]},'HorizontalAlignment','right');

% plot absolute time differences
ah(3)=subplot(3,3,7);
try
  thresh2=500;
  [c,lain2,gof2,stda,ttt{1}]=cosmoh(differ(rows)*tmulti{1,1},ah(3),thresh2,tmulti{1,2});
  c.FaceColor = [0.8500 0.3250 0.0980];
  lain2.Color = 'red';
  t(3)=title(sprintf('obs-pred (x1e-6)'));
  xl(2)=xlabel(sprintf('obs-pred (x1e-6)'));
end
nolabels(ah(3),2)

% plot relative time differences
ah(4)=subplot(3,3,8);
try
  thresh1=500;
  [b,lain1,gof1,stdr,ttt{2}]=cosmoh(1e6*rel(rows),ah(4),thresh1,tmulti{1,2});        
  b.FaceColor = [0.400 0.6667 0.8431];
  lain1.Color = 'blue';
  t(4)=title('Relative Time Differences');
  xl(3)=xlabel(sprintf('obs-pred (relative, x1e-6)'));
end
nolabels(ah(4),2)

ah(5)=subplot(3,3,9);
skp=1000;
zex=d.utme(1:skp:end);
zwi=d.utmn(1:skp:end);
refx=min(d.utme);
refy=min(d.utmn);
sclx=1000;
scly=1000;
plot((d.utme-refx)/sclx,(d.utmn-refy)/scly,'LineWidth',1/2,'Color','k');
hold on
scatter((zex-refx)/sclx,(zwi-refy)/scly,5,...
        'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')
wgs84 = wgs84Ellipsoid('meter');
[lat,lon,~] = ecef2geodetic(wgs84,xyz(1),xyz(2),xyz(3));
warning off MATLAB:nargchk:deprecated
[dogx,dogy,~] = deg2utm(lat,mod(lon,360));
warning on MATLAB:nargchk:deprecated

box on
grid on
longticks([])
xl(4)=xlabel('easting [km]');
yl(3)=ylabel('northing [km]');
set(gca,'YAxisLocation','right')
openup(ah(5),5,10);
openup(ah(5),6,10);
scatter((dogx-refx)/sclx,(dogy-refy)/scly,20,...
        'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r')
hold off

tt=supertit(ah([1 3]),sprintf('Distance from Truth = [%g %s %g %s %g %s] = |%3.3f %s|, True Sound Speed = %g m/s\nSound Speed Error = %g m/s, GPS Perturbations = +/-[%g %s %g %s %g %s], Depth = %g m',xyzn(1),xmulti{1,2},xyzn(2),xmulti{1,2},xyzn(3),xmulti{1,2},sqrt(xyzn(1)^2+xyzn(2)^2+xyzn(3)^2),xmulti{1,2},v0,vg-v0,xstd,permulti{1,2},ystd,permulti{1,2},zstd,permulti{1,2},oceandep));
tt.FontSize = 12;
movev(tt,0.325);

%figdisp(sprintf('experiment_%g_%g_%g',xyzn(1),xyzn(2),xyzn(3)),[],[],2,[],'epstopdf')

% sage gage
set(ttt{1},'FontSize',6)
set(ttt{2},'FontSize',6)
set(tx,'FontSize',6)
set(ah,'FontSize',7)
set([xl yl],'FontSize',7)
delete([tt t(1) t(2) t(3) t(4)])
movev(ah(2),0.05)
movev(ah(3:5),0.03);


figdisp([],[],[],2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cosmot(t)
try
  xlim([t(1) t(end)])
catch
  xlim(datenum([t(1) t(end)]))
end
grid on
longticks([],3)

function [b,lain,gof,stdd,t]=cosmoh(data,ax,thresh,multi)
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
longticks([])
stdd=std(data);
nstd=3;
yel=[0 max(b.YData)+0.1*max(b.YData)];
ylim(yel)
grid on
xlabel(sprintf('residuals [%s]',multi))
t=text(ax.XLim(1)+0.05*abs(ax.XLim(2)-ax.XLim(1)),0.75*ax.YLim(2),...
       sprintf('N = %3.0f\nstd = %3.2f\nmed = %3.0f\navg = %3.0f',...
                    length(data),stdd,median(data),mean(data)));
hold on
pd = fitdist(data,'Normal');
xvals = b.XData; 
yvals = pdf(pd,xvals);
area = sum(b.YData)*diff(b.XData(1:2));
lain = plot(xvals,yvals*area,'LineWidth',1);
if gof > thresh | abs(mean(data)) > 2
    lain.LineStyle = '--';
end
hold off
