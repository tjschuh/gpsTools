function sagegage(d,tmax,xyz,xyzn,v,vn)

% default plotting is off
defval('plt',0)
    
% C-DOG location [x,y,z] [m]
oceandep = 5225;
[x,y,z] = gps2dep([1.977967 -5.073198 3.3101016]*1e6,oceandep);
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

% add uncertainty to ship locations
% |std| = sqrt(xstd^2 + ystd^2 + zstd^2)
% assuming xstd = ystd, zstd = 2*xstd, |std| = 2 (from previous research)
xstd = 0.8165; ystd = 0.8165; zstd = 1.6330;
perturbs = [xstd ystd zstd];
permulti{1,1} = 1e-2; permulti{1,2} = 'cm'; 
d.x0 = d.x; d.y0 = d.y; d.z0 = d.z;
d.x = d.x0 + xstd.*randn(length(d.x0),1).*permulti{1,1};
d.y = d.y0 + ystd.*randn(length(d.y0),1).*permulti{1,1};
d.z = d.z0 + zstd.*randn(length(d.z0),1).*permulti{1,1};

% forward model: creating perturbed data
%hsr = sqrt((d.x0-xyzg(1)).^2 + (d.y0-xyzg(2)).^2 + (d.z0-xyzg(3)).^2);
hsr = sqrt((d.x-xyzg(1)).^2 + (d.y-xyzg(2)).^2 + (d.z-xyzg(3)).^2);
hst = hsr./vg;
% error between slant times and predicted slant times
differ = st-hst;
% calculate relative time difference (unitless quantity)
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

stda=std(differ(rows)*tmulti{1,1});

% plotting
figure(1)
clf
lwidth = 1.5;

% plot st and hst
ah(1)=subplot(3,3,[1 3]);
p(1)=plot(d.t,hst,'-','color',[0.8500 0.3250 0.0980],'LineWidth',lwidth);
hold on
p(2)=plot(d.t,st,'b--','LineWidth',lwidth);
hold off
datetick('x','HH')
xticklabels([])
cosmot(d.t)
ylim([3.25 9])
ylabel('slant range time [s]')
%title(sprintf('hst and st, ocean depth = %g m',oceandep))
title(sprintf('hst and st'))
legend({'hst','st'})

% plot absolute time differences
ah(2)=subplot(3,3,[4 6]);
p(3)=plot(d.t,differ*tmulti{1,1},'color',[0.8500 0.3250 0.0980],'LineWidth',lwidth);
hold on
p(4)=plot(d.t,movmean(differ*tmulti{1,1},3600,'omitnan'),'b','LineWidth',lwidth);
hold off
datetick('x','HH')
cosmot(d.t)
if nansum(differ) ~= 0
    ylim(1.1*[-max(abs(differ*tmulti{1,1})) max(abs(differ*tmulti{1,1}))])
end
title('Difference between st and hst')
ylabel(sprintf('slant range time [%s]',tmulti{1,2}))
xlabel('time [h]')
% if pval >= 1e-15
%     text(ah(2).XLim(1)+0.01*abs(ah(2).XLim(2)-ah(2).XLim(1)),0.9*ah(2).YLim(1),...
%          sprintf('dw = %.3f, p = %.3g',dw,pval));
% else
%     text(ah(2).XLim(1)+0.01*abs(ah(2).XLim(2)-ah(2).XLim(1)),0.9*ah(2).YLim(1),...
%          sprintf('dw = %.3f, p = 0',dw));
% end
% if pval >= pthresh
%     text(ah(2).XLim(2)-0.25*abs(ah(2).XLim(2)-ah(2).XLim(1)),0.9*ah(2).YLim(1),...
%          sprintf('ACCEPTED'))
% else
%     text(ah(2).XLim(2)-0.25*abs(ah(2).XLim(2)-ah(2).XLim(1)),0.9*ah(2).YLim(1),...
%          sprintf('REJECTED'))
% end

% plot absolute time differences
ah(3)=subplot(3,3,7);
try
    thresh2=500;
    [c,lain2,gof2,stda]=cosmoh(differ(rows)*tmulti{1,1},ah(3),thresh2,tmulti{1,2});
    c.FaceColor = [0.8500 0.3250 0.0980];
    lain2.Color = 'red';
    title('Absolute Time Differences')
catch
end

% plot relative time differences
ah(4)=subplot(3,3,8);
try
    thresh1=500;
    [b,lain1,gof1,stdr]=cosmoh(rel(rows),ah(4),thresh1,tmulti{1,2});        
    b.FaceColor = [0.400 0.6667 0.8431];
    lain1.Color = 'blue';
    title('Relative Time Differences')
    xlabel('residuals')
catch
end

ah(5)=subplot(3,3,9);
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
[lat,lon,~] = ecef2geodetic(wgs84,xyz(1),xyz(2),xyz(3));
[dogx,dogy,~] = deg2utm(lat,mod(lon,360));
box on
grid on
longticks([],2)
xl(2)=xlabel('easting [km]');
yl(2)=ylabel('northing [km]');
openup(ah(5),5,10);
openup(ah(5),6,10);
scatter((dogx-refx)/sclx,(dogy-refy)/scly,20,...
        'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r')
hold off

tt=supertit(ah([1 3]),sprintf('Distance from Truth = [%g %s %g %s %g %s] = |%3.3f %s|, True Sound Speed = %g m/s\nSound Speed Error = %g m/s, GPS Perturbations = +/-[%g %s %g %s %g %s], Depth = %g m',xyzn(1),xmulti{1,2},xyzn(2),xmulti{1,2},xyzn(3),xmulti{1,2},sqrt(xyzn(1)^2+xyzn(2)^2+xyzn(3)^2),xmulti{1,2},v0,vg-v0,xstd,permulti{1,2},ystd,permulti{1,2},zstd,permulti{1,2},oceandep));
tt.FontSize = 12;
movev(tt,0.325);

%figdisp(sprintf('experiment_%g_%g_%g',xyzn(1),xyzn(2),xyzn(3)),[],[],2,[],'epstopdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cosmot(t)
xlim([t(1) t(end)])
grid on
longticks([],3)

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
lain = plot(xvals,yvals*area,'LineWidth',2);
if gof > thresh | abs(mean(data)) > 2
    lain.LineStyle = '--';
end
hold off