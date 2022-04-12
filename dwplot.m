function dwplot(d,tmax,data,thresh)
% DWPLOT(data,thresh)
%
% take in Durbin-Watson test statistics and p-values, and standard
% deviations corresponding to perturbations in x,y,z from a true
% C-DOG location and generate 2D contour plots where p-value > thresh
%
% INPUT:
%
% d        x,y,z,t ship positions
% tmax     start and end times from d.t
% data     actual data containing x,y,z,dw,p-value,std
% thresh   p-value threshold [default: 0.05]
%
% Originally written by tschuh-at-princeton.edu, 04/05/2022
% Last modified by tschuh-at-princeton.edu, 04/12/2022

% to do
% use contourf on top of plots to show p > 0.05, p > 0.025, etc
% caxis([]) sets colorbar axis
% maybe add thresh value to title

% eventually make these inputs
v0 = 1500; dv = 0;

% p-value threshold
defval('thresh',0.05)

% size needs to equal 2n+1 x 2n+1 where n=x*y (mesh(x,y))
% in other words, reshape size = length(data)^(1/3)
errsz=nthroot(length(data),3);

figure(1)

% change this to a for loop
ah(1)=subplot(2,2,1);
ddxy = data(data(:,3) == 0,:);
xycornx = reshape(ddxy(:,1),errsz.*[1 1]);
xycorny = reshape(ddxy(:,2),errsz.*[1 1]);
errellipxy=reshape([ddxy(:,5) > thresh].*ddxy(:,6),errsz.*[1 1]);
% need to set 0s to NaNs to have colorbar ignore them
%errellipxy(errellipxy==0)=NaN;
imagesc([xycornx(1,1) xycornx(1,end)],[xycorny(1,1) xycorny(end,1)],errellipxy)
%title('y vs x')
xlabel('x [mm]')
ylabel('y [mm]')
axis tight equal
hold off
cosmo(ah(1),data(:,1),data(:,2))

ah(2)=subplot(2,2,2);
ddyz = data(data(:,1) == 0,:);
yzcornx = reshape(ddyz(:,2),errsz.*[1 1]);
yzcorny = reshape(ddyz(:,3),errsz.*[1 1]);
errellipyz=reshape([ddyz(:,5) > thresh].*ddyz(:,6),errsz.*[1 1]);
%errellipyz(errellipyz==0)=NaN;
imagesc([yzcornx(1,1) yzcornx(1,end)],[yzcorny(1,1) yzcorny(end,1)],errellipyz)
xlabel('y [mm]')
ylabel('z [mm]')
axis tight equal
hold off
cosmo(ah(2),data(:,2),data(:,3))

ah(3)=subplot(2,2,3);
ddxz = data(data(:,2) == 0,:);
xzcornx = reshape(ddxz(:,1),errsz.*[1 1]);
xzcorny = reshape(ddxz(:,3),errsz.*[1 1]);
errellipxz=reshape([ddxz(:,5) > thresh].*ddxz(:,6),errsz.*[1 1]);
%errellipxz(errellipxz==0)=NaN;
imagesc([xzcornx(1,1) xzcornx(1,end)],[xzcorny(1,1) xzcorny(end,1)],errellipxz)
xlabel('x [mm]')
ylabel('z [mm]')
axis tight equal
hold off
cosmo(ah(3),data(:,1),data(:,3))

% need to enble this to take in data without d.utm
% plot ship trajectory w/ C-DOG
% need xyz0 from gps2fwd if you want to plot C-DOG
ah(4)=subplot(2,2,4);
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
%wgs84 = wgs84Ellipsoid('meter');
%[lat,lon,~] = ecef2geodetic(wgs84,xyz0(1),xyz0(2),xyz0(3));
%[dogx,dogy,~] = deg2utm(lat,mod(lon,360));
box on
grid on
longticks([],2)
xl(2)=xlabel('easting [km]');
yl(2)=ylabel('northing [km]');
openup(ah(4),5,10);
openup(ah(4),6,10);
%scatter((dogx-refx)/sclx,(dogy-refy)/scly,10,...
%        'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r')
%[lat,lon,~] = ecef2geodetic(wgs84,xyzg(1),xyzg(2),xyzg(3));
%[fakex,fakey,~] = deg2utm(lat,mod(lon,360));
%scatter((fakex-refx)/sclx,(fakey-refy)/scly,10,...
%        'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','r')
hold off

tt=supertit(ah([1 2]),sprintf('True Sound Speed = %g m/s, Sound Speed Error = %g m/s\nGPS Perturbations = +/-[2 cm 2 cm 2 cm]',v0,dv));
tt.FontSize = 11;
movev(tt,0.2)
keyboard
% save figure as pdf
figdisp(sprintf('dwdata-plot'),[],[],2,[],'epstopdf')

function cosmo(ah,datax,datay)
a = colorbar;
a.Label.String = 'std [\mus]';
a.FontSize = 11;
%xlim([min(datax) max(datax)])
%ylim([min(datay) max(datay)])
ah.YDir = 'normal';    
grid on
longticks