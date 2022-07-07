function dwplot(d,data,thresh)
% DWPLOT(d,data,thresh)
%
% take in Durbin-Watson test statistics and p-values, and standard
% deviations corresponding to perturbations in x,y,z from a true
% C-DOG location and generate 2D contour plots where p-value > thresh
%
% INPUT:
%
% d        x,y,z,t ship positions
% data     actual data containing x,y,z,dw,p-value,std
% thresh   p-value threshold [default: 0.05]
%
% OUTPUT:
%
% 2x2 subplot with error ellipses and ship trajectory
%    
% EXAMPLES:
%
% run gps2syn.m first
% dwplot(d,xyzdwp,0.05)
%    
% Originally written by tschuh-at-princeton.edu, 04/05/2022
% Last modified by tschuh-at-princeton.edu, 06/06/2022

% to do:
% need to make Perturbations an input somehow
% add ocean height to title --> make ocean height input somehow
% make this more general so can be called by gps2syn
% make different cases within for loop to plot
% results in different ways

% eventually make these inputs
v0 = 1500; dv = 0;

% p-value threshold
defval('thresh',0.05)

% size needs to equal 2n+1 x 2n+1 where n=x*y (mesh(x,y))
% in other words, reshape size = length(data)^(1/3)
errsz=nthroot(length(data),3);

figure(gcf)
clf

% make lookup table
index = [3 1 2; 1 2 3; 2 1 3];
labels = {'x' 'y'; 'y' 'z'; 'x' 'z'};

% plot error ellipse for each coordinate pair
for i=1:3
    ah(i)=subplot(2,2,i);
    % hold 1 variable constant 0 and work with other 2
    dd = data(data(:,index(i,1)) == 0,:);
    % define corners for mesh grid
    cornx = reshape(dd(:,index(i,2)),errsz.*[1 1]);
    corny = reshape(dd(:,index(i,3)),errsz.*[1 1]);
    % error ellipse are rows where p-value > thresh
    errellip = reshape([dd(:,5) > thresh].*dd(:,6),errsz.*[1 1]);
    % need to set 0s to NaNs and then cosmo will turn them white
    errellip(errellip==0) = NaN;
    % actually make colormap of error ellipse
    isc = imagesc([cornx(1,1) cornx(1,end)],[corny(1,1) corny(end,1)],errellip);
    xlabel(sprintf('%s [mm]',labels{i,1}))
    ylabel(sprintf('%s [mm]',labels{i,2}))
    tint=1;
    xticks(min(data(index(i,2),:)):tint:max(data(index(i,2),:)));
    yticks(min(data(index(i,3),:)):tint:max(data(index(i,3),:)));
    axis tight equal
    hold on
    cosmo(ah(i),isc,errellip)
    % add contour line(s)
    % error contour(s) to overlay on errellip
    errcont = reshape((dd(:,5)),errsz.*[1 1]);
    % actually overlay contour line(s)
    contour(cornx,corny,errcont,[thresh thresh],'k','LineWidth',1.25);
    % caxis([ ]) sets limits on colorbar
    caxis([min(errellip,[],'all') max(errellip,[],'all')])
    hold off
end

fourplt = 1;
if fourplt == 1
% plot ship trajectory w/ C-DOG
% need xyz0 from gps2fwd if you want to plot C-DOG location
ah(4)=subplot(2,2,4);
if isfield(d,'utme') == 1
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
else
    sclx=1000;
    scly=1000;
    plot(d.x/sclx,d.y/scly,'LineWidth',1,'Color','k')
    hold on
    box on
    grid on
    longticks([],2)
    xl(2)=xlabel('easting [km]');
    yl(2)=ylabel('northing [km]');
    openup(ah(4),5,10);
    openup(ah(4),6,10);
    scatter(0,0,20,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','r')
    hold off
end
end

% need to make Perturbations an input somehow
tt=supertit(ah([1 2]),sprintf('True Sound Speed = %g m/s, Sound Speed Error = %g m/s, Ocean Depth = 5225 m\nGPS Perturbations = +/-[0.8165 cm 0.8165 cm 1.633 cm], p-value > %g',v0,dv,thresh));
tt.FontSize = 11;
movev(tt,0.2)

% save figure as pdf
figdisp(sprintf('dwdata-plot-%03d',1000*thresh),[],[],2,[],'epstopdf')

function cosmo(ah,isc,iscdata)
a = colorbar;
a.Label.String = 'std [\mus]';
a.FontSize = 11;
ah.YDir = 'normal';    
grid on
longticks
% set NaNs to white when plotting
if ndims(iscdata) == 2
    set(isc,'AlphaData',~isnan(iscdata))
elseif ndims(iscdata) == 3
    set(isc,'AlphaData',~isnan(iscdata(:,:,1)))
end