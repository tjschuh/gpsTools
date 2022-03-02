function varargout = gps2inv(st,d,vg,xyzg)
% xyzmat = GPS2INV(st,d,vg,xyzg)
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
% OUTPUT:
%
% xyzmat
%
% EXAMPLE:
%
% load Unit1234-camp.mat
% [st,xyz,v] = gps2fwd(d,tmax,[1.977967 -5.073198 3.3101016]*1e6,1500);
% xyzmat = gps2inv(st,d,v,xyz);
%
% Originally written by tschuh-at-princeton.edu, 02/16/2022

% will clean this up later
% just need working version for talk
% add noise to ship locations [+/- 10 cm]

% need to get rid of rows with NaNs
% need to combine d.t and d.xyz into 1 matrix to remove all NaN rows
%d.t(any(isnan(d.t),2),:)=[];
%d.xyz(any(isnan(d.xyz),2),:)=[];
%st(any(isnan(st),2),:)=[];
d.x(any(isnan(d.x),2),:)=[];
d.y(any(isnan(d.y),2),:)=[];
d.z(any(isnan(d.z),2),:)=[];
st(any(isnan(st),2),:)=[];

% constant sound speed profile for now [m/s]
defval('vg',1500)

% guess the correct location
defval('xyzg',[1.977967 -5.073198 3.3101016]*1e6)

v0 = vg;
xyz0 = xyzg;

% mesh/grid
% need to make these constants
% 0.01 --> 0.005
int = 0.01;
xn = -0.1:int:0.1;
yn = -0.1:int:0.1;
zn = -0.1:int:0.1;

% int = 1;
% xn = -10:int:10;
% yn = -10:int:10;
% zn = -10:int:10;

xyzmat = zeros(length(xn)*length(yn)*length(zn),4);
counter = 1;

for i=1:length(xn)
    for j=1:length(yn)
        for k=1:length(zn)
            sol(1) = xyzg(1) + xn(i);
            sol(2) = xyzg(2) + yn(j);
            sol(3) = xyzg(3) + xn(k);

            % The forward model is based on a guess
            hsr = sqrt((d.x-sol(1)).^2 + (d.y-sol(2)).^2 + (d.z-sol(3)).^2);
            % simple forward model (could have used GPS2FWD again, more complicated forward models, etc.)
            hst = hsr./vg;
            rmse = norm(st - hst);

            xyzmat(counter,1) = xn(i);
            xyzmat(counter,2) = yn(j);
            xyzmat(counter,3) = zn(k);
            xyzmat(counter,4) = rmse;

            counter = counter + 1;
        end
    end
end

% find what xn,yn,zn min(rmse) corresponds to
minerr = xyzmat(find(xyzmat(:,4) == min(xyzmat(:,4))),:);

% plotting
figure(1);
%ah(1)=subplot(2,2,[1 3]);
sz=50;
scatter3(100*xyzmat(:,1),100*xyzmat(:,2),100*xyzmat(:,3),sz,1000*xyzmat(:,4),'filled')
xlabel(sprintf('x distance from truth [cm]'))
ylabel(sprintf('y distance from truth [cm]'))
zlabel(sprintf('z distance from truth [cm]'))
title(sprintf('Quality of Inversion, v_{diff} = %g m/s\nxyz_{truth} = [%g m %g m %g m], min = [%g cm %g cm %g cm %g ms]',vg-1500,xyz0(1),xyz0(2),xyz0(3),100*minerr(1),100*minerr(2),100*minerr(3),1000*minerr(4)))
a = colorbar;
a.Label.String = 'rmse [ms]';
a.FontSize = 11;
xlim([110*min(xn) 110*max(xn)])
ylim([110*min(yn) 110*max(yn)])
zlim([110*min(zn) 110*max(zn)])
grid on
longticks

%figdisp(sprintf('gps2inv-cube'),[],'',2,[],'epstopdf')
%close

% also need to plot slice at z = 0 because we know ocean depth very well
xymat = xyzmat(find(xyzmat(:,3) == 0),:);

g=figure;
%ah(2)=subplot(2,2,2);
scatter(100*xymat(:,1),100*xymat(:,2),sz,1000*xymat(:,4),'filled')
xlabel(sprintf('x distance from truth [cm]'))
ylabel(sprintf('y distance from truth [cm]'))
title(sprintf('Quality of Inversion, v_{diff} = %g m/s\nslice at z = 0 cm',vg-1500))
a = colorbar;
a.Label.String = 'rmse [ms]';
a.FontSize = 11;
xlim([110*min(xn) 110*max(xn)])
ylim([110*min(yn) 110*max(yn)])
grid on
longticks

%figdisp(sprintf('gps2inv-slice'),[],'',2,[],'epstopdf')
%close

% plot only rows where rmse < something
thresh = 15;
xyzmatell = xyzmat(find(xyzmat(:,4) <= (thresh/1000)),:);

%ah(3)=subplot(2,2,4);
h=figure;
scatter3(100*xyzmatell(:,1),100*xyzmatell(:,2),100*xyzmatell(:,3),sz,1000*xyzmatell(:,4),'filled')
xlabel(sprintf('x distance from truth [cm]'))
ylabel(sprintf('y distance from truth [cm]'))
zlabel(sprintf('z distance from truth [cm]'))
title(sprintf('rmse < %g ms, v_{diff} = %g m/s\nxyz_{truth} = [%g m %g m %g m], min = [%g cm %g cm %g cm %g ms]',thresh,vg-1500,xyz0(1),xyz0(2),xyz0(3),100*minerr(1),100*minerr(2),100*minerr(3),1000*minerr(4)))
a = colorbar;
a.Label.String = 'rmse [ms]';
a.FontSize = 11;
xlim([110*min(xn) 110*max(xn)])
ylim([110*min(yn) 110*max(yn)])
zlim([110*min(zn) 110*max(zn)])
grid on
longticks

%figdisp(sprintf('gps2inv-ellipsoid'),[],'',2,[],'epstopdf')
%close

% need to plot a histogram of the rmse values
% not centered on 0 because rmse >= 0
%histogram(xyzmat(:,4))
% need to get acoustic curve from dog2mat.m

% optional output
varns={xyzmat};
varargout=varns(1:nargout);