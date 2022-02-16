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
% [st,xyz,v] = gps2fwd(d,tmax,[2e6 -4.5e6 3e6],1500);
% xyzmat = gps2inv(st,d,v,xyz);
%
% Originally written by tschuh-at-princeton.edu, 02/16/2022

% will clean this up later
% just need working version for talk

% need to get rid of rows with NaNs
% need to combine d.t and d.xyz into 1 matrix to remove all NaN rows
%d.t(any(isnan(d.t),2),:)=[];
d.xyz(any(isnan(d.xyz),2),:)=[];
st(any(isnan(st),2),:)=[];

% constant sound speed profile for now [m/s]
defval('vg',1500)

% guess the correct 
defval('xyzg',[2e6 -4.5e6 3e6])

v0 = vg;
xyz0 = xyzg;

xn = -0.1:0.01:0.1;
yn = -0.1:0.01:0.1;
zn = -0.1:0.01:0.1;

xyzmat = zeros(length(xn)*length(yn)*length(zn),4);
counter = 1;

for i=1:length(xn)
    for j=1:length(yn)
        for k=1:length(zn)
            sol(1) = xyzg(1) + xn(i);
            sol(2) = xyzg(2) + yn(j);
            sol(3) = xyzg(3) + xn(k);

            % The forward model is based on a guess
            hsr = sqrt((d.xyz(:,1)-sol(1)).^2 + (d.xyz(:,2)-sol(2)).^2 + (d.xyz(:,3)-sol(3)).^2);
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

% plotting
f=figure;
f.Position = [675 281 955 680];
sz=50;

scatter3(xyzmat(:,1),xyzmat(:,2),xyzmat(:,3),sz,xyzmat(:,4),'filled')
% change this to cm
xlabel(sprintf('x distance from truth [m]'))
ylabel(sprintf('y distance from truth [m]'))
zlabel(sprintf('z distance from truth [m]'))
%title(sprintf(''))
a = colorbar;
a.Label.String = 'rmse [s]';
a.FontSize = 11;
xlim([-0.11 0.11])
ylim([-0.11 0.11])
zlim([-0.11 0.11])
grid on
longticks
% add v0 label
% add xyz0/truth label
% also find min rmse and what xn,yn,zn it corresponds to

% also need to plot some slices from the cube plot

% also need to plot a histogram of the rmse values

% optional output
varns={xyzmat};
varargout=varns(1:nargout);