function varargout = rng2pos(dxyz,st,vg,xyzg,expnum)
% [sol,hst,rmse] = RNG2POS(dxyz,st,vg,xyzg,expnum)
%
% Given Precise Point Position time series of a certain configuration
% (e.g. an average four-station set like from GPS2RNG), and a corresponding
% sequence of "observed" slant times, do a least-squares inversion to
% calculate the predicted location of a beacon on the seafloor
%
% INPUT:
%
% dxyz          time series of "source" locations, e.g. an average GNSS series
% st            slant times/arrival times of the signals from the
%               "source" to the receiver whose position we want to find
%               here. These are from the DOGS, or, synthetic, from GPS2RNG
% vg            sound speed guess [m/s]
% vn            uncertainty in sound speed [m/s]
% xyzg          initial beacon location guess [x y z] [m]
% xyzn          uncertainty in xyz [m]
% expnum        number of experiments you want to run [default: 1000]
%
% OUTPUT:
%
% sol           calculated [x y z] beacon location [m] best-fitting st 
% hst
% rmse
%
% EXAMPLE:
%
% [st,dxyz,sr,xyzg,vg]=gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'});
% [st,dxyz,sr,xyzg,vg]=gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'},...
%   [],[2e6 -4.5e6 3e6]);
% [sol,hst,rmse]=rng2pos(dxyz,st,vg,xyzg,[]);
%
% [st,dxyz,sr,xyzg,vg]=gps2rng({'spiral.mat'},[],[0 0 0]);
% [sol,hst,rmse]=rng2pos(dxyz,st,vg,xyzg,[]);
%
% Originally written by tschuh-at-princeton.edu, 02/07/2022
% Last modified by fjsimons-at-princeton.edu, 02/08/2022
% Last modified by tschuh-at-princeton.edu, 02/14/2022

% Idea Behind Code:
% given a true beacon location using a forward model
% guess a perturbed location and perturbed velocity
% and use the true slant times to see just how bad
% our perturbed guesses were

% add true noise to both slant times and dxyz positions

% need to get rid of rows with NaNs
dxyz(any(isnan(dxyz),2),:)=[];
st(any(isnan(st),2),:)=[];

% constant sound speed profile for now [m/s]
defval('vg',1500)
v0 = vg;
xyz0 = xyzg;

% generate a random decimal value between -10 m/s and 10 m/s
%defval('vn',randi([-10 9]) + rand)

% generate a random decimal value between -10 cm and 10 cm
%defval('xyzn',[(randi(201)-101)./1000 (randi(201)-101)./1000 (randi(201)-101)./1000])

% number of trials to run
defval('expnum',1000)

vmat = zeros(expnum,6);

for i=1:expnum
    % generate a random decimal value between -10 m/s and 10 m/s
    vn = randi([-10 9]) + rand;

    % add noise to vg
    vg = v0 + vn;

    % generate a random decimal value between -10 cm and 10 cm
    xyzn = [(randi(201)-101)./1000 (randi(201)-101)./1000 (randi(201)-101)./1000];

    % take magnitude of xyzn
    mag = sqrt(xyzn(1)^2 + xyzn(2)^2 + xyzn(3)^2);

    % add noise to xyzg
    %dxyz = dxyz + xyzn;
    xyzg = xyz0 + xyzn;

    % Solution begins with the guess
    sol = xyzg;

    % The forward model is based on a guess
    hsr = sqrt((dxyz(:,1)-sol(1)).^2 + (dxyz(:,2)-sol(2)).^2 + (dxyz(:,3)-sol(3)).^2);
    % simple forward model (could have used GPS2RNG again, more complicated forward models, etc.)
    hst = hsr./vg;
    rmse = norm(st - hst);

    vmat(i,1) = vn;
    vmat(i,2) = xyzn(1);
    vmat(i,3) = xyzn(2);
    vmat(i,4) = xyzn(3);
    vmat(i,5) = mag;
    vmat(i,6) = rmse;
end

% plot rmse with a color bar and vn on x-axis, mag(xyzn) on y-axis
% this works, but we next want to complete inversion part to reduce
% rmse values and do this again
% may also want to change how the noise is being incorporated, I may have done this wrong
f=figure;
f.Position = [675 281 955 680];
sz=50;

subplot(2,2,1)
scatter(vmat(:,1),vmat(:,2),sz,vmat(:,6),'filled')
xlabel(sprintf('velocity perturbation from %i [m/s]',v0))
%ylabel(sprintf('distance from truth [%i %i %i] [m]',x0,y0,z0))
ylabel(sprintf('x distance from truth [m]'))
%title(sprintf(''))
ylim([-0.1 0.1])
cosmo1

subplot(2,2,2)
scatter(vmat(:,1),vmat(:,3),sz,vmat(:,6),'filled')
xlabel(sprintf('velocity perturbation from %i [m/s]',v0))
%ylabel(sprintf('distance from truth [%i %i %i] [m]',x0,y0,z0))
ylabel(sprintf('y distance from truth [m]'))
%title(sprintf(''))
ylim([-0.1 0.1])
cosmo1

subplot(2,2,3)
scatter(vmat(:,1),vmat(:,4),sz,vmat(:,6),'filled')
xlabel(sprintf('velocity perturbation from %i [m/s]',v0))
%ylabel(sprintf('distance from truth [%i %i %i] [m]',x0,y0,z0))
ylabel(sprintf('z distance from truth [m]'))
%title(sprintf(''))
ylim([-0.1 0.1])
cosmo1

subplot(2,2,4)
scatter(vmat(:,1),vmat(:,5),sz,vmat(:,6),'filled')
xlabel(sprintf('velocity perturbation from %i [m/s]',v0))
%ylabel(sprintf('distance from truth [%i %i %i] [m]',x0,y0,z0))
ylabel(sprintf('absolute distance from truth [m]'))
%title(sprintf(''))
ylim([0 0.2])
cosmo1

figdisp(sprintf('rng2pos-plt1'),[],'',2,[],'epstopdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vn = randi([-10 9]) + rand;
vg = v0 + vn;

xyzmat = zeros(expnum,4);

for i=1:expnum
xyzn = [(randi(201)-101)./1000 (randi(201)-101)./1000 (randi(201)-101)./1000];
xyzg = xyz0 + xyzn;
sol = xyzg;
hsr = sqrt((dxyz(:,1)-sol(1)).^2 + (dxyz(:,2)-sol(2)).^2 + (dxyz(:,3)-sol(3)).^2);
hst = hsr./vg;
rmse = norm(st - hst);
xyzmat(i,1) = xyzn(1);
xyzmat(i,2) = xyzn(2);
xyzmat(i,3) = xyzn(3);
xyzmat(i,4) = rmse;
end

g=figure;
g.Position = [675 281 955 680];
sz=50;

subplot(2,2,1)
scatter(xyzmat(:,1),xyzmat(:,2),sz,xyzmat(:,4),'filled')
xlabel(sprintf('x distance from truth [m]'))
ylabel(sprintf('y distance from truth [m]'))
%title(sprintf(''))
cosmo2

subplot(2,2,2)
scatter(xyzmat(:,1),xyzmat(:,3),sz,xyzmat(:,4),'filled')
xlabel(sprintf('x distance from truth [m]'))
ylabel(sprintf('z distance from truth [m]'))
%title(sprintf(''))
cosmo2

subplot(2,2,3)
scatter(xyzmat(:,2),xyzmat(:,3),sz,xyzmat(:,4),'filled')
xlabel(sprintf('y distance from truth [m]'))
ylabel(sprintf('z distance from truth [m]'))
%title(sprintf(''))
cosmo2

% would be cool to plot trajectory of ship
%subplot(2,2,4)
%scatter()

figdisp(sprintf('rng2pos-plt2'),[],'',2,[],'epstopdf')

% optional output
varns={sol,hst,rmse};
varargout=varns(1:nargout);

function cosmo1()
a = colorbar;
a.Label.String = 'rmse [s]';
a.FontSize = 11;
xlim([-10 10])
grid on
longticks

function cosmo2()
a = colorbar;
a.Label.String = 'rmse [s]';
a.FontSize = 11;
xlim([-0.1 0.1])
ylim([-0.1 0.1])
grid on
longticks

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inversion part works, but will probably become a different function

% Inversion Problem Flow:
% compute hsr and hst from inital guess
% rmse=norm(st - hst)
% if rmse is bad, refine sol through least-squares inversion
% sol=sol+something
% use sol to get updated hst
% rmse=norm(st-hst)
% repeat if needed

% try again
% Now the inversion needs to start

% how small does rmse need to be?
while rmse > .1
    % calculate sensitivity matrix G
    % Gi1 = d(tpre)/dx, Gi2 = d(tpre)/dy, Gi3 = d(tpre)/dz
    Gi1 = (sol(1) - dxyz(:,1))./(vg.*sqrt((sol(1) - dxyz(:,1)).^2 + (sol(2) - dxyz(:,2)).^2 + (sol(3) - dxyz(:,3)).^2));
    Gi2 = (sol(2) - dxyz(:,2))./(vg.*sqrt((sol(1) - dxyz(:,1)).^2 + (sol(2) - dxyz(:,2)).^2 + (sol(3) - dxyz(:,3)).^2));
    Gi3 = (sol(3) - dxyz(:,3))./(vg.*sqrt((sol(1) - dxyz(:,1)).^2 + (sol(2) - dxyz(:,2)).^2 + (sol(3) - dxyz(:,3)).^2));
    %Gi4 = -sqrt((sol(1) - dxyz(:,1)).^2 + (sol(2) - dxyz(:,2)).^2 + (sol(3) - dxyz(:,3)).^2)/(vg^2);

    G = [Gi1 Gi2 Gi3];

    % calculate difference between observational slant times and predicted slant times
    deltad = st - hst;

    % chi^2 = (deltad - Gdeltam)^2
    % this is what we are trying to minimize
    % so we differentiate and set equal to zero and solve for deltam
    deltam = inv(transpose(G)*G)*transpose(G)*deltad;

    % deltam contains our adjustments to our predicted solution (x,y,z)
    % add deltam to our initial xyzg to get a new sol
    sol = xyzg + deltam';

    % compute hsr and hst again
    hsr = sqrt((dxyz(:,1)-sol(1)).^2 + (dxyz(:,2)-sol(2)).^2 + (dxyz(:,3)-sol(3)).^2);
    hst = hsr./vg;

    % calculate norm again
    rmse = norm(st - hst)
end

% optional output
varns={sol,hst,rmse};
varargout=varns(1:nargout);