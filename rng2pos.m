function varargout = rng2pos(dxyz,st,v,xyzg)
% [sol,hst] = GPS2INV(dxyz,st,vguess,xyzg)
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
% vguess        trial sound speed [m/s]
% xyzg          initial beacon location guess x,y,z
%
% OUTPUT:
%
% sol           calculated [x y z] beacon location [m] best-fitting st 
%
% EXAMPLE:
%
% [st,dxyz,sr,xyzg,vg]=...
%     gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'});
% sol=rng2pos(dxyz,st,vg,xyzg);
%
% Originally written by tschuh-at-princeton.edu, 02/07/2022
% Last modified by tschuh-at-princeton.edu, 02/08/2022
% Last modified by fjsimons-at-princeton.edu, 02/08/2022

% To Do:
% figure out how to speed up
% what's a good phi to use?
% what's another metric to determine when to stop iterating?
% add noise to slant times and receiver positions (covariance matrix)

% initial guess
defval('sol',[2e6 -4.5e6 3e6])

% constant sound speed profile for now [m/s]
defval('v',1500)

% combine station positions with time arrivals
% actual data observations
dxyzt = [dxyz st];

% Solution begins with the guess
sol=xyzg;

% Beginning is bad
% rmse=norm(st);
% while norm bad
% change your sol
% sol=sol+something
% The forward model is based on a guess
hsr=sqrt((dxyz(:,1)-sol(1)).^2+(dxyz(:,2)-sol(2)).^2+(dxyz(:,3)-sol(3)).^2);
% simple forward model (could have used GPS2RNG again, more complicated
% forward models, etc.)
hst = hsr./v;
rmse=norm(st(~isnan(st))-hst(~isnan(st)));
% try again

% optional output
varns={sol,hst,rmse};
varargout=varns(1:nargout);

return


% Now the inversion needs to start

% calculate sensitivity matrix G
% Gi1 = d(tpre)/dx, Gi2 = d(tpre)/dy,
% Gi3 = d(tpre)/dz, Gi4 = d(tpre)/dt
Gi1 = (v.*(xyzg(1) - dxyzt(:,1)))./sqrt((xyzg(1) - dxyzt(:,1)).^2 + (xyzg(2) - dxyzt(:,2)).^2 + (xyzg(3) - dxyzt(:,3)).^2);
Gi2 = (v.*(xyzg(2) - dxyzt(:,2)))./sqrt((xyzg(1) - dxyzt(:,1)).^2 + (xyzg(2) - dxyzt(:,2)).^2 + (xyzg(3) - dxyzt(:,3)).^2);
Gi3 = (v.*(xyzg(3) - dxyzt(:,3)))./sqrt((xyzg(1) - dxyzt(:,1)).^2 + (xyzg(2) - dxyzt(:,2)).^2 + (xyzg(3) - dxyzt(:,3)).^2);

G = [Gi1 Gi2 Gi3];

% calculate difference between observed travel time and predicted travel time
deltad = dxyzt(:,4) - tpre;

% chi^2 = (deltad - Gdeltam)^2
% this is what we are trying to minimize
% so we differentiate and set equal to zero and solve for deltam
deltam = inv(transpose(G)*G)*transpose(G)*deltad;

% deltam contains our adjustments to our predicted solution (x,y,z,t)
% add deltam to our initial xyzg to get a new xyzg m
m = xyzg + deltam';

% calculate phi (sum of squares of residuals)
% to determine how good how prediction is
phi = (sum(deltad))^2;

% calculate uncertainty between tru and prediction
% this probably wrong
phi2 = sqrt((m(1) - tru(1))^2 + (m(2) - tru(2))^2 + (m(3) - tru(3))^2 + (m(4) - tru(4))^2);

% go again and keep refining until phi < thresh and we are happy
thresh = 1e10;
while phi2 >= 1e5
    tpre = m(4) + sqrt((m(1) - dxyzt(:,1)).^2 + (m(2) - dxyzt(:,2)).^2 + (m(3) - dxyzt(:,3)).^2)/v;

    Gi1 = (v.*(m(1) - dxyzt(:,1)))./sqrt((m(1) - dxyzt(:,1)).^2 + (m(2) - dxyzt(:,2)).^2 + (m(3) - dxyzt(:,3)).^2);
    Gi2 = (v.*(m(2) - dxyzt(:,2)))./sqrt((m(1) - dxyzt(:,1)).^2 + (m(2) - dxyzt(:,2)).^2 + (m(3) - dxyzt(:,3)).^2);
    Gi3 = (v.*(m(3) - dxyzt(:,3)))./sqrt((m(1) - dxyzt(:,1)).^2 + (m(2) - dxyzt(:,2)).^2 + (m(3) - dxyzt(:,3)).^2);
    Gi4 = ones(length(dxyzt),1);

    G = [Gi1 Gi2 Gi3 Gi4];

    deltad = dxyzt(:,4) - tpre;

    deltam = inv(transpose(G)*G)*transpose(G)*deltad;

    m = m + deltam';

    phi = (sum(deltad))^2;

    phi2 = sqrt((m(1) - tru(1))^2 + (m(2) - tru(2))^2 + (m(3) - tru(3))^2 + (m(4) - tru(4))^2);
end

sol = m;

