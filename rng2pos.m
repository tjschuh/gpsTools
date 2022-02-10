function varargout = rng2pos(dxyz,st,vg,xyzg)
% [sol,hst] = RNG2POS(dxyz,st,vg,xyzg)
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
% xyzg          initial beacon location guess x,y,z
%
% OUTPUT:
%
% sol           calculated [x y z] beacon location [m] best-fitting st 
%
% EXAMPLE:
%
% [dxyz,sr,st,xyzg,vg]=gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'});
% sol=rng2pos(dxyz,st,vg,xyzg);
%
% Originally written by tschuh-at-princeton.edu, 02/07/2022
% Last modified by fjsimons-at-princeton.edu, 02/08/2022
% Last modified by tschuh-at-princeton.edu, 02/10/2022

% Problem Flow:
% compute hsr and hst from inital guess
% rmse=norm(st - hst)
% if rmse is bad, refine sol through least-squares inversion
% sol=sol+something
% use sol to get updated hst and run rmse=norm(st-hst)

% constant sound speed profile for now [m/s]
defval('vg',1500)

% Solution begins with the guess
sol = xyzg;

% The forward model is based on a guess
hsr = sqrt((dxyz(:,1)-sol(1)).^2 + (dxyz(:,2)-sol(2)).^2 + (dxyz(:,3)-sol(3)).^2);
% simple forward model (could have used GPS2RNG again, more complicated forward models, etc.)
hst = hsr./vg;
rmse = norm(st(~isnan(st)) - hst(~isnan(hst)));

% optional output
varns={sol,hst,rmse};
varargout=varns(1:nargout);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% try again
% Now the inversion needs to start

%while rmse > ??

% calculate sensitivity matrix G
% Gi1 = d(tpre)/dx, Gi2 = d(tpre)/dy, Gi3 = d(tpre)/dz
Gi1 = (vg.*(sol(1) - dxyz(:,1)))./sqrt((sol(1) - dxyz(:,1)).^2 + (sol(2) - dxyz(:,2)).^2 + (sol(3) - dxyz(:,3)).^2);
Gi2 = (vg.*(sol(2) - dxyz(:,2)))./sqrt((sol(1) - dxyz(:,1)).^2 + (sol(2) - dxyz(:,2)).^2 + (sol(3) - dxyz(:,3)).^2);
Gi3 = (vg.*(sol(3) - dxyz(:,3)))./sqrt((sol(1) - dxyz(:,1)).^2 + (sol(2) - dxyz(:,2)).^2 + (sol(3) - dxyz(:,3)).^2);

G = [Gi1 Gi2 Gi3];

% calculate difference between observational slant times and predicted slant times
% hst is right, but st might be wrong?
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
rmse = norm(st(~isnan(st)) - hst(~isnan(hst)));

%end
