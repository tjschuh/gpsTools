function varargout = rng2pos(dxyz,st,vg,vn,xyzg,xyzn)
% [sol,hst,rmse] = RNG2POS(dxyz,st,vg,vn,xyzg,xyzn)
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
%
% OUTPUT:
%
% sol           calculated [x y z] beacon location [m] best-fitting st 
%
% EXAMPLE:
%
% [st,dxyz,sr,xyzg,vg]=gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'});
% [st,dxyz,sr,xyzg,vg]=gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'},...
%   [],[2e6 -4.5e6 3e6]);
% for i=1:10    
% [sol,hst,rmse]=rng2pos(dxyz,st,vg,[],xyzg,[]);
% end
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

% generate a random decimal value between -10 m/s and 10 m/s
defval('vn',randi([-10 9]) + rand)

% add noise to vg
vg = vg + vn;

% generate a random decimal value between -10 cm and 10 cm
defval('xyzn',[(randi(201)-101)./1000 (randi(201)-101)./1000 (randi(201)-101)./1000])

% take magnitude of xyzn
mag = sqrt(xyzn(1)^2 + xyzn(2)^2 + xyzn(3)^2);

% add noise to xyzg
xyzg = xyzg + xyzn;

% Solution begins with the guess
sol = xyzg;

% The forward model is based on a guess
hsr = sqrt((dxyz(:,1)-sol(1)).^2 + (dxyz(:,2)-sol(2)).^2 + (dxyz(:,3)-sol(3)).^2);
% simple forward model (could have used GPS2RNG again, more complicated forward models, etc.)
hst = hsr./vg;
rmse = norm(st(~isnan(st)) - hst(~isnan(hst)));

% plot rmse with a color bar and vn on x-axis, mag(xyzn) on y-axis
% this works, but we next want to complete inversion part to reduce
% rmse values and do this again
% may also want to change how the noise is being incorporated, I may have done this wrong
sz = 50;
scatter(vn,mag,sz,rmse,'filled')
a = colorbar;
a.Label.String = 'rmse';
a.FontSize = 11;
xlim([-10 10])
ylim([0 0.2])
grid on
longticks
xlabel('Velocity Noise [m/s]')
ylabel('XYZ Noise [m]')
hold on

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
