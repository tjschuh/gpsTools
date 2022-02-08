function varargout = gps2inv(files,st,meth,v,guess)
% [sol,tru] = GPS2INV(files,st,meth,v,guess)
%
% Given Precise Point Position time series of four different units including
% slant times, consolidate them into 1 dataset of receivers with x,y,z,t
% to then do a least-squares inversion to calculate the predicted location
% of a beacon on the seafloor
%
% INPUT:
%
% recfiles      cell with MAT-filename strings containing data structures
% st            slant times/arrival times corresponding to the recfiles
% meth          method to convert recfiles into one dataset [default: 'ave']
% v             sound speed [m/s]
% guess         inital beacon location guess x,y,z,t 
%
% OUTPUT:
%
% sol           calculated [x y z t] beacon location [m]
% tru           actual [x y z 0] beacon location [m]
%
% EXAMPLE:
%
% [sol,tru] = gps2inv({'0001-05340.mat','0002-05340.mat','0003-05340.mat','0004-05340.mat'})
%
% Originally written by tschuh-at-princeton.edu, 02/07/2022
% Last modified by tschuh-at-princeton.edu, 02/08/2022
% Last modified by fjsimons-at-princeton.edu, 02/08/2022

% To Do:
% figure out how to speed up
% what's a good phi to use?
% what's another metric to determine when to stop iterating?
% add noise to slant times and receiver positions (covariance matrix)

% true beacon location (x,y,z,t)
tru = [1.979e6 -5.074e6 3.30385e6 0];

% initial guess
defval('guess',[2e6 -4.5e6 3e6 2])

defval('meth','ave')
% constant sound speed profile for now [m/s]
defval('v',1500)

% slant times
defval('st','st.mat')
load('st')

% set beacon location on seafloor
dogx=tru(1);
dogy=tru(2);
dogz=tru(3);

% Get extension from first filename
[~,fname,~]=fileparts(files{1});
filex=suf(fname,'-');

% combine all datasets into 1 large matrix with no plotting
[d,tmax]=mat2mod(files);

% Averaging method, or keeping them invidual, or taking only one,
switch meth
  case 'ave'
    % more or less one line version of mat2com.m
    dxyz=squeeze(nanmean(reshape(cat(1,d(:).xyz),size(d(1).xyz,1),length(d),3),2));
end

% combine station positions with time arrivals
% actual data observations
dxyzt = [dxyz st];

% calculate predicted travel times using initial guess
tpre = guess(4) + sqrt((guess(1) - dxyzt(:,1)).^2 + (guess(2) - dxyzt(:,2)).^2 + (guess(3) - dxyzt(:,3)).^2)/v;

% calculate sensitivity matrix G
% Gi1 = d(tpre)/dx, Gi2 = d(tpre)/dy,
% Gi3 = d(tpre)/dz, Gi4 = d(tpre)/dt
Gi1 = (v.*(guess(1) - dxyzt(:,1)))./sqrt((guess(1) - dxyzt(:,1)).^2 + (guess(2) - dxyzt(:,2)).^2 + (guess(3) - dxyzt(:,3)).^2);
Gi2 = (v.*(guess(2) - dxyzt(:,2)))./sqrt((guess(1) - dxyzt(:,1)).^2 + (guess(2) - dxyzt(:,2)).^2 + (guess(3) - dxyzt(:,3)).^2);
Gi3 = (v.*(guess(3) - dxyzt(:,3)))./sqrt((guess(1) - dxyzt(:,1)).^2 + (guess(2) - dxyzt(:,2)).^2 + (guess(3) - dxyzt(:,3)).^2);
Gi4 = ones(length(dxyzt),1);

G = [Gi1 Gi2 Gi3 Gi4];

% calculate difference between observed travel time and predicted travel time
deltad = dxyzt(:,4) - tpre;

% chi^2 = (deltad - Gdeltam)^2
% this is what we are trying to minimize
% so we differentiate and set equal to zero and solve for deltam
deltam = inv(transpose(G)*G)*transpose(G)*deltad;

% deltam contains our adjustments to our predicted solution (x,y,z,t)
% add deltam to our initial guess to get a new guess m
m = guess + deltam';

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

% optional output
varns={sol,tru};
varargout=varns(1:nargout);
