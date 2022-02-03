function gps2rng(files)
% GPS2RNG(files)
%
% Given Precise Point Position time series of four different units, compute
% a new average ship position time series from the four units and then
% compute the slant range to a fixed, imaginary beacon on the seafloor
%
% INPUT:
%
% files        cell with MAT-filename strings containing data structures
%
% EXAMPLE:
%
% gps2rng({'0001-05340.mat','0002-05340.mat','0003-05340.mat','0004-05340.mat'})
% gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'})
%
% Originally written by tschuh-at-princeton.edu, 11/24/2021
% Last modified by tschuh-at-princeton.edu 02/03/2022

% set beacon location on seafloor
% these come from file 05040 first xyz positions, but z-5km
% eventually will need to be a changing value until it reaches seafloor
% after that will not be known and will be deduced from timestamps sent back
dogx = 1.979e6;
dogy = -5.074e6;
dogz = 3.30385e6; %~5 km below surface

% combine all datasets into 1 with no plotting
[~,fname,~] = fileparts(files{1});
fname = sprintf('slantrange-%s',suf(fname,'-'));

d = mat2com(files,0);

% calculate slant range between ship and beacon for each second in meters
sr = sqrt((d.x(:) - dogx).^2 + (d.y(:) - dogy).^2 + (d.z(:) - dogz).^2);

% constant sound speed profile for now [m/s]
v = 1500;

% calculate slant time from slant range and v [s]
% currently not doing anything with this
st = sr./v;

% plot the distance [km] vs time
% would be nice to eventually add error bars
f=figure(1); clf
f.Position=[250 500 1100 600];

plot(d.t,sr*1e-3,'LineWidth',2)
xlim([d.t(1) d.t(end)])
buff = 0.01;
ylim([min(sr*1e-3)-buff*min(sr*1e-3) max(sr*1e-3)+buff*max(sr*1e-3)])
ylabel('Slant Range [km]')
title(sprintf('Raw Slant Range Measurements: %s-%s',datestr(d.t(1)),datestr(d.t(end))))
grid on
longticks([],2)

% save figure as pdf
figdisp(fname,[],'',2,[],'epstopdf')

% close figure
close