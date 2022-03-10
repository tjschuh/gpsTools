function [x0,y0,z0]=gps2depth(dxyz,depth)
% [x0,y0,z0]=GPS2DEPTH(dxyz,depth)
%
% Given a point set of GPS coordinates and a depth, finds a point at depth
% below that point.
%
% INPUT:
%
% dxyz         The GPS points
% depth        The nominal depth
%
% OUTPUT:
%
% x0,y0,z0     A sensible subsurface location guess
%
% Last modified by fjsimons-at-alum.mit.edu, 02/07/2022

% Likely water depth from PrincetonSeafloorGeodesy-SURVEY3.pdf ORIGIN
%defval('depth',gebco(68+42/60,-(31+27/60)));

% The GPS points in spherical coordinates
[AZ,EL,r]=cart2sph(dxyz(:,1),dxyz(:,2),dxyz(:,3));
% All possible xyz points for every surface point except at depth
[x0,y0,z0]=sph2cart(AZ,EL,r-abs(depth));

