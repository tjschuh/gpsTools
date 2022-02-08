function [x0,y0,z0]=gps2gus(dxyz,depth)
% [x0,y0,z0]=GPS2GUS(dxyz,depth)
%
% Given a set of GPS coordinates and a depth, finds a likely drop point by
% determining the point of closest approach corresponding to that depth.
% This does not work great at all for complicated trajectories.
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
defval('depth',gebco(68+42/60,-(31+27/60)));

% Subsample for the search
sint=300;
x=dxyz(1:sint:end,1);
y=dxyz(1:sint:end,2);
z=dxyz(1:sint:end,3);

% All the subsampled GPS points in spherical coordinates
[AZ,EL,r]=cart2sph(x,y,z);
% All possible xyz points for every surface point except at depth
[x0,y0,z0]=sph2cart(AZ,EL,r-abs(depth));

% All the possible distances - could use XXPDIST
sr=sqrt((x(:)-x0(:)').^2+(y-y0(:)').^2+(z-z0(:)').^2);

% The minimum distance - there may be multiple, it's not that great
[d,k]=min(sr(:));
[i,j]=ind2sub(size(sr),k);

% The initial guess
x0=x0(j);
y0=y0(j);
z0=z0(j);

