function gps2dis(files,protype)
% GPS2DIS(files,protype)
%
% Given Precise Point Position time series of four different units, computes
% their pairwise distances, and plots them
%
% INPUT:
% 
% files        cell with MAT-filename strings containing data structures
% protype      type of prd file ('ppp' or 'rtk')
%
% EXAMPLE:
%
% gps2his({'0001-05340.mat','0002-05340.mat','0003-05340.mat','0004-05340.mat'})
%
% Last modified by fjsimons-at-princeton.edu, 02/08/2022



