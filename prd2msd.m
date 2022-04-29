function prd2msd(prdfile,xyzpos,freq)
% PRD2MSD(prdfile,xyzpos,freq)
%
% convert a .prd file to .mat and then to mseed with the correct name
%
% INPUT:
%
% prdfile
% xyzpos
% freq
%
% OUTPUT:
%
% .mat file
% mseed file
%
% EXAMPLES:
%
% prd2msd('kin_2011070_0842.prd',[-3904422.6794  3484842.8454  3633777.2141],1)
%
% Originally written by tschuh-at-princeton.edu, 04/29/2022

% save .prd as .mat file and as struct d
d = kinpro(prdfile,xyzpos,0,0);

% encoding format:
% 11 = Steim-2 compression (data converted to int32)
EF = 11;
% multiplier to make data play well with Juypter Notebooks
multi = 10000;
% data
LXE = multi.*d.enu(:,1);
LXN = multi.*d.enu(:,2);
LXU = multi.*d.enu(:,3);
% data record length must by power of 2 >= 256
sze = whos('LXE');
RLE = sze.bytes;
RLE = pow2(ceil(log2(RLE)));
szn = whos('LXN');
RLN = szn.bytes;
RLN = pow2(ceil(log2(RLN)));
szu = whos('LXU');
RLU = szu.bytes;
RLU = pow2(ceil(log2(RLU)));
if RLE < 256 || RLN < 256 || RLU < 256
    error('record length < 256 bytes')
end
% finally create mseed file for each component of data
mkmseed('GN.0842.00.LXE',LXE,d.t,freq,EF,RLE)
mkmseed('GN.0842.00.LXN',LXN,d.t,freq,EF,RLN)
mkmseed('GN.0842.00.LXU',LXU,d.t,freq,EF,RLU)