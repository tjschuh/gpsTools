function prd2msd(prdfile,xpos,ypos,zpos,freq)
% PRD2MSD(prdfile,xyzpos,freq)
%
% convert a .prd file to .mat and then to .mseed with the correct name
%
% INPUT:
%
% prdfile         *.prd output file produced via PRIDE PPP-AR
% xyzpos          approx position of GPS receiver from RINEX header
% freq            sampling interval of GPS data
%
% OUTPUT:
%
% .mat file      struct containing time series of East, North, & Up components of GPS station
% .mseed file    file containing time series of East, North, & Up components of GPS station
%
% EXAMPLES:
%
% prd2msd('kin_2011070_0842.prd',[-3904422.6794  3484842.8454  3633777.2141],1)
%
% Originally written by tschuh-at-princeton.edu, 04/29/2022
% Last modified by tschuh-at-princeton.edu, 06/23/2022

% save .prd as .mat file and as struct d
xyzpos=[xpos ypos zpos];
d = kin2pro(prdfile,xyzpos,0,0,0);
% need to separate d into multiple days
% non-array variables to exclude from the tabling procedure
drem={'timezone','ENUunit'};
% multiplier to signify level of precision (nm?)
multi = 1e6;
for i=1:3
	d.enu(:,i)=floor(multi.*d.enu(:,i));
end
% fill in missing data with a constant
d=retimes(d,drem,{'secondly','fillwithmissing'});
con = -99999999;
d.enu(isnan(d.enu))=con;
d.nsats(isnan(d.nsats))=con;
d.pdop(isnan(d.pdop))=con;
%d=retimes(d,drem,{'secondly','fillwithconstant','Constant',con});
% find total number of days by adding in last 15 minutes (900 seconds)
% and dividing by number of seconds in a day (86400 seconds)
sind = 86400;
numdays=(length(d.t)+899)/sind;

fifmin = con.*ones(899,1);
tims = d.t(end) + seconds(1:899);

% make mseed files for each day
for i=1:numdays
    % encoding format:
    % 11 = Steim-2 compression (data converted to int32)
    EF = 11;
    % data
    % last day is missing last 15 minutes so formula breaks down
    if i < numdays
        LXE = d.enu((sind*(i-1))+1:sind*i,1);
        LXN = d.enu((sind*(i-1))+1:sind*i,2);
        LXZ = d.enu((sind*(i-1))+1:sind*i,3);
        % accounting only for total num of sats
        DOP = d.pdop((sind*(i-1))+1:sind*i,1);
        SAT = d.nsats((sind*(i-1))+1:sind*i,1);
    else
        LXE = [d.enu((sind*(i-1))+1:end,1); fifmin];
        LXN = [d.enu((sind*(i-1))+1:end,2); fifmin];
        LXZ = [d.enu((sind*(i-1))+1:end,3); fifmin];
        DOP = [d.pdop((sind*(i-1))+1:end,1); fifmin];
        SAT = [d.nsats((sind*(i-1))+1:end,1); fifmin];
    end
    % data record length must by power of 2 >= 256
    % set to 1024 for now
    %sze = whos('LXE');
    %RLE = sze.bytes; RLE = pow2(ceil(log2(RLE)));
    RLE = 4096;
    %szn = whos('LXN');
    %RLN = szn.bytes; RLN = pow2(ceil(log2(RLN)));
    RLN = 4096;
    %szu = whos('LXZ');
    %RLU = szu.bytes; RLU = pow2(ceil(log2(RLU)));
    RLU = 4096;
    if RLE < 256 || RLN < 256 || RLU < 256
        error('record length < 256 bytes')
    end
    % finally create mseed file for each component of data

    if i < numdays
        mkmseed('GN.0842.00.LXE',LXE,d.t((sind*(i-1))+1:sind*i),freq,EF,RLE)
        mkmseed('GN.0842.00.LXN',LXN,d.t((sind*(i-1))+1:sind*i),freq,EF,RLN)
        mkmseed('GN.0842.00.LXZ',LXZ,d.t((sind*(i-1))+1:sind*i),freq,EF,RLU)
        mkmseed('GN.0842.00.DOP',DOP,d.t((sind*(i-1))+1:sind*i),freq,EF,RLU)
        mkmseed('GN.0842.00.SAT',SAT,d.t((sind*(i-1))+1:sind*i),freq,EF,RLU)
    else
        mkmseed('GN.0842.00.LXE',LXE,[d.t((sind*(i-1))+1:end); tims'],freq,EF,RLE)
        mkmseed('GN.0842.00.LXN',LXN,[d.t((sind*(i-1))+1:end); tims'],freq,EF,RLN)
        mkmseed('GN.0842.00.LXZ',LXZ,[d.t((sind*(i-1))+1:end); tims'],freq,EF,RLU)
	mkmseed('GN.0842.00.DOP',DOP,[d.t((sind*(i-1))+1:end); tims'],freq,EF,RLU)
        mkmseed('GN.0842.00.SAT',SAT,[d.t((sind*(i-1))+1:end); tims'],freq,EF,RLU)
    end
end
