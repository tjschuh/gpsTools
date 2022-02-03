function gps2std(unit1file,unit2file,unit3file,unit4file)
% GPS2STD(unit1file,unit2file,unit3file,unit4file)
%
% compute distances between 6 sets of receiver pairs
% calculate linear polyfit and residuals
% iteratively remove outliers from datasets to improve
% standard deviations and plot results
%
% INPUT:
%
% unit1file     mat file containing data collected by unit 1
% unit2file     mat file containing data collected by unit 2
% unit3file     mat file containing data collected by unit 3
% unit4file     mat file containing data collected by unit 4
%
% Originally written by tschuh-at-princeton.edu, 01/14/2022
% Last modified by tschuh-at-princeton.edu, 02/03/2022

% need to change structure to match gps2his.m

% use mat2mod to convert data to all be same time spans with no time gaps
[d1,d2,d3,d4] = mat2mod(unit1file,unit2file,unit3file,unit4file);
[~,fname,~] = fileparts(unit1file);

% compute distances between receivers
dist12 = sqrt((d1.xyz(:,1)-d2.xyz(:,1)).^2 + (d1.xyz(:,2)-d2.xyz(:,2)).^2 + (d1.xyz(:,3)-d2.xyz(:,3)).^2);
dist13 = sqrt((d1.xyz(:,1)-d3.xyz(:,1)).^2 + (d1.xyz(:,2)-d3.xyz(:,2)).^2 + (d1.xyz(:,3)-d3.xyz(:,3)).^2);
dist14 = sqrt((d1.xyz(:,1)-d4.xyz(:,1)).^2 + (d1.xyz(:,2)-d4.xyz(:,2)).^2 + (d1.xyz(:,3)-d4.xyz(:,3)).^2);
dist23 = sqrt((d2.xyz(:,1)-d3.xyz(:,1)).^2 + (d2.xyz(:,2)-d3.xyz(:,2)).^2 + (d2.xyz(:,3)-d3.xyz(:,3)).^2);
dist24 = sqrt((d2.xyz(:,1)-d4.xyz(:,1)).^2 + (d2.xyz(:,2)-d4.xyz(:,2)).^2 + (d2.xyz(:,3)-d4.xyz(:,3)).^2);
dist34 = sqrt((d3.xyz(:,1)-d4.xyz(:,1)).^2 + (d3.xyz(:,2)-d4.xyz(:,2)).^2 + (d3.xyz(:,3)-d4.xyz(:,3)).^2);

% find rows where nsats <= 4 and/or where pdop >= 15 or = 0
nthresh = 4; pthresh = 15;

% redefine pdop and nsats so they are easier to work with
p1 = d1.pdop; p2 = d2.pdop; p3 = d3.pdop; p4 = d4.pdop;
n1 = d1.nsats(:,1); n2 = d2.nsats(:,1); n3 = d3.nsats(:,1); n4 = d4.nsats(:,1);

% find good (g) and bad (b) data so we're only working with non-greyed data
% [g b] = dist
good12 = dist12; good13 = dist13; good14 = dist14; good23 = dist23; good24 = dist24; good34 = dist34; 
good12(p1>=pthresh | p1==0 | n1<=nthresh | p2>=pthresh | p2==0 | n2<=nthresh) = NaN;
good13(p1>=pthresh | p1==0 | n1<=nthresh | p3>=pthresh | p3==0 | n3<=nthresh) = NaN;
good14(p1>=pthresh | p1==0 | n1<=nthresh | p4>=pthresh | p4==0 | n4<=nthresh) = NaN;
good23(p2>=pthresh | p2==0 | n2<=nthresh | p3>=pthresh | p3==0 | n3<=nthresh) = NaN;
good24(p2>=pthresh | p2==0 | n2<=nthresh | p4>=pthresh | p4==0 | n4<=nthresh) = NaN;
good34(p3>=pthresh | p3==0 | n3<=nthresh | p4>=pthresh | p4==0 | n4<=nthresh) = NaN;

% remove any rows containing NaNs
d12 = rmNaNrows(good12);
d13 = rmNaNrows(good13);
d14 = rmNaNrows(good14);
d23 = rmNaNrows(good23);
d24 = rmNaNrows(good24);
d34 = rmNaNrows(good34);

% use d to find residuals (e)
p = polyfit([1:length(d12)]',d12,1); a12 = 1000*p(1); b12 = p(2);
x12 = (a12/1000).*[1:length(d12)]' + b12; e12 = 1000*(x12 - d12);
p = polyfit([1:length(d13)]',d13,1); a13 = 1000*p(1); b13 = p(2);
x13 = (a13/1000).*[1:length(d13)]' + b13; e13 = 1000*(x13 - d13);
p = polyfit([1:length(d14)]',d14,1); a14 = 1000*p(1); b14 = p(2);
x14 = (a14/1000).*[1:length(d14)]' + b14; e14 = 1000*(x14 - d14);
p = polyfit([1:length(d23)]',d23,1); a23 = 1000*p(1); b23 = p(2);
x23 = (a23/1000).*[1:length(d23)]' + b23; e23 = 1000*(x23 - d23);
p = polyfit([1:length(d24)]',d24,1); a24 = 1000*p(1); b24 = p(2);
x24 = (a24/1000).*[1:length(d24)]' + b24; e24 = 1000*(x24 - d24);
p = polyfit([1:length(d34)]',d34,1); a34 = 1000*p(1); b34 = p(2);
x34 = (a34/1000).*[1:length(d34)]' + b34; e34 = 1000*(x34 - d34);

% iteratively remove outliers to improve the "normality" of the data
pd12 = fitdist(e12,'Normal');
var12(1,1) = length(e12); var12(1,2) = pd12.sigma;
thresh = 2.5;
counter = 2;
while length(e12) > length(d1.t)/3
    e12 = rmoutliers(e12,'percentiles',[thresh 100-thresh]);
    thresh = thresh + 2.5;
    pd12 = fitdist(e12,'Normal');
    var12(counter,1) = length(e12); var12(counter,2) = pd12.sigma;
    counter = counter + 1;
end

pd13 = fitdist(e13,'Normal');
var13(1,1) = length(e13); var13(1,2) = pd13.sigma;
thresh = 2.5;
counter = 2;
while length(e13) > length(d1.t)/3
    e13 = rmoutliers(e13,'percentiles',[thresh 100-thresh]);
    thresh = thresh + 2.5;
    pd13 = fitdist(e13,'Normal');
    var13(counter,1) = length(e13); var13(counter,2) = pd13.sigma;
    counter = counter + 1;
end

pd14 = fitdist(e14,'Normal');
var14(1,1) = length(e14); var14(1,2) = pd14.sigma;
thresh = 2.5;
counter = 2;
while length(e14) > length(d1.t)/3
    e14 = rmoutliers(e14,'percentiles',[thresh 100-thresh]);
    thresh = thresh + 2.5;
    pd14 = fitdist(e14,'Normal');
    var14(counter,1) = length(e14); var14(counter,2) = pd14.sigma;
    counter = counter + 1;
end

pd23 = fitdist(e23,'Normal');
var23(1,1) = length(e23); var23(1,2) = pd23.sigma;
thresh = 2.5;
counter = 2;
while length(e23) > length(d1.t)/3
    e23 = rmoutliers(e23,'percentiles',[thresh 100-thresh]);
    thresh = thresh + 2.5;
    pd23 = fitdist(e23,'Normal');
    var23(counter,1) = length(e23); var23(counter,2) = pd23.sigma;
    counter = counter + 1;
end

pd24 = fitdist(e24,'Normal');
var24(1,1) = length(e24); var24(1,2) = pd24.sigma;
thresh = 2.5;
counter = 2;
while length(e24) > length(d1.t)/3
    e24 = rmoutliers(e24,'percentiles',[thresh 100-thresh]);
    thresh = thresh + 2.5;
    pd24 = fitdist(e24,'Normal');
    var24(counter,1) = length(e24); var24(counter,2) = pd24.sigma;
    counter = counter + 1;
end

pd34 = fitdist(e34,'Normal');
var34(1,1) = length(e34); var34(1,2) = pd34.sigma;
thresh = 2.5;
counter = 2;
while length(e34) > length(d1.t)/3
    e34 = rmoutliers(e34,'percentiles',[thresh 100-thresh]);
    thresh = thresh + 2.5;
    pd34 = fitdist(e34,'Normal');
    var34(counter,1) = length(e34); var34(counter,2) = pd34.sigma;
    counter = counter + 1;
end

maxstd = max([var12(1,2) var13(1,2) var14(1,2) var23(1,2) var24(1,2) var34(1,2)]);

% plotting
f=figure;
f.Position = [250 500 1100 800];

ah(1) = subplot(3,2,1);
plot(var12(:,1),var12(:,2))
cosmo(gca,'GPS Pair 1-2','# of Data Points','Std [mm]',e12,d1.t,maxstd)

ah(2) = subplot(3,2,2);
plot(var13(:,1),var13(:,2))
cosmo(gca,'GPS Pair 1-3','# of Data Points','Std [mm]',e13,d1.t,maxstd)

ah(3) = subplot(3,2,3);
plot(var14(:,1),var14(:,2))
cosmo(gca,'GPS Pair 1-4','# of Data Points','Std [mm]',e14,d1.t,maxstd)

ah(4) = subplot(3,2,4);
plot(var23(:,1),var23(:,2))
cosmo(gca,'GPS Pair 2-3','# of Data Points','Std [mm]',e23,d1.t,maxstd)

ah(5) = subplot(3,2,5);
plot(var24(:,1),var24(:,2))
cosmo(gca,'GPS Pair 2-4','# of Data Points','Std [mm]',e24,d1.t,maxstd)

ah(6) = subplot(3,2,6);
plot(var34(:,1),var34(:,2))
cosmo(gca,'GPS Pair 3-4','# of Data Points','Std [mm]',e34,d1.t,maxstd)

% finishing touches
tt=supertit(ah([1 2]),sprintf('Std vs # of Data Points (Ship Data from %s to %s)',datestr(d1.t(1)),datestr(d1.t(end))));
movev(tt,0.3)

%figdisp(sprintf('std-%s',fname),[],'',2,[],'epstopdf')

%close

% cosmetics
function cosmo(ax,titl,xlab,ylab,minlen,maxlen,maxstd)
set(ax,'XDir','reverse')
title(titl)
xlabel(xlab)
ylabel(ylab)
xlim([length(minlen) length(maxlen)])
ylim([0 maxstd+0.05*maxstd])
grid on
longticks