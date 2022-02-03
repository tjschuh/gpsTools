function pppvrtk(pppfile1,pppfile2,pppfile3,pppfile4,rtkfile1,rtkfile2,rtkfile3,rtkfile4)
% PPPVRTK(pppfile1,pppfile2,pppfile3,pppfile4,rtkfile1,rtkfile2,rtkfile3,rtkfile4)
%
% compare PPP and RTK by taking post-processed datasets from both methods
% across the same time spans of the the residuals of the distances between
% the 4 GPS receivers on-board R/V Atlantic Explorer and (1) plotting the 
% histograms back-to-back of the residuals, (2) plotting q-q plots showing 
% how close to a normal distribution each dataset is, (3) plotting the 
% standard deviation of each dataset as a function of the number of points 
% in the datasets, and (4) plotting the difference in standard deviation 
% as a function of number of points in the datasets to compapre PPP and RTK  
%
% INPUT:
%
% pppfile1     ppp mat file containing data collected by unit 1
% pppfile2     ppp mat file containing data collected by unit 2
% pppfile3     ppp mat file containing data collected by unit 3
% pppfile4     ppp mat file containing data collected by unit 4
% rtkfile1     rtk mat file containing data collected by unit 1
% rtkfile2     rtk mat file containing data collected by unit 2
% rtkfile3     rtk mat file containing data collected by unit 3
% rtkfile4     rtk mat file containing data collected by unit 4
%
% OUTPUT:
%
% subplot of histograms back-to-back of both datasets
% subplot of q-q plots for visual inspection of goodness of fit
% subplot of standard deviation vs number of points in both datasets
% subplot of PPP vs RTK using change in std from previous plot
%
% EXAMPLE:
%
% pppvrtk('0001-05340.mat','0002-05340.mat','0003-05340.mat','0004-05340.mat','0001-F089.mat','0002-F089.mat','0003-F089.mat','0004-F089.mat')
%
% Originally written by tschuh-at-princeton.edu, 01/17/2022
% Last modified by tschuh-at-princeton.edu, 02/03/2022

% need serious edits:
% split into multiple functions (at least 4)
% copy changes to gps2his.m

% make it so you dont need to be in directory with data
% this works, but it's commented out for now
% diro = '/home/tschuh/pppvsrtk/camp/';
% pppfile1 = strcat(diro,pppfile1);
% pppfile2 = strcat(diro,pppfile2);
% pppfile3 = strcat(diro,pppfile3);
% pppfile4 = strcat(diro,pppfile4);
% rtkfile1 = strcat(diro,rtkfile1);
% rtkfile2 = strcat(diro,rtkfile2);
% rtkfile3 = strcat(diro,rtkfile3);
% rtkfile4 = strcat(diro,rtkfile4);

% which leg is this hour from
trip = ['leg2'];

% (i=1) go through all 4 ppp datasets and compute pvar datasets
% which are the residuals from the GPS receiver distances, then
% (i=2) go through all 4 rtk datasets and compute rvar datasets
for i = 1:2
    if i == 1 % ppp
        % use mat2mod to convert ppp data to all be same time spans with no time gaps
        [d1,d2,d3,d4] = mat2mod(pppfile1,pppfile2,pppfile3,pppfile4);
        [~,pfname,~] = fileparts(pppfile1);
        ptime = d1.t;
    else % rtk
        % use mat2mod to convert rtk data to all be same time spans with no time gaps
        [d1,d2,d3,d4] = mat2mod(rtkfile1,rtkfile2,rtkfile3,rtkfile4);
        [~,rfname,~] = fileparts(rtkfile1);
        rtime = d1.t;
    end
        
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

    % use d to find residuals e
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

    if i == 1 % ppp
        % iteratively remove outliers to improve the "normality" of the data
        % make copies of e --> ee so we can still work with e
        pe12 = e12;
        pd12 = fitdist(e12,'Normal');
        pvar12(1,1) = length(e12); pvar12(1,2) = pd12.sigma;
        thresh = 2.5;
        counter = 2;
        while length(e12) > length(d1.t)/3
            e12 = rmoutliers(e12,'percentiles',[thresh 100-thresh]);
            thresh = thresh + 1;
            pd12 = fitdist(e12,'Normal');
            pvar12(counter,1) = length(e12); pvar12(counter,2) = pd12.sigma;
            counter = counter + 1;
        end

        pe13 = e13;
        pd13 = fitdist(e13,'Normal');
        pvar13(1,1) = length(e13); pvar13(1,2) = pd13.sigma;
        thresh = 2.5;
        counter = 2;
        while length(e13) > length(d1.t)/3
            e13 = rmoutliers(e13,'percentiles',[thresh 100-thresh]);
            thresh = thresh + 1;
            pd13 = fitdist(e13,'Normal');
            pvar13(counter,1) = length(e13); pvar13(counter,2) = pd13.sigma;
            counter = counter + 1;
        end

        pe14 = e14;
        pd14 = fitdist(e14,'Normal');
        pvar14(1,1) = length(e14); pvar14(1,2) = pd14.sigma;
        thresh = 2.5;
        counter = 2;
        while length(e14) > length(d1.t)/3
            e14 = rmoutliers(e14,'percentiles',[thresh 100-thresh]);
            thresh = thresh + 1;
            pd14 = fitdist(e14,'Normal');
            pvar14(counter,1) = length(e14); pvar14(counter,2) = pd14.sigma;
            counter = counter + 1;
        end

        pe23 = e23;
        pd23 = fitdist(e23,'Normal');
        pvar23(1,1) = length(e23); pvar23(1,2) = pd23.sigma;
        thresh = 2.5;
        counter = 2;
        while length(e23) > length(d1.t)/3
            e23 = rmoutliers(e23,'percentiles',[thresh 100-thresh]);
            thresh = thresh + 1;
            pd23 = fitdist(e23,'Normal');
            pvar23(counter,1) = length(e23); pvar23(counter,2) = pd23.sigma;
            counter = counter + 1;
        end

        pe24 = e24;
        pd24 = fitdist(e24,'Normal');
        pvar24(1,1) = length(e24); pvar24(1,2) = pd24.sigma;
        thresh = 2.5;
        counter = 2;
        while length(e24) > length(d1.t)/3
            e24 = rmoutliers(e24,'percentiles',[thresh 100-thresh]);
            thresh = thresh + 1;
            pd24 = fitdist(e24,'Normal');
            pvar24(counter,1) = length(e24); pvar24(counter,2) = pd24.sigma;
            counter = counter + 1;
        end

        pe34 = e34;
        pd34 = fitdist(e34,'Normal');
        pvar34(1,1) = length(e34); pvar34(1,2) = pd34.sigma;
        thresh = 2.5;
        counter = 2;
        while length(e34) > length(d1.t)/3
            e34 = rmoutliers(e34,'percentiles',[thresh 100-thresh]);
            thresh = thresh + 1;
            pd34 = fitdist(e34,'Normal');
            pvar34(counter,1) = length(e34); pvar34(counter,2) = pd34.sigma;
            counter = counter + 1;
        end
    else % rtk
        % iteratively remove outliers to improve the "normality" of the data
        re12 = e12;
        pd12 = fitdist(e12,'Normal');
        rvar12(1,1) = length(e12); rvar12(1,2) = pd12.sigma;
        thresh = 2.5;
        counter = 2;
        while length(e12) > length(d1.t)/3
            e12 = rmoutliers(e12,'percentiles',[thresh 100-thresh]);
            thresh = thresh + 1;
            pd12 = fitdist(e12,'Normal');
            rvar12(counter,1) = length(e12); rvar12(counter,2) = pd12.sigma;
            counter = counter + 1;
        end

        re13 = e13;
        pd13 = fitdist(e13,'Normal');
        rvar13(1,1) = length(e13); rvar13(1,2) = pd13.sigma;
        thresh = 2.5;
        counter = 2;
        while length(e13) > length(d1.t)/3
            e13 = rmoutliers(e13,'percentiles',[thresh 100-thresh]);
            thresh = thresh + 1;
            pd13 = fitdist(e13,'Normal');
            rvar13(counter,1) = length(e13); rvar13(counter,2) = pd13.sigma;
            counter = counter + 1;
        end

        re14 = e14;
        pd14 = fitdist(e14,'Normal');
        rvar14(1,1) = length(e14); rvar14(1,2) = pd14.sigma;
        thresh = 2.5;
        counter = 2;
        while length(e14) > length(d1.t)/3
            e14 = rmoutliers(e14,'percentiles',[thresh 100-thresh]);
            thresh = thresh + 1;
            pd14 = fitdist(e14,'Normal');
            rvar14(counter,1) = length(e14); rvar14(counter,2) = pd14.sigma;
            counter = counter + 1;
        end

        re23 = e23;
        pd23 = fitdist(e23,'Normal');
        rvar23(1,1) = length(e23); rvar23(1,2) = pd23.sigma;
        thresh = 2.5;
        counter = 2;
        while length(e23) > length(d1.t)/3
            e23 = rmoutliers(e23,'percentiles',[thresh 100-thresh]);
            thresh = thresh + 1;
            pd23 = fitdist(e23,'Normal');
            rvar23(counter,1) = length(e23); rvar23(counter,2) = pd23.sigma;
            counter = counter + 1;
        end

        re24 = e24;
        pd24 = fitdist(e24,'Normal');
        rvar24(1,1) = length(e24); rvar24(1,2) = pd24.sigma;
        thresh = 2.5;
        counter = 2;
        while length(e24) > length(d1.t)/3
            e24 = rmoutliers(e24,'percentiles',[thresh 100-thresh]);
            thresh = thresh + 1;
            pd24 = fitdist(e24,'Normal');
            rvar24(counter,1) = length(e24); rvar24(counter,2) = pd24.sigma;
            counter = counter + 1;
        end

        re34 = e34;
        pd34 = fitdist(e34,'Normal');
        rvar34(1,1) = length(e34); rvar34(1,2) = pd34.sigma;
        thresh = 2.5;
        counter = 2;
        while length(e34) > length(d1.t)/3
            e34 = rmoutliers(e34,'percentiles',[thresh 100-thresh]);
            thresh = thresh + 1;
            pd34 = fitdist(e34,'Normal');
            rvar34(counter,1) = length(e34); rvar34(counter,2) = pd34.sigma;
            counter = counter + 1;
        end
    end
end

maxstd = max([pvar12(1,2) pvar13(1,2) pvar14(1,2) pvar23(1,2) pvar24(1,2) pvar34(1,2) ...
              rvar12(1,2) rvar13(1,2) rvar14(1,2) rvar23(1,2) rvar24(1,2) rvar34(1,2)]);

% plotting
% plotting histograms first
f=figure;
f.Position = [250 500 1100 800];

ah1(1) = subplot(3,2,1);
% remove outliers to get better results
try
    pee12 = rmoutliers(pe12,'gesd');
    ree12 = rmoutliers(re12,'gesd');
catch
    pee12 = rmoutliers(pe12,'percentiles',[5 95]);
    ree12 = rmoutliers(re12,'percentiles',[5 95]);
end
% calculate nbins using Freedman-Diaconis Method
pnbins12 = round((max(pee12) - min(pee12))/(2*iqr(pee12)*(length(pee12))^(-1/3)));
rnbins12 = round((max(ree12) - min(ree12))/(2*iqr(ree12)*(length(ree12))^(-1/3)));
% force both histograms to have same number of bins/std (bpstd)
pbpstd12 = pnbins12/std(pee12);
rbpstd12 = rnbins12/std(ree12);
% use mean of bpstd to find new nbins
pnbins12 = round(std(pee12)*(pbpstd12+rbpstd12)/2);
rnbins12 = round(std(ree12)*(pbpstd12+rbpstd12)/2);
% calculate goodness of fit compared to normal distribution (gof)
[~,~,pstats12] = chi2gof(pee12,'NBins',pnbins12);
[~,~,rstats12] = chi2gof(ree12,'NBins',rnbins12);
pgof12 = pstats12.chi2stat/pstats12.df;
rgof12 = rstats12.chi2stat/rstats12.df;
% make histograms and overlaid curves using histcounts, bar, and plot
[pyvals12,pedges12] = histcounts(pee12,pnbins12);
pbarc12 = 0.5*(pedges12(1:end-1) + pedges12(2:end));
pbar12 = bar(pbarc12,pyvals12,'BarWidth',1);
hold on
[ryvals12,redges12] = histcounts(ree12,rnbins12);
rbarc12 = 0.5*(redges12(1:end-1) + redges12(2:end));
rbar12 = bar(rbarc12,ryvals12,'BarWidth',1);
[pline12,rline12] = cosmo1(gca,ah1(1),sprintf('PPP = %i/%i, GPS Pair 1-2, RTK = %i/%i',length(pee12),length(ptime),length(ree12),length(rtime)),'Residuals [mm]','Counts',pee12,pbar12,pgof12,ree12,rbar12,rgof12);

ah1(2) = subplot(3,2,2);
% remove outliers to get better results
try
    pee13 = rmoutliers(pe13,'gesd');
    ree13 = rmoutliers(re13,'gesd');
catch
    pee13 = rmoutliers(pe13,'percentiles',[5 95]);
    ree13 = rmoutliers(re13,'percentiles',[5 95]);
end
pnbins13 = round((max(pee13) - min(pee13))/(2*iqr(pee13)*(length(pee13))^(-1/3)));
rnbins13 = round((max(ree13) - min(ree13))/(2*iqr(ree13)*(length(ree13))^(-1/3)));
pbpstd13 = pnbins13/std(pee13);
rbpstd13 = rnbins13/std(ree13);
pnbins13 = round(std(pee13)*(pbpstd13+rbpstd13)/2);
rnbins13 = round(std(ree13)*(pbpstd13+rbpstd13)/2);
[~,~,pstats13] = chi2gof(pee13,'NBins',pnbins13);
[~,~,rstats13] = chi2gof(ree13,'NBins',rnbins13);
pgof13 = pstats13.chi2stat/pstats13.df;
rgof13 = rstats13.chi2stat/rstats13.df;
[pyvals13,pedges13] = histcounts(pee13,pnbins13);
pbarc13 = 0.5*(pedges13(1:end-1) + pedges13(2:end));
pbar13 = bar(pbarc13,pyvals13,'BarWidth',1);
hold on
[ryvals13,redges13] = histcounts(ree13,rnbins13);
rbarc13 = 0.5*(redges13(1:end-1) + redges13(2:end));
rbar13 = bar(rbarc13,ryvals13,'BarWidth',1);
[pline13,rline13] = cosmo1(gca,ah1(2),sprintf('PPP = %i/%i, GPS Pair 1-3, RTK = %i/%i',length(pee13),length(ptime),length(ree13),length(rtime)),'Residuals [mm]','Counts',pee13,pbar13,pgof13,ree13,rbar13,rgof13);

ah1(3) = subplot(3,2,3);
% remove outliers to get better results
try
    pee14 = rmoutliers(pe14,'gesd');
    ree14 = rmoutliers(re14,'gesd');
catch
    pee14 = rmoutliers(pe14,'percentiles',[5 95]);
    ree14 = rmoutliers(re14,'percentiles',[5 95]);
end
pnbins14 = round((max(pee14) - min(pee14))/(2*iqr(pee14)*(length(pee14))^(-1/3)));
rnbins14 = round((max(ree14) - min(ree14))/(2*iqr(ree14)*(length(ree14))^(-1/3)));
pbpstd14 = pnbins14/std(pee14);
rbpstd14 = rnbins14/std(ree14);
pnbins14 = round(std(pee14)*(pbpstd14+rbpstd14)/2);
rnbins14 = round(std(ree14)*(pbpstd14+rbpstd14)/2);
[~,~,pstats14] = chi2gof(pee14,'NBins',pnbins14);
[~,~,rstats14] = chi2gof(ree14,'NBins',rnbins14);
pgof14 = pstats14.chi2stat/pstats14.df;
rgof14 = rstats14.chi2stat/rstats14.df;
[pyvals14,pedges14] = histcounts(pee14,pnbins14);
pbarc14 = 0.5*(pedges14(1:end-1) + pedges14(2:end));
pbar14 = bar(pbarc14,pyvals14,'BarWidth',1);
hold on
[ryvals14,redges14] = histcounts(ree14,rnbins14);
rbarc14 = 0.5*(redges14(1:end-1) + redges14(2:end));
rbar14 = bar(rbarc14,ryvals14,'BarWidth',1);
[pline14,rline14] = cosmo1(gca,ah1(3),sprintf('PPP = %i/%i, GPS Pair 1-4, RTK = %i/%i',length(pee14),length(ptime),length(ree14),length(rtime)),'Residuals [mm]','Counts',pee14,pbar14,pgof14,ree14,rbar14,rgof14);

ah1(4) = subplot(3,2,4);
% remove outliers to get better results
try
    pee23 = rmoutliers(pe23,'gesd');
    ree23 = rmoutliers(re23,'gesd');
catch
    pee23 = rmoutliers(pe23,'percentiles',[5 95]);
    ree23 = rmoutliers(re23,'percentiles',[5 95]);
end
pnbins23 = round((max(pee23) - min(pee23))/(2*iqr(pee23)*(length(pee23))^(-1/3)));
rnbins23 = round((max(ree23) - min(ree23))/(2*iqr(ree23)*(length(ree23))^(-1/3)));
pbpstd23 = pnbins23/std(pee23);
rbpstd23 = rnbins23/std(ree23);
pnbins23 = round(std(pee23)*(pbpstd23+rbpstd23)/2);
rnbins23 = round(std(ree23)*(pbpstd23+rbpstd23)/2);
[~,~,pstats23] = chi2gof(pee23,'NBins',pnbins23);
[~,~,rstats23] = chi2gof(ree23,'NBins',rnbins23);
pgof23 = pstats23.chi2stat/pstats23.df;
rgof23 = rstats23.chi2stat/rstats23.df;
[pyvals23,pedges23] = histcounts(pee23,pnbins23);
pbarc23 = 0.5*(pedges23(1:end-1) + pedges23(2:end));
pbar23 = bar(pbarc23,pyvals23,'BarWidth',1);
hold on
[ryvals23,redges23] = histcounts(ree23,rnbins23);
rbarc23 = 0.5*(redges23(1:end-1) + redges23(2:end));
rbar23 = bar(rbarc23,ryvals23,'BarWidth',1);
[pline23,rline23] = cosmo1(gca,ah1(4),sprintf('PPP = %i/%i, GPS Pair 2-3, RTK = %i/%i',length(pee23),length(ptime),length(ree23),length(rtime)),'Residuals [mm]','Counts',pee23,pbar23,pgof23,ree23,rbar23,rgof23);

ah1(5) = subplot(3,2,5);
% remove outliers to get better results
try
    pee24 = rmoutliers(pe24,'gesd');
    ree24 = rmoutliers(re24,'gesd');
catch
    pee24 = rmoutliers(pe24,'percentiles',[5 95]);
    ree24 = rmoutliers(re24,'percentiles',[5 95]);
end
pnbins24 = round((max(pee24) - min(pee24))/(2*iqr(pee24)*(length(pee24))^(-1/3)));
rnbins24 = round((max(ree24) - min(ree24))/(2*iqr(ree24)*(length(ree24))^(-1/3)));
pbpstd24 = pnbins24/std(pee24);
rbpstd24 = rnbins24/std(ree24);
pnbins24 = round(std(pee24)*(pbpstd24+rbpstd24)/2);
rnbins24 = round(std(ree24)*(pbpstd24+rbpstd24)/2);
[~,~,pstats24] = chi2gof(pee24,'NBins',pnbins24);
[~,~,rstats24] = chi2gof(ree24,'NBins',rnbins24);
pgof24 = pstats24.chi2stat/pstats24.df;
rgof24 = rstats24.chi2stat/rstats24.df;
[pyvals24,pedges24] = histcounts(pee24,pnbins24);
pbarc24 = 0.5*(pedges24(1:end-1) + pedges24(2:end));
pbar24 = bar(pbarc24,pyvals24,'BarWidth',1);
hold on
[ryvals24,redges24] = histcounts(ree24,rnbins24);
rbarc24 = 0.5*(redges24(1:end-1) + redges24(2:end));
rbar24 = bar(rbarc24,ryvals24,'BarWidth',1);
[pline24,rline24] = cosmo1(gca,ah1(5),sprintf('PPP = %i/%i, GPS Pair 2-4, RTK = %i/%i',length(pee24),length(ptime),length(ree24),length(rtime)),'Residuals [mm]','Counts',pee24,pbar24,pgof24,ree24,rbar24,rgof24);

ah1(6) = subplot(3,2,6);
% remove outliers to get better results
try
    pee34 = rmoutliers(pe34,'gesd');
    ree34 = rmoutliers(re34,'gesd');
catch
    pee34 = rmoutliers(pe34,'percentiles',[5 95]);
    ree34 = rmoutliers(re34,'percentiles',[5 95]);
end
pnbins34 = round((max(pee34) - min(pee34))/(2*iqr(pee34)*(length(pee34))^(-1/3)));
rnbins34 = round((max(ree34) - min(ree34))/(2*iqr(ree34)*(length(ree34))^(-1/3)));
pbpstd34 = pnbins34/std(pee34);
rbpstd34 = rnbins34/std(ree34);
pnbins34 = round(std(pee34)*(pbpstd34+rbpstd34)/2);
rnbins34 = round(std(ree34)*(pbpstd34+rbpstd34)/2);
[~,~,pstats34] = chi2gof(pee34,'NBins',pnbins34);
[~,~,rstats34] = chi2gof(ree34,'NBins',rnbins34);
pgof34 = pstats34.chi2stat/pstats34.df;
rgof34 = rstats34.chi2stat/rstats34.df;
[pyvals34,pedges34] = histcounts(pee34,pnbins34);
pbarc34 = 0.5*(pedges34(1:end-1) + pedges34(2:end));
pbar34 = bar(pbarc34,pyvals34,'BarWidth',1);
hold on
[ryvals34,redges34] = histcounts(ree34,rnbins34);
rbarc34 = 0.5*(redges34(1:end-1) + redges34(2:end));
rbar34 = bar(rbarc34,ryvals34,'BarWidth',1);
[pline34,rline34] = cosmo1(gca,ah1(6),sprintf('PPP = %i/%i, GPS Pair 3-4, RTK = %i/%i',length(pee34),length(ptime),length(ree34),length(rtime)),'Residuals [mm]','Counts',pee34,pbar34,pgof34,ree34,rbar34,rgof34);

% finishing touches
tt=supertit(ah1([1 2]),sprintf('Demeaned Residuals of Ship Data from %s to %s',datestr(rtime(1)),datestr(rtime(end))));
movev(tt,0.3)
b = annotation('textbox',[0.485 0.075 0 0],'String',trip,'FitBoxToText','on');
b.FontSize = 12;

figdisp(sprintf('histo-%s-%s',pfname,rfname),[],'',2,[],'epstopdf')

close

% plot the qq plots second
g=figure;
g.Position = [250 500 1100 800];

ah2(1) = subplot(3,2,1);
pqq12 = qqplot(pee12);
hold on
rqq12 = qqplot(ree12);
cosmo2('GPS Pair 1-2',pqq12,rqq12)

ah2(2) = subplot(3,2,2);
pqq13 = qqplot(pee13);
hold on
rqq13 = qqplot(ree13);
cosmo2('GPS Pair 1-3',pqq13,rqq13)

ah2(3) = subplot(3,2,3);
pqq14 = qqplot(pee14);
hold on
rqq14 = qqplot(ree14);
cosmo2('GPS Pair 1-4',pqq14,rqq14)

ah2(4) = subplot(3,2,4);
pqq23 = qqplot(pee23);
hold on
rqq23 = qqplot(ree23);
cosmo2('GPS Pair 2-3',pqq23,rqq23)

ah2(5) = subplot(3,2,5);
pqq24 = qqplot(pee24);
hold on
rqq24 = qqplot(ree12);
cosmo2('GPS Pair 2-4',pqq24,rqq24)

ah2(6) = subplot(3,2,6);
pqq34 = qqplot(pee34);
hold on
rqq34 = qqplot(ree34);
cosmo2('GPS Pair 3-4',pqq34,rqq34)

% finishing touches
tt=supertit(ah2([1 2]),sprintf('QQ Plots of Residuals vs Standard Normals (%s to %s)',datestr(d1.t(1)),datestr(d1.t(end))));
movev(tt,0.3)
b = annotation('textbox',[0.485 0.075 0 0],'String',trip,'FitBoxToText','on');
b.FontSize = 12;

figdisp(sprintf('qq-%s-%s',pfname,rfname),[],'',2,[],'epstopdf')

close

% plot the std curves third
h=figure;
h.Position = [250 500 1100 800];

lblue = [0.4 0.6667 0.8431];
lgreen = [0.466 0.674 0.188];

ah3(1) = subplot(3,2,1);
plot(pvar12(:,1),pvar12(:,2),'Color',lblue,'LineWidth',2)
hold on
plot(rvar12(:,1),rvar12(:,2),'Color',lgreen,'LineWidth',2)
cosmo3(gca,'GPS Pair 1-2','# of Data Points','Std [mm]',e12,d1.t,maxstd)

ah3(2) = subplot(3,2,2);
plot(pvar13(:,1),pvar13(:,2),'Color',lblue,'LineWidth',2)
hold on
plot(rvar13(:,1),rvar13(:,2),'Color',lgreen,'LineWidth',2)
cosmo3(gca,'GPS Pair 1-3','# of Data Points','Std [mm]',e13,d1.t,maxstd)

ah3(3) = subplot(3,2,3);
plot(pvar14(:,1),pvar14(:,2),'Color',lblue,'LineWidth',2)
hold on
plot(rvar14(:,1),rvar14(:,2),'Color',lgreen,'LineWidth',2)
cosmo3(gca,'GPS Pair 1-4','# of Data Points','Std [mm]',e14,d1.t,maxstd)

ah3(4) = subplot(3,2,4);
plot(pvar23(:,1),pvar23(:,2),'Color',lblue,'LineWidth',2)
hold on
plot(rvar23(:,1),rvar23(:,2),'Color',lgreen,'LineWidth',2)
cosmo3(gca,'GPS Pair 2-3','# of Data Points','Std [mm]',e23,d1.t,maxstd)

ah3(5) = subplot(3,2,5);
plot(pvar24(:,1),pvar24(:,2),'Color',lblue,'LineWidth',2)
hold on
plot(rvar24(:,1),rvar24(:,2),'Color',lgreen,'LineWidth',2)
cosmo3(gca,'GPS Pair 2-4','# of Data Points','Std [mm]',e24,d1.t,maxstd)

ah3(6) = subplot(3,2,6);
plot(pvar34(:,1),pvar34(:,2),'Color',lblue,'LineWidth',2)
hold on
plot(rvar34(:,1),rvar34(:,2),'Color',lgreen,'LineWidth',2)
cosmo3(gca,'GPS Pair 3-4','# of Data Points','Std [mm]',e34,d1.t,maxstd)

% finishing touches
tt=supertit(ah3([1 2]),sprintf('Std vs # of Data Points (Ship Data from %s to %s)',datestr(rtime(1)),datestr(rtime(end))));
movev(tt,0.3)
b = annotation('textbox',[0.485 0.075 0 0],'String',trip,'FitBoxToText','on');
b.FontSize = 12;

figdisp(sprintf('std-%s-%s',pfname,rfname),[],'',2,[],'epstopdf')

close

% plot the ppp vs rtk 2 graph fourth
% + --> PPP better, - --> RTK better
% this works sometimes, it depends on lengths of pvar and rvar
% take difference between rvar and pvar values
try
diff12 = rvar12(:,2)-pvar12(:,2);
diff13 = rvar13(:,2)-pvar13(:,2);
diff14 = rvar14(:,2)-pvar14(:,2);
diff23 = rvar23(:,2)-pvar23(:,2);
diff24 = rvar24(:,2)-pvar24(:,2);
diff34 = rvar34(:,2)-pvar34(:,2);

% find global max and min
maxdiff = max([diff12 diff13 diff14 diff23 diff24 diff34],[],'all');
mindiff = min([diff12 diff13 diff14 diff23 diff24 diff34],[],'all');

k=figure;
k.Position = [250 500 1100 800];

ah4(1) = subplot(3,2,1);
plot((pvar12(:,1)+rvar12(:,1))./2,diff12,'LineWidth',2)
cosmo4(gca,'GPS Pair 1-2, Positive --> PPP, Negative --> RTK','# of Data Points','Change in Std',e12,d1.t,mindiff,maxdiff)

ah4(2) = subplot(3,2,2);
plot((pvar13(:,1)+rvar13(:,1))./2,diff13,'LineWidth',2)
cosmo4(gca,'GPS Pair 1-3, Positive --> PPP, Negative --> RTK','# of Data Points','Change in Std',e13,d1.t,mindiff,maxdiff)

ah4(3) = subplot(3,2,3);
plot((pvar14(:,1)+rvar14(:,1))./2,diff14,'LineWidth',2)
cosmo4(gca,'GPS Pair 1-4, Positive --> PPP, Negative --> RTK','# of Data Points','Change in Std',e14,d1.t,mindiff,maxdiff)

ah4(4) = subplot(3,2,4);
plot((pvar23(:,1)+rvar23(:,1))./2,diff23,'LineWidth',2)
cosmo4(gca,'GPS Pair 2-3, Positive --> PPP, Negative --> RTK','# of Data Points','Change in Std',e23,d1.t,mindiff,maxdiff)

ah4(5) = subplot(3,2,5);
plot((pvar24(:,1)+rvar24(:,1))./2,diff24,'LineWidth',2)
cosmo4(gca,'GPS Pair 2-4, Positive --> PPP, Negative --> RTK','# of Data Points','Change in Std',e24,d1.t,mindiff,maxdiff)

ah4(6) = subplot(3,2,6);
plot((pvar34(:,1)+rvar34(:,1))./2,diff34,'LineWidth',2)
cosmo4(gca,'GPS Pair 3-4, Positive --> PPP, Negative --> RTK','# of Data Points','Change in Std',e34,d1.t,mindiff,maxdiff)

% finishing touches
tt=supertit(ah4([1 2]),sprintf('PPP vs RTK (Ship Data from %s to %s)',datestr(rtime(1)),datestr(rtime(end))));
movev(tt,0.3)
b = annotation('textbox',[0.485 0.075 0 0],'String',trip,'FitBoxToText','on');
b.FontSize = 12;

figdisp(sprintf('comp-%s-%s',pfname,rfname),[],'',2,[],'epstopdf')

close
catch
end

% cosmetics for histogram plots
function [pline,rline] = cosmo1(ax,ah,titl,xlab,ylab,pdata,pbar,pgof,rdata,rbar,rgof)
% plot rtk histogram below x-axis
rbar.YData = -1*rbar.YData;
ax.XGrid = 'on';
ax.YGrid = 'off';
ax.GridColor = [0 0 0];
ax.TickLength = [0 0];
title(titl)
xlabel(xlab)
ylabel(ylab)
if abs(max(pbar.YData)) >= abs(min(rbar.YData))
    ylim([-max(pbar.YData)-0.1*max(pbar.YData) max(pbar.YData)+0.1*max(pbar.YData)])
elseif abs(max(pbar.YData)) < abs(min(rbar.YData))
    ylim([min(rbar.YData)+0.1*min(rbar.YData) -min(rbar.YData)-0.1*min(rbar.YData)])
end
if std(pdata)>=std(rdata)
    xlim([round(-3*std(pdata),2) round(3*std(pdata),2)])
    xticks([round(-3*std(pdata),2) round(-2*std(pdata),2) round(-std(pdata),2) 0 round(std(pdata),2) round(2*std(pdata),2) round(3*std(pdata),2)])
    xticklabels({round(-3*std(pdata),0),round(-2*std(pdata),0),round(-std(pdata),0),0,round(std(pdata),0),round(2*std(pdata),0),round(3*std(pdata),0)})
    text(1.9*std(pdata),ah.YLim(2)-0.45*ah.YLim(2),sprintf('std = %.0f\nmed = %.0f\nmean = %.0f\ngof = %.0f',std(pdata),median(pdata),mean(pdata),pgof),'FontSize',9)
    text(1.9*std(pdata),ah.YLim(1)-0.45*ah.YLim(1),sprintf('std = %.0f\nmed = %.0f\nmean = %.0f\ngof = %.0f',std(rdata),median(rdata),mean(rdata),rgof),'FontSize',9)
    ppct = (length(pdata(pdata<=round(3*std(pdata)) & pdata>=round(-3*std(pdata))))/length(pdata))*100;
    text(-2.9*std(pdata),ah.YLim(2)-0.4*ah.YLim(2),sprintf('%05.2f%%\nmin = %.0f\nmax = %.0f',ppct,min(pdata),max(pdata)),'FontSize',9)
    rpct = (length(rdata(rdata<=round(3*std(rdata)) & rdata>=round(-3*std(rdata))))/length(rdata))*100;
    text(-2.9*std(pdata),ah.YLim(1)-0.4*ah.YLim(1),sprintf('%05.2f%%\nmin = %.0f\nmax = %.0f',rpct,min(rdata),max(rdata)),'FontSize',9)
elseif std(pdata)<std(rdata)
    xlim([round(-3*std(rdata),2) round(3*std(rdata),2)])
    xticks([round(-3*std(rdata),2) round(-2*std(rdata),2) round(-std(rdata),2) 0 round(std(rdata),2) round(2*std(rdata),2) round(3*std(rdata),2)])
    xticklabels({round(-3*std(rdata),0),round(-2*std(rdata),0),round(-std(rdata),0),0,round(std(rdata),0),round(2*std(rdata),0),round(3*std(rdata),0)})
    text(1.9*std(rdata),ah.YLim(2)-0.45*ah.YLim(2),sprintf('std = %.0f\nmed = %.0f\nmean = %.0f\ngof = %.0f',std(pdata),median(pdata),mean(pdata),pgof),'FontSize',9)
    text(1.9*std(rdata),ah.YLim(1)-0.45*ah.YLim(1),sprintf('std = %.0f\nmed = %.0f\nmean = %.0f\ngof = %.0f',std(rdata),median(rdata),mean(rdata),rgof),'FontSize',9)
    ppct = (length(pdata(pdata<=round(3*std(pdata)) & pdata>=round(-3*std(pdata))))/length(pdata))*100;
    text(-2.9*std(rdata),ah.YLim(2)-0.4*ah.YLim(2),sprintf('%05.2f%%\nmin = %.0f\nmax = %.0f',ppct,min(pdata),max(pdata)),'FontSize',9)
    rpct = (length(rdata(rdata<=round(3*std(rdata)) & rdata>=round(-3*std(rdata))))/length(rdata))*100;
    text(-2.9*std(rdata),ah.YLim(1)-0.4*ah.YLim(1),sprintf('%05.2f%%\nmin = %.0f\nmax = %.0f',rpct,min(rdata),max(rdata)),'FontSize',9)
end
longticks([],2)
pbar.FaceColor = [0.4 0.6667 0.8431]; %light blue bars
rbar.FaceColor = [0.466 0.674 0.188]; %lime green bars
% plot vertical lines at respective medians
line([median(pdata) median(pdata)],[0 max(pbar.YData)],'Color','k','LineStyle','--','LineWidth',1.5);
line([median(rdata) median(rdata)],[0 min(rbar.YData)],'Color','k','LineStyle','--','LineWidth',1.5);
% plot horizontal line at 0
yline(0,'k','LineWidth',1.5);
% plot pd curves on top of histograms
ppd = fitdist(pdata,'Normal');
pxvals = [-3*std(pdata):3*std(pdata)];
pyvals = pdf(ppd,pxvals);
parea = sum(pbar.YData)*diff(pbar.XData(1:2));
pline = plot(pxvals,pyvals*parea,'r','LineWidth',2);
% if gof > 4 make red curve dotted b/c data is not normal
if pgof > 4
    pline.LineStyle = '--';
end
rpd = fitdist(rdata,'Normal');
rxvals = [-3*std(rdata):3*std(rdata)];
ryvals = pdf(rpd,rxvals);
rarea = sum(rbar.YData)*diff(rbar.XData(1:2));
rline = plot(rxvals,ryvals*rarea,'r','LineWidth',2);
% if gof > 4 make red curve dotted b/c data is not normal
if rgof > 4
    rline.LineStyle = '--';
end

% cosmetics for qq plots
function cosmo2(titl,pqq,rqq)
grid on
longticks([],2)
title(titl)
ylabel("Residual Quantiles")
pqq(1).Marker = '.';
pqq(1).MarkerSize = 8;
pqq(1).MarkerEdgeColor = [0 0 1];
pqq(2).LineWidth = 2;
pqq(2).LineStyle = '--';
pqq(2).Color = [0 0.447 0.741];
pqq(3).Color = [0 0.447 0.741];
pqq(3).LineWidth = 2;
pqq(3).LineStyle = '--';
rqq(1).Marker = '.';
rqq(1).MarkerSize = 8;
rqq(1).MarkerEdgeColor = [0 0.5 0];
rqq(2).LineWidth = 2;
rqq(2).LineStyle = '--';
rqq(2).Color = [0.466 0.647 0.188];
rqq(3).Color = [0.466 0.647 0.188];
rqq(3).LineWidth = 2;
rqq(3).LineStyle = '--';
legend([pqq(1) rqq(1)],{'PPP','RTK'},'Location','northwest')

% cosmetics for std plots
function cosmo3(ax,titl,xlab,ylab,minlen,maxlen,maxstd)
set(ax,'XDir','reverse')
title(titl)
xlabel(xlab)
ylabel(ylab)
xlim([length(minlen) length(maxlen)])
ylim([0 maxstd+0.05*maxstd])
legend('PPP','RTK')
grid on
longticks([],2)

% cosmetics for ppp vs rtk plots
function cosmo4(ax,titl,xlab,ylab,minlen,maxlen,mindiff,maxdiff)
set(ax,'XDir','reverse')
title(titl)
xlabel(xlab)
ylabel(ylab)
xlim([length(minlen) length(maxlen)])
if mindiff >= 0
    ylim([0 maxdiff+0.1*maxdiff])
elseif maxdiff <= 0
    ylim([mindiff-0.1*mindiff 0])
else
    ylim([mindiff+0.1*mindiff maxdiff+0.1*maxdiff])
end
grid on
longticks([],2)