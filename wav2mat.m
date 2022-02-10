function wav2mat(firstfile,lastfile,rseg,plotver)
% WAV2MAT(firstfile,lastfile,rseg,plotver)
%
% Load in Acoustic data file, reshape rows,
% look for any data jumps, and plot the results
% WORKING VERSION FOR NOW
% NEEDS TO BE MERGED WITH CHPLOT
%
% INPUT:
%
% firstfile    the running number of the first file, e.g. 0, 1, 99, 362
% lastfile     the running number of the last file
% rseg         size of segment length in seconds [default: 2]
% plotver      switch (1 or 0) to turn on/off plotting feature [default: 1]
%
% OUTPUT:
%
% .mat file    output file saved as mat file to working directory
%
% EXAMPLE:
%
% wav2mat(697,697)
%
% Originally written by tschuh@princeton.edu, 02/10/2022

% By default, plotting feature is turned on
defval('plotver',1)
% Assume the record length in seconds
rlens = 60;
% Define a random segment length for inspection in seconds
defval('rseg',2)
% Assume the sampling rate
Fs = 400000;
% Moving average in seconds
mas = 1;
% Moving average left and right in units of samples
maslr = ([mas mas]/2)*Fs;

for file = firstfile:lastfile
    disp(sprintf('Working on file %3.4i (%3.4i of %3.4i)',file,file-firstfile+1,lastfile-firstfile+1))
    % Open file, turn it into matrix, and then
    % reshape it into a 4 "channel" matrix
    % Here is the filename
    fname = sprintf('file%d.data',file);
    sname=sprintf('file%d.mat',file);
    if exist(sname) ~= 2
      % Read the data
      fid = fopen(fname);
      FourChan = reshape(fread(fid,inf,'int16'),4,[]);
      fclose(fid);

      % Allocate the channels and identify jumps
      [FourChan,jumps] = challocate(FourChan);
      % Save it so next time you neither have to read it nor allocate it
      save(sname,'FourChan','jumps')
      fid = fopen(sname,'w');
      fwrite(fid,FourChan,'int16');
      fclose(fid);
    else
      disp(sprintf('loading file %s',sname))
      % fopen .mat file to get FourChan
      % not sure how to extract jumps
      % but not a priority right now
      fid = fopen(sname);
      FourChan = reshape(fread(fid,inf,'int16'),4,[]);
      fclose(fid);
    end
    
    % Define a random segment of length = rseg
    lowbound = randi(size(FourChan,2)-rseg*Fs-1);
    upbound = lowbound + rseg*Fs-1;
    sub = FourChan(:,lowbound:upbound);
    bot = Fs*ceil(lowbound/Fs);
    top = Fs*floor(upbound/Fs);
    tensec = [bot:Fs:top];

    for i = 1:size(FourChan,1)
        avg(i) = mean(FourChan(i,:));
    end

    % Open a figure
    f=figure;
    f.Position = [250 500 1050 550];
    
    % Colors
    color1 = [0 0 0];
    color2 = [0.4660 0.6740 0.1880];
    color3 = [0.6350 0.0780 0.1840];

    % Plot FourChan
    xtixl2 = 0:5:60;
    xtix2 = xtixl2*Fs;
    xlimit2 = [0 rlens*Fs];
    titl2 = {'PPS & Timestamp', 'Pre-Amped 400 kHz Acoustic Data', ...
             '12.5 kHz Bandpassed Acoustic Data', 'Low-Frequency Hydrophone'};

    ah(1)=subplot(3,2,1);
    plot(FourChan(1,:),'color',color1)
    ylim([4900 5700])
    yticks([5000 5600])
    cosmo(gca,titl2{1},xlimit2,xtix2,[])

    ah(2)=subplot(3,2,3);
    plot(FourChan(2,:),'color',color2)
    movev(ah(2),0.01)
    ylim([min(FourChan(2,:))-abs(.05*min(FourChan(2,:))) ...
          max(FourChan(2,:))+(.05*max(FourChan(2,:)))])
    yticks([min(FourChan(2,:)) round(avg(2)) max(FourChan(2,:))])
    cosmo(gca,titl2{2},xlimit2,xtix2,[])

    ah(3)=subplot(3,2,5);
    plot(FourChan(3,:),'color',color3)
    movev(ah(3),0.02)
    ylim([min(FourChan(3,:))-abs(.01*min(FourChan(3,:))) ...
          max(FourChan(3,:))+(.01*max(FourChan(3,:)))])
    yticks([min(FourChan(3,:)) round(avg(3)) max(FourChan(3,:))])
    cosmo(gca,titl2{3},xlimit2,xtix2,xtixl2)
    
    % Plot sub
    xtixl3 = tensec/Fs;
    xtix3 = (bot-lowbound):Fs:((10*Fs)-1-(upbound-top));
    xlimit3 = [0 rseg*Fs];
    titl3 = {rseg + " Second Segment"};

    ah(4)=subplot(3,2,2);
    plot(sub(1,:),'color',color1)
    ylim([4900 5700])
    yticks([5000 5600])
    cosmo(gca,titl3{1},xlimit3,xtix3,[])

    ah(5)=subplot(3,2,4);
    plot(sub(2,:),'color',color2)
    movev(ah(5),0.01)
    cosmo(gca,[],xlimit3,xtix3,[])
    
    ah(6)=subplot(3,2,6);
    plot(sub(3,:),'color',color3)
    movev(ah(6),0.02)
    cosmo(gca,[],xlimit3,xtix3,xtixl3)
    
    tt=supertit(ah([1 4]),sprintf('BIOS Unit2, Minute %d',file));
    movev(tt,0.2)

    % need to avoid annotation and use textbox...
    % not a priority right now
    %if exist('jumps','var') ~= 0
    %    a = annotation('textbox',[0.77 0.94 0 0],'String',['# of jumps = ' num2str(jumps)],'FitBoxToText','on');
    %    a.FontSize = 12;
    %else
    %    %jumps doesnt exist as a variable, so just ignore it
    %end

    % Save a PDF
    figdisp(sprintf('file%3.4i',file),[],[],2,[],'epstopdf')

    if firstfile == lastfile || file == lastfile
        %if working with 1 file, or working on the last file of a set, don't clf
    else
        clf
    end
end

% Cosmetics
function cosmo(ax,titl,xlimit,xtix,xtixl)
ax.XGrid = 'on';
ax.YGrid = 'off';
ax.GridColor = [0 0 0];
ax.TickLength = [0 0];
title(titl,'FontSize',14)
xlim(xlimit)
xticks(xtix)
xticklabels(xtixl)
%longticks([],4)