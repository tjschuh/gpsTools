function dog2rng(fname)
% DOG2RNG(fname)
%
% Works on file(s) made by DOG2MAT and turns it into a range measurement
%
% INPUT:
%
% fname         file(s) made by DOG2MAT with time tags (up to 4)
%
% EXAMPLE:
%
% dog2rng({'DOG1-camp.mat','DOG2-camp.mat','DOG3-camp.mat','DOG4-camp.mat'})
%
% Originally written by fjsimons-at-princeton.edu, 02/07/2022
% Last modified by tschuh-at-princeton.edu, 07/20/2022

% make this much cleaner and lose duplicate code
% play with cosmetics of figure

% DOG 1-4 approx positions
xyz=[1.977967 -5.073198 3.3101016; 1.9765409 -5.074798706 3.308558817; ...
     1.980083296 -5.073417583 3.308558817; 1.978732384 -5.075186123 3.30666684].*1e6;
ofx1=[5823; 1; 1; 2250]; ofx2=[4800; 12500; 12000; 1]; ofx3=[NaN; 86000; NaN; 77500];
% ofy values must be integers for some reason
ofy1=[0; 1; 1; 2]; ofy2=[-4; -4; -3; -3]; ofy3=[0; 0; 0; -1];
ofxy1=[73606; 81000; 80000; 60000]; ofxy2=[1; 1; 1; 22000]; ofxy3=[1; 1; 1; 41000];
n=[1; 1; 1; 1]; stoff=[0; 0; 1000; 9000];

% create a figure
figure(1)
clf

for i=1:4
    % generate GNSS slant times (st) using DOG approx position and GPS2RNG
    v=1500; depth=5225;
    st=gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'},'ave',xyz(i,:),v,depth);

    % plot GNSS slant times from ship positions
    subplot(4,1,i)
    plot(st,'b'); hold on

    % Load the dog tags
    load(char(fname(i)))

    % do the alignment of the acoustic data with the GNSS data

    % always want to unwrap before anything else:
    % correct where the big negative jumps in the dependent
    % variable are, they need an extra second, and where they are big
    % positive they stop needing that, so that's unwrapping; the result is in seconds
    % Unwrap disregards NaN's but maybe it shouldn't
    % Could INTERPOLATE over the NaNs for the unwrapping, then undo the interpolation
    tags(:,2)=unwrap(tags(:,2)/1e9*2*pi)/2/pi;
    % First, you eyeball an offset and cut it off
    tags=tags(ofx1(i,1):end,1:2)+ofy1(i,1);
    % plot once to see everything and get a baseline
    p(1)=plot(tags(:,2),'g.','MarkerSize',3);
    % Now need to stick in extra NaNs:
    % Shots are every seconds so difference between arrival times cannot be
    % much more than that, since the boat doesn't go that fast, and the
    % velocity of the water doesn't change that much. So whenever the second
    % skips by more than n, there must be a missed detection
    ntags=retimei(tags,n(i,1));
    % speed corrections:
    % end bit between ofxy1 and end
    ntags(ofxy1(i,1)+1:end,2)=ntags(ofxy1(i,1)+1:end,2)+ofy2(i,1);
    % middle bit between ofxy2 and ofxy3
    ntags(ofxy2(i,1):ofxy3(i,1),2)=ntags(ofxy2(i,1):ofxy3(i,1),2)+ofy3(i,1);
    % Eyeball a new offset and cut that off
    ntags=ntags(ofx2(i,1):end,1:2);
    % re-plot the new tags and you should be good
    if isnan(ofx3(i,1)) == 0
        p(2)=plot([1:ofx3(i,1)]+stoff(i,1),ntags(1:ofx3(i,1),2),'r.','MarkerSize',3);
    else
        p(2)=plot([1:length(ntags)]+stoff(i,1),ntags(:,2),'r.','MarkerSize',3);
    end
    % delete original plot after everything looks good    
    delete(p(1))
    hold off

    % cosmetics
    Nh=24;
    xlim([0 3600*Nh])
    ylim([0.85*min(st) 1.05*max(st)])
    longticks(gca,2)
    grid on
    xlabel('time [h]')
    nh=3;
    xticks(0:nh*3600:size(ntags,1))
    xticklabels(0:nh:floor(size(ntags,1)/3600))
    ylabel('slant range time [s]')
    % plot legend
    hold on
    pl(1)=plot(-1,-1,'b-');
    pl(2)=plot(-1,-1,'r-');
    hold off
    legs=legend(pl,'GNSS','Acoustic','FontSize',6,'Location','northeast');
end

keyboard

% print pdf
figdisp(sprintf('dog2rng-plot'),[],'',2,[],'epstopdf')