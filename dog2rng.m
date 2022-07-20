function dog2rng(fname,num)
% DOG2RNG(fname,num)
%
% Works on a file made by DOG2MAT and turns it into a range measurement
%
% INPUT:
%
% fname         file made by DOG2MAT with time tags
% num           dog number corresponding to file (temporary input)
%
% EXAMPLE:
%
% dog2rng({'DOG1-camp.mat'},1)
%
% Originally written by fjsimons-at-princeton.edu, 02/07/2022
% Last modified by tschuh-at-princeton.edu, 07/13/2022

% make this much cleaner and lose duplicate code

% Load the ranges with the "dead pixels" nan-ified - singles or jointly
if num == 1
    % Note this is about right with respect to DOG1
    xyz=[1.977967 -5.073198 3.3101016]*1e6;

    ofx1=5823;
    % This next one is dependent on n below
    ofx2=4800; n=1;
    % Another range for a different ofy
    ofxy1=73606;
    ofy1=0;
    ofy2=-4; stoff=0;
    ofy3=0;
    ofxy2=1; ofxy3=1;
elseif num == 2
    % DOG2 approx sea-surface xyz position (this is about right)
    xyz=[1.9765409 -5.074798706 3.308558817]*1e6;

    ofx1=1; ofy1=1;
    ofx2=12500; ofy2=-4;
    n=1; ofxy1=81000; ofx3=86000;
    stoff=0;
    ofy3=0;
    ofxy2=1; ofxy3=1;
elseif num == 3
    % DOG3 approx position (this is about right)
    xyz=[1.980083296 -5.073417583 3.308558817]*1e6;

    ofx1=1; ofy1 = 1;
    ofx2=12000; ofy2=-3;
    n=1; ofxy1=80000;
    stoff=1000;
    ofy3=0;
    ofxy2=1; ofxy3=1;
elseif num == 4
    % DOG4 approx position (this is about right)
    xyz=[1.978732384 -5.075186123 3.30666684]*1e6;

    ofx1=2250; ofy1=2; %ofy1 must be integer for some reason
    ofx2=1; ofy2=-3;
    n=1; ofxy1=60000;
    ofxy2=22000; ofxy3=41000;
    ofy3=-1;
    ofx3=77500; stoff=9000;
end

% generate GNSS slant times (st) using DOG approx position and GPS2RNG
v=1500; depth=5225;
st=gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'},'ave',xyz,v,depth);

% create a figure
figure(1)
clf

% plot GNSS slant times from ship positions
plot(st,'b'); hold on

% Load the dog tags
load(char(fname))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do the alignment of the acoustic data with the GNSS data

% always want to unwrap before anything else:
% correct where the big negative jumps in the dependent
% variable are, they need an extra second, and where they are big
% positive they stop needing that, so that's unwrapping; the result is in seconds
% Unwrap disregards NaN's but maybe it shouldn't
% Could INTERPOLATE over the NaNs for the unwrapping, then undo the interpolation
tags(:,2)=unwrap(tags(:,2)/1e9*2*pi)/2/pi;
% First, you eyeball an offset and cut it off
tags=tags(ofx1:end,1:2)+ofy1;
% plot once to see everything and get a baseline
p(1)=plot(tags(:,2),'g.','MarkerSize',3);
% Now need to stick in extra NaNs:
% Shots are every seconds so difference between arrival times cannot be
% much more than that, since the boat doesn't go that fast, and the
% velocity of the water doesn't change that much. So whenever the second
% skips by more than n, there must be a missed detection
ntags=retimei(tags,n);
% speed corrections:
% end bit between ofxy1 and end
ntags(ofxy1+1:end,2)=ntags(ofxy1+1:end,2)+ofy2;
% middle bit between ofxy2 and ofxy3
ntags(ofxy2:ofxy3,2)=ntags(ofxy2:ofxy3,2)+ofy3;
% Eyeball a new offset and cut that off
ntags=ntags(ofx2:end,1:2);
% re-plot the new tags and you should be good
if exist('ofx3','var') == 1
    p(2)=plot([1:ofx3]+stoff,ntags(1:ofx3,2),'r.','MarkerSize',3);
else
    p(2)=plot([1:length(ntags)]+stoff,ntags(:,2),'r.','MarkerSize',3);
end
% delete original plot after everything looks good    
delete(p(1))
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
keyboard
% print pdf
figdisp(sprintf('unit%i-st',num),[],'',2,[],'epstopdf')