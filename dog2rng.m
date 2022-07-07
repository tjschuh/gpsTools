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
% dog2rng({'DOG1-camp.mat'})
%
% Originally written by fjsimons-at-princeton.edu, 02/07/2022
% Last modified by tschuh-at-princeton.edu, 07/07/2022

% Unwrap disregards NaN's but maybe it shouldn't
% Could INTERPOLATE over the NaNs for the unwrapping, then undo the interpolation

% Load the ranges with the "dead pixels" nan-ified - singles or jointly
% Note this is about right with respect to DOG1
st=gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'});

figure(1)
clf

% plot GNSS slant times from ship positions
plot(st,'b'); hold on

% Load the dog tags
load(char(fname))

if num == 1
    ofx1=5823;
    % This next one is dependent on n below
    ofx2=4800; n=1;
    % Another range for a different ofy
    ofxy1=73606;
    ofy1=3.1;
    ofy2=-0.925;

    % First, you eyeball an offset and cut it off
    tags=tags(ofx1:end,1:2);
    % And now you correct where the big negative jumps in the dependent
    % variable are, they need an extra second, and where they are big
    % positive they stop needing that, so that's unwrapping; the result is in seconds
    tags(:,2)=unwrap(tags(:,2)/1e9*2*pi)/2/pi;

    p(1)=plot(tags(:,2)+ofy1,'g.','MarkerSize',3);
    
    % Now need to stick in extra NaNs

    % Shots are every seconds so difference between arrival times cannot be
    % much more than that, since the boat doesn't go that fast, and the
    % velocity of the water doesn't change that much. So whenever the second
    % skips by more than n, there must be a missed detection
    ntags=retimei(tags,n);

    % Then, you eyeball a speed correction
    ntags(1:ofxy1,2)=ntags(1:ofxy1,2)+ofy1;
    ntags(ofxy1+1:end,2)=ntags(ofxy1+1:end,2)+ofy2;

    % Eyeball a new offset and cut that off
    ntags=ntags(ofx2:end,1:2);

    % And replot those new tags
    p(2)=plot(ntags(:,2),'r.','MarkerSize',3);
    keyboard
    delete(p(1))
    hold off
elseif num == 2
    tags(:,2)=unwrap(tags(:,2)/1e9*2*pi)/2/pi;
    p(1)=plot(tags(:,2),'g.','MarkerSize',3);

    ntags=retimei(tags,1);

    p(2)=plot(ntags(:,2),'r.','MarkerSize',3);
    keyboard
elseif num == 3
    ofx1=1; ofy1 = 1;
    ofx2=12000; ofy2=-3;
    n=1; ofxy1=80000;
    
    % unwrap first
    tags(:,2)=unwrap(tags(:,2)/1e9*2*pi)/2/pi;

    % then cut off beginning by eyeballing
    tags=tags(ofx1:end,1:2)+ofy1;

    p(1)=plot(tags(:,2),'g.','MarkerSize',3);

    % then retime
    ntags=retimei(tags,n);

    ntags(1:ofxy1,2)=ntags(1:ofxy1,2);
    ntags(ofxy1+1:end,2)=ntags(ofxy1+1:end,2)+ofy2;
    
    ntags=ntags(ofx2:end,1:2);

    p(2)=plot(ntags(:,2),'r.','MarkerSize',3);
    keyboard
    delete(p(1))
    hold off
    keyboard
elseif num == 4
    % load in DOG4
end
    
% cosmetics
Nh=24;
xlim([0 3600*Nh])
ylim([0.85*min(ntags(:,2)) 1.05*max(ntags(:,2))])
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
figdisp('dog2mat2',[],'',2,[],'epstopdf')    

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % These were bad ideas
% good=0;
% if good
%   inta=floor(st);
%   plot(st)
%   hold on
%   % This switches the tags to the left
%   ofs=0
%   for ofs=9750
%     stags=ntags([ofs:min(size(ntags,1),size(inta,1))],:);
%     p=plot(stags(1:end,2)/1e9+inta(1:size(stags,1)),'k.','MarkerS',3);
%     title(num2str(ofs))                              
%     pause
%     delete(p)
%   end

%   % This switches to the tags to the right
%   ofs=0
%   for ofs=1:10:3000
%     sinta=inta([ofs:min(size(ntags,1),size(inta,1))]);
%     p=plot(ntags(1:length(sinta),2)/1e9+sinta,'k.','MarkerS',3);
%     title(num2str(ofs))                              
%     pause
%     delete(p)
%   end
% end

