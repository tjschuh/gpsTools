function tags=dog2rng(fname)
% tags=DOG2RNG(fname)
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
% Last modified by tschuh-at-princeton.edu, 07/26/2022

% add offset locations (start and end points where NaNs creep in)
% have this function call gps2syn/dwplot

% if not 4 DOGS used in input, dont run code
%if length(fname) != 4
%    error('Must have exactly 4 sets of DOG data')
%end

% DOG 1-4 approx positions
xyz=[1.977967 -5.073198 3.3101016; 1.9765409 -5.074798706 3.308558817; ...
     1.980083296 -5.073417583 3.308558817; 1.978732384 -5.075186123 3.30666684].*1e6;
% ofx1=how much to cut off at the beginning of the data before retimei is used
% ofx2=how much more to cut off at the beginning of the data after retimei is used
% ofx3=how much to cut off at the end of the data during plotting
ofx1=[5823; 1; 1; 2250]; ofx2=[4800; 12500; 12000; 1]; ofx3=[NaN; 86000; NaN; 77500];
% ofy values must be integers for some reason
% ofy values are integer shifts in the y-direction of sections of the data
ofy1=[0; 1; 1; 2]; ofy2=[-4; -4; -3; -3]; ofy3=[0; 0; 0; -1];
% ofxy values are locations in the data to separate at and perform ofy shifts on
ofxy1=[72856; 80678; 79625; 59490]; ofxy2=[1; 1; 1; 22000]; ofxy3=[1; 1; 1; 41000];
% n=threshold for NaN replacements-->if more than n consecutive missed detections, add NaNs
% stoff=st offset to be used during plotting so GNSS and acoustic data line up
n=[1; 1; 1; 1]; stoff=[-100; 0; 900; 8800];
stoff2=[17100; 17200; 18500; 17100];

% create a figure
figure(1)
clf

for i=1:4
    % generate GNSS slant times (st) using DOG approx position and GPS2RNG
    v=1500; depth=5225;
    st(:,i)=gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'},'ave',xyz(i,:),v,depth);

    % plot GNSS slant times from ship positions
    ah(i)=subplot(4,1,i);
    plot(st(:,i),'b'); hold on

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
    p(i,1)=plot(tags(:,2),'g.','MarkerSize',3);
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
    % add in st offset for plot alignment
    ntags(:,1)=ntags(:,1)-stoff2(i,1);
    % if we have non-NaN ofx3 value, build that in
    if isnan(ofx3(i,1)) == 0
        ntags=ntags(1:ofx3(i,1),1:2);
    end
    % re-plot the new tags and you should be good
    p(i,2)=plot(ntags(:,1),ntags(:,2),'r.','MarkerSize',3);
    % delete original plot after everything looks good    
    delete(p(i,1))
    hold off

    % to normalize the ntags data and have each dataset "count" the same way
    % take the actual offset spot (ofxy1) and add the st offset used during plotting (stoff)
    % and then subtract the bit used to cut off the end of the original ntags (ofx2)
    % this number signifies where NaNs (approximately) end and data resumes (in hrs)
    ofst(i,1) = (ofxy1(i,1) + stoff(i,1) - ofx2(i,1))/3600;
    % these numbers are what we want ofst to equal
    %ofst=[68056; 68178; 68625; 68489];
    % these numbers are where NaNs begin and data cuts out (not implemented yet)
    %[?; ?; ?; 68859];
    
    % cosmetics
    t(i)=title(sprintf('C-DOG approx position: x = %4.3f km, y = %4.3f km, z = %4.3f km',xyz(i,1)*1e-3,xyz(i,2)*1e-3,xyz(i,3)*1e-3),'FontSize',8);
    Nh=24;
    xlim([0 3600*Nh])
    ylim([0.85*min(st(:,i)) 1.05*max(st(:,i))])
    longticks(gca,4)
    grid on
    nh=3;
    xticks(0:nh*3600:size(ntags,1))    
    xticklabels({})    
    if i == 4
        xlabel('time [h]')
        xticklabels(0:nh:floor(size(ntags,1)/3600))
    end
    %ylabel('slant range time [s]')
    txt(i)=text(900,0.9*ah(i).YLim(2),sprintf('offset at %2.2f hrs ',ofst(i,1)),'FontSize',7);
    % plot legend
    hold on
    pl(1)=plot(-1,-1,'b-');
    pl(2)=plot(-1,-1,'r-');
    hold off
    legs=legend(pl,'GNSS','Acoustic','FontSize',4,'Location','northeast');
    % get subplot positions for single ylabel workaround
    g(i,:)=get(ah(i),'position');
end

% workaround to make a single ylabel on plot
height=g(1,2)+g(1,4)-g(4,2)-1.3*g(3,2);
width=g(4,1)+g(4,3)-g(3,1);
gg=axes('position',[g(3,1) g(3,2) width height],'visible','off'); 
gg.XLabel.Visible='off';
gg.YLabel.Visible='on';
axes(gg)
ylabel('slant range time [s]')

keyboard

% print pdf
figdisp(sprintf('dog2rng-plot'),[],'',2,[],'epstopdf')