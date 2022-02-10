function tags=dog2mat(diro)
% tags=DOG2MAT(diro)
%
% Reads in all the *.txt files with data from the geodesy board
%
% INPUT
%
% diro        Some directory
% 
% OUTPUT:
%
% tags        The tags in two-column format
%
% EXAMPLE:
%
% Ran this, then saved the outputs in a *mat
%
% Last modified by fjsimons-at-alum.mit.edu, 02/07/2022

% Finds all the files
files=ls2cell(diro);

tags=[];
for index=1:length(files)
  d=load(files{index});
  tags=[tags ; d];
end

% Save to MAT

% Make a plot if there is no output? Watch out for implicit vector sums
%plot([tags(:,1)-tags(1,1)]*1e9+tags(:,2)-[1:size(tags(:,1),1)]')
%plot(tags(:,1)*1e9+tags(:,2))

% Optional output
varns={tags};
varargout=varns(1:nargout);

d e f i n i t e l y n o t r e a d y y e t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unwrap disregards NaN's but maybe it shouldn't
% Could INTERPOLATE over the NaNs for the unwrapping, then undo the interpolation

% Load the ranges with the "dead pixels" nan-ified - singles or jointly
% Note this is about right with respect to DOG1
st=gps2rng({'Unit1-camp.mat','Unit2-camp.mat','Unit3-camp.mat','Unit4-camp.mat'});

subplot(211)
plot(st,'b'); hold on
% Load the dog tags
load DOG1-camp
ofx1=5823;
% This next one is dependent on n below
ofx2=4550; n=2;
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

% Then you plot and look
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

% Then get rid of the old
delete(p(1))
hold off

Nh=24;
xlim([0 3600*Nh])
ylim([3 10])
longticks(gca,2)
grid on
xlabel('time [h]')
nh=3;
xticks(0:nh*3600:size(ntags,1))
xticklabels(0:nh:floor(size(ntags,1)/3600))
ylabel('slant range time [s]')

hold on
pl(1)=plot(-1,-1,'b-');
pl(2)=plot(-1,-1,'r-');
hold off
legs=legend(pl,'GNSS','Acoustic','FontSize',6,'Location','northeast');
moveh(legs,0.101)
movev(legs,0.01)

figdisp('dog2mat',[],[],2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These were bad ideas
good=0;
if good
  inta=floor(st);
  plot(st)
  hold on
  % This switches the tags to the left
  ofs=0
  for ofs=9750
    stags=ntags([ofs:min(size(ntags,1),size(inta,1))],:);
    p=plot(stags(1:end,2)/1e9+inta(1:size(stags,1)),'k.','MarkerS',3);
    title(num2str(ofs))                              
    pause
    delete(p)
  end

  % This switches to the tags to the right
  ofs=0
  for ofs=1:10:3000
    sinta=inta([ofs:min(size(ntags,1),size(inta,1))]);
    p=plot(ntags(1:length(sinta),2)/1e9+sinta,'k.','MarkerS',3);
    title(num2str(ofs))                              
    pause
    delete(p)
  end
end
