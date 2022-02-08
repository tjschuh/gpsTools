function gps2plt(d,col,ovr,ext)
% GPS2PLT(d,col,ovr,ext)
%
% Plots a GNSS data structure, e.g., as created by PRD2MAT.
%
% INPUT:
%
% d        A GNSS data structure
% col      A RGB color
% ovr      A name override for the plot
% ext      An extension name for the plot
%
% Originally written by tschuh-at-princeton.edu, 11/12/2021
% Last modified by tschuh-at-princeton.edu, 02/03/2022
% Last modified by fjsimons-at-alum.mit.edu, 02/08/2022

% Use (RE)TIME(TABLE) to fill in time skips/data gaps with NaNs
drem={'xyzunit','lonlatunit','utmunit','heightunit','satlabels','utmzone'};
d=retimes(d,drem,{'secondly','fillwithmissing'});

% keep rows where nsats > nthresh and pdop < pthres and pdop~=0
nthresh = 4; pthresh = 15;

% Symbol size
sz=10;

% plotting interval
pint=20;
deltat=seconds(d.t(pint+1)-d.t(1));

% Average speed
vave=nanmean(sqrt([diff(d.utmeasting).^2+diff(d.utmnorthing).^2])./...
	     seconds(diff(d.t)))*3600/1000;

% subsample first, then look for condition
dkp={'t','xyz','lat','lon','utmeasting','utmnorthing','height','nsats','pdop'};
for index=1:length(dkp)
  % Don't forget that some are multicolumn, and keep the last one
  ssam=unique([1:pint:length(d.(dkp{index})) length(d.(dkp{index}))]);
  d.(dkp{index})=d.(dkp{index})(ssam,:);
end

% find the good data condition
cond=d.pdop<pthresh & d.pdop~=0 & d.nsats(:,1)>nthresh;

% Just open a figure and fiddle with positions outide
figure(1)
clf
% plot utm coordinates referenced to the minima
x=d.utmeasting -min(d.utmeasting);
y=d.utmnorthing-min(d.utmnorthing);
% and refer the times to the beginning of the sequence
t=d.t-d.t(1);
% Tick marks
ntix=3;
ttix=[1 round([1:ntix-1]*length(d(1).t)/(ntix-1))];
tixl=datestr(d.t(ttix),'HH:MM:SS');

% Only good data
gx=x; gx(~cond)=NaN;
gy=y; gy(~cond)=NaN;
gh=d.height; gh(~cond)=NaN;
% Only bad data
bx=x; bx(cond)=NaN;
by=y; by(cond)=NaN;
bh=d.height; bh(cond)=NaN;

% First panel - the ship track %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ah(1)=subplot(2,2,[1 3]);
% First the good data
sg=scatter(gx',gy',sz,seconds(t),'filled'); hold on
% Then the bad data
sb=scatter(bx',by',sz,grey,'filled'); hold off
xlabel('relative easting [m]')
ylabel('relative northing [m]')
tl(1)=title(sprintf('Ship Location (Every %dth Point)',pint));
grid on; box on
axis equal tight
% open it up a smidgen
xel=xlim; yel=ylim;
ah(1).XLim=xel+[-1 1]*range(xel)/20;
ah(1).YLim=yel+[-1 1]*range(yel)/20;
% Some minor annotations
tx(1)=text(0,0,sprintf('dop < %.0f\nsat > %.0f\n%s = %3.1f km/h\n%st = %i s',...
		       pthresh,nthresh,'v',vave,'\Delta',deltat),...
	   'VerticalAlignment','bottom');


longticks(ah)
% Then the color scale
colormap(jet)
cb=colorbar('southoutside','Ticks',seconds(t(ttix)),'TickLabels',tixl);
%shrink(cb,1,1.5)
longticks(cb)
% Give the color bar an xlabel with the day!
dat1=datestr(d.t(1),'dd mmm yyyy');
dat2=datestr(d.t(end),'dd mmm yyyy');
if dat1==dat2
  cb.XLabel.String=sprintf('%s',dat1);
else
  cb.XLabel.String=sprintf('%s to %s',dat1,dat2);
end

% Second panel - the elevation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ah(2)=subplot(2,2,2);
% First the good data
plot(d.t,gh,'color',col); hold on
% Then the bad data
plot(d.t,bh,'color',grey); hold off
xlim([d.t(1) d.t(end)])
xticks(d.t(ttix))
xticklabels(tixl)

% remove outliers so plotting looks better
% htout = rmoutliers(d.height,'mean');
% How much was removed, in percent
% outpct = (length(d.height)-length(htout))/length(d.height)*100;
% ylim([min(htout,[],'all')-0.005*abs(min(htout,[],'all')) max(htout,[],'all')+0.005*abs(max(htout,[],'all'))])
grid on
longticks
ylabel('height above WGS84 [m]')
tl(2)=title(sprintf('Ship Height (Every %dth Point)',pint));

% Third panel, nsats and pdop  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ah(3)=subplot(2,2,4);
yyaxis left
ps=plot(d.t,d.nsats(:,1),'b');
yticks([min(d.nsats(:,1))-1:max(d.nsats(:,1))+1])
ylim([min(d.nsats(:,1))-0.5 max(d.nsats(:,1))+0.5])
% ylabel('number of observed satellites')
ylabel('number of satellites')
yyaxis right
pp=plot(d.t,d.pdop,'r');
ylim([min(d.pdop)-0.25 max(d.pdop)+0.25])
xlim([d.t(1) d.t(end)])
xticks(d.t(ttix))
xticklabels(tixl)
%ylabel('position dilution of precision')
ylabel('dilution of precision')
% can only turn grid on for left axis
grid on
longticks
tl(3)=title('Total Number of Satellites and PDOP');

tl(4)=supertit(ah([1 2]),sprintf('Ship Data from %s to %s',d.t(1),d.t(end)));
movev(tl(4),0.3)

% Get rid of the ugly titles
delete(tl)
% Print to file
defval('ovr',mfilename)
figdisp(ovr,ext,[],2)
