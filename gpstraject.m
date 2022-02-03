function gpstraject
% GPSTRAJECT
%
% Make a beautiful diagram of the seafloor geodetic campaign
%
% INPUT:
%
% OUTPUT:
%
% Last modified by fjsimons-at-princeton.edu, 01/22/2022

% Top topo map in [0 360] and [-90 90]
xels=[278 296];
yels=[25 35];

% Receiver symbols and colors
sym={'.','o','.'}; 
me={'r','k','g'};
mf={'r','k','g'};
ms=[2 2 2];

% Get one set of data; execute on ariel
diro='/data1/seafloorgeodesy/SwiftNavData/BIOSCruise/Unit2';
% What the data file is, actually
fname=sprintf('%s_%s.mat',mfilename,suf(diro,'/'));

try 
  % Get it
  load(fname)
  % Make sure that's saved as well (wasn't at first)
  defval('dirs',{'leg1','camp','leg2',});
catch
  % This for Ariel
  addpath('/home/fjsimons/PROGRAMS/MFILES/slepian_alpha')
  % Make it
  for index=1:length(dirs)
    mats=ls2cell(fullfile(diro,dirs{index},'*.mat'));
    % Every round minute is our target, alternative would be every minute apart
    bign=60*length(mats);
    % Initialize 
    t.(dirs{index})=NaT(bign,1);
    hwgs.(dirs{index})=nan(bign,1);
    lonlat.(dirs{index})=nan(bign,2);
    utmeno.(dirs{index})=nan(bign,2);
    for ondex=1:length(mats)
      % Load the file which renders the structure 'd'
      load(fullfile(diro,dirs{index},mats{ondex}))
      % Running index
      runin=(ondex-1)*60+1:ondex*60;
      % We're only saving by the round minute, how about that
      lojik=second(d.t)==0;
      runon=1:sum(lojik);
      % Proper assignment with makeups for nan
      t.(dirs{index})(runin(runon))=d.t(lojik);
      hwgs.(dirs{index})(runin(runon))=d.height(lojik);
      lonlat.(dirs{index})(runin(runon),:)=[d.lon(lojik) d.lat(lojik)];
      utmeno.(dirs{index})(runin(runon),:)=[d.utmeasting(lojik) d.utmnorthing(lojik)];
    end
    % After that, remove all the nans? or don't bother
  end
  % Save as MAT file, careful, be explicit with extension
  save(fname,'t','hwgs','lonlat','utmeno','dirs')
end

% Make the figure
clf
ah=krijetem(subnum(2,1));
aha(1,:)=getpos(ah(1));
clf
ah=krijetem(subnum(2,2));
aha(2,:)=getpos(ah(3));
aha(3,:)=getpos(ah(4));
clf
clear ah
for index=1:size(aha,1)
  ah(index)=axes('position',aha(index,:));
end

% Overview map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(1))

% Get GEBCO bathymetry with one interpolated color bar
[z,C11,CMN]=getco(xels,yels);

% Plot GEBCO bathymetry
cax=[-6000 500];
[~,cm]=cax2dem(cax,'ver');
[cb(1),xcb(1)]=plotz(z,C11,CMN,cax,cm);
cbarticks(cb(1),cax,unique([0 cax(1):1000:cax(2) cax(2)]),'vert');
xcb(1).String='GEBCO depth (m)';

% Skip
skp=60;

% Don't forget to return to the main axes
axes(ah(1))
% Only plot the legs, not the camp
for index=[1 3]
  hold on
  pl(index)=plot(lonlat.(dirs{index})(1:skp:end,1),lonlat.(dirs{index})(1:skp:end,2),...
       'Marker',sym{index},'MarkerFaceColor',mf{index},'MarkerEdgeColor',me{index},...
		 'MarkerSize',ms(index)+2,'LineStyle','none');
end
[~,pc]=plotcont(C11+[-1 1]*5,CMN+[1 -1]*5);
% Maybe adjust to common field of view regardless of resolution even if
% that means cutting pixels off
axis([xels yels])
% xlim([C11(1) CMN(1)]); ylim([CMN(2) C11(2)])
hold off
delete(pc)
xl(1)=xlabel(sprintf('longitude (%s)',176));
yl(1)=ylabel(sprintf('latitude (%s)',176));

% Detail map 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(3))

% Skip
skp=10;

% Just the campaign
index=2;
zex=utmeno.(dirs{index})(1:skp:end,1);
zwi=utmeno.(dirs{index})(1:skp:end,2);
refx=min(zex);
refy=min(zwi);
sclx=1000;
scly=1000;
pl(index)=plot((zex-refx)/sclx,(zwi-refy)/scly,...
	       'Marker',sym{index},'MarkerFaceColor',mf{index},'MarkerEdgeColor',me{index},...
	       'MarkerSize',ms(index),'LineStyle','-','Color',mf{index});
box on

grid on
xl(2)=xlabel('easting (km)');
yl(2)=ylabel('northing (km)');
axis equal tight
openup(ah(3),5,10)
openup(ah(3),6,10)
% Plot box on first graph, didn't bother to save zone
B=axis; utmz='19 R';
% Put the reference and the scaling back
B(1:2)=B(1:2)*sclx;
B(3:4)=B(3:4)*scly;
B=B+[refx refx refy refy];
warning off MATLAB:nargchk:deprecated
[lat,lon]=utm2deg(B([1 2 2 1 1]),B([4 4 3 3 4]),repmat(utmz,5,1));
warning on MATLAB:nargchk:deprecated
axes(ah(1)); hold on
pb=plot(lon+360,lat,'LineWidth',1,'Color',mf{2});
hold off

% Detail map 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(ah(2))

% Get GEBCO bathymetry
[z,C11,CMN]=getco(minmax(lon)+360,minmax(lat));

% Plot GEBCO bathymetry with a different color bar
cax=round(minmax(z)/10)*10;
cax=[-5320 -5180];

[cb(2),xcb(2)]=plotz(z,C11,CMN,cax,kelicol);
% You don't have complete liberty to rescale as the scale works on the
% orginal ticks, and there sometimes is a missing tick that you need to
% make appear again...
cbarticks(cb(2),cax,20,'vert');
ylim(ylim+[-1 1]*range(ylim)/1000)
xcb(2).String='GEBCO depth (m)';
axes(ah(2))
xl(2)=xlabel(sprintf('longitude (%s)',176));
yl(2)=ylabel(sprintf('latitude (%s)',176));

% Only plot the camp
index=2; hold on
plo=plot(lonlat.(dirs{index})(1:skp:end,1),lonlat.(dirs{index})(1:skp:end,2),...
	 'Marker',sym{index},'MarkerFaceColor',mf{index},'MarkerEdgeColor',me{index},...
	 'MarkerSize',ms(index),'LineStyle','-','Color',mf{index});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cosmetics
longticks(ah(1),2)
longticks(ah(2:3))

ah(3).YAxisLocation='right';

moveh(ah(2),.105)
moveh(cb(2),.105)
delete(xcb(2))

% Print to PDF
figdisp([],[],[],2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z,C11,CMN]=getco(xels,yels)
% Get the bathymetry
C11=[xels(1) yels(2)];
CMN=[xels(2) yels(1)];
spc=1/60/4; vers=2019; npc=20;
% Make the save file
defval('savefile',fullfile(getenv('IFILES'),'TOPOGRAPHY','EARTH','GEBCO','HASHES',...
			   sprintf('%s.mat',hash([C11 CMN vers npc spc],'SHA-1'))))
if exist(savefile)~=2
  % I cannot remember why LAT should count backwards according to GEBCO
  [LON,LAT]=meshgrid([C11(1):spc:CMN(1)],[C11(2):-spc:CMN(2)]);
  LON=LON-360;
  z=gebco(LON,LAT,vers,npc);
  save(savefile,'z')
else
  load(savefile)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cb,xcb]=plotz(z,C11,CMN,cax,cm);
% Different than in POLYNESIA for two colorbars
% See IMAGEF for comments on pixel registration, which these are

imagefnan(C11,CMN,z,cm,cax)
% Or else this
[cb,xcb]=addcb('ver',cax,cax,cm);
cb.YAxisLocation='right';

shrink(cb,1.5,1)
moveh(cb,0.04)
longticks(cb,2)
