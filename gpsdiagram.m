function gpsdiagram(dist12,dist13,dist14,dist23,dist24,dist34)
% GPSDIAGRAM(dist12,dist13,dist14,dist23,dist24,dist34)
%
% Make a simple diagram of 4 GPS locations given 6 distances between them
%
% INPUT:
%
% dist12      approximate distance between receivers 1 and 2 [m]
% dist13      approximate distance between receivers 1 and 3 [m]
% dist14      approximate distance between receivers 1 and 4 [m]
% dist23      approximate distance between receivers 2 and 3 [m]
% dist24      approximate distance between receivers 2 and 4 [m]
% dist34      approximate distance between receivers 3 and 4 [m]
%
% OUTPUT:
%
% plot of diagram
%
% Originally written by tschuh-at-princeton.edu, 01/12/2022
% Last modified by fjsimons-at-princeton.edu, 01/20/2022

% Values from R/V Atlantic Explorer
defval('dist12', 4.8)
defval('dist13',12.2)
defval('dist14',10.1)
defval('dist23',10.3)
defval('dist24',11.3)
defval('dist34', 7.1)

% Determine GPS relative coordinates relative to receiver 4
gps(1,:) = [0      dist14];
gps(2,:) = [dist12 dist14];
gps(3,:) = [dist34 0     ];
gps(4,:) = [0      0     ];

% Find out where the diagonals (1,3) and (2,4) cross and use that as a reference
% https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
D =[gps(1,1)-gps(3,1)]*[gps(2,2)-gps(4,2)]-[gps(1,2)-gps(3,2)]*[gps(2,1)-gps(4,1)];
xo=([gps(1,1)*gps(3,2)-gps(1,2)*gps(3,1)]*[gps(2,1)-gps(4,1)]...
    -[gps(1,1)-gps(3,1)]*[gps(2,1)*gps(4,2)-gps(2,2)*gps(4,1)])/D;
yo=([gps(1,1)*gps(3,2)-gps(1,2)*gps(3,1)]*[gps(2,2)-gps(4,2)]...
    -[gps(1,2)-gps(3,2)]*[gps(2,1)*gps(4,2)-gps(2,2)*gps(4,1)])/D;
        
% Percentage threshold for calculation approximation in if statement
pct = 0.05;

% Receiver symbols and colors
sym='o'; 
me='k';
mf=grey;

% Font size of text boxes
fs=10;
% Offset for these text boxes
to=[0 1 ; 0 1 ; 0 -1 ; 0 -1]/1.5;
% Offset for axis limits
offset = [1.5 0.75];

% Special box strings
ts={'bow','stern','port','starboard'};
% Offset for these boxes
os=[0 3 ; 0 -5 ; -1.5 0 ; 1.25 0];
% Alignment for these boxes
as={'center','center','right','left'};

% if pythagorean thm holds for the given distances, then this
% must be the configuration we want so go ahead and plot
if abs((dist24)^2 - ((dist12)^2 + (dist14)^2)) <= pct*(dist24)^2 ...
        && abs((dist13)^2 - ((dist14)^2 + (dist34)^2)) <= pct*(dist13)^2
    clf
    hold on
    % Plot the straight connecting lines
    ps=plot(gps([1 2 3 4 1],1),gps([1 2 3 4 1],2),'k');
    % Plot the diagonal connecting lines
    pd=plot(gps([1 3 2 4],1),gps([1 3 2 4],2),'k');
    
    % Plotting of the GPS statioons, on top!
    hold on
    for index=1:4
      pg=plot(gps(index,1),gps(index,2),sym,'MarkerEdgeColor',me,'MarkerFaceColor',mf);
    end
    
    % Text boxes identifying the receivers 
    for index=1:size(gps,1)
      xa(index)=boxtex(gps(index,:)+to(index,:),[],index,fs,2);
    end
    hold off

    % grid on
    % if grid on neex to use BOXTEX but take the boxes away... so the
    % white covers the grid lines....

    % Now identify the nautical terms relative to everything else
    for index=1:length(ts)
      ta(index)=text(xo+os(index,1),yo+os(index,2),ts{index});
      set(ta(index),'HorizontalAlignment',as{index},'FontSize',fs+2)
    end    
    
    % Cosmetics and cleanup
    xlim([0 gps(3,1)]+[-1 1]*offset(2))
    ylim([0 gps(1,2)]+[-1 1]*offset(1))
    xlabel('distance [m]')
    ylabel('distance [m]')
    longticks(gca,2)
    box on
    
    % Didn't like the first and last two at some point
    % xt=get(gca,'XTick'); set(gca,'XTick',xt(2:end-1));
else
    error('Cannot compute configuration currently')
end

% Print to PDF
figdisp([],[],[],2,[],'epstopdf')