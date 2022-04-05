function dwplot(data,thresh)
% DWPLOT(data,thresh)
%
% description
%
% INPUT:
%
% data
% thresh
%
% Originally written by tschuh-at-princeton.edu, 04/05/2022

% p-value threshold
defval('thresh',0.05)

sz = 50;

% plot ellipsoid
figure(1)
for i=1:length(data)          
    if data(i,4) > thresh
        scatter3(data(i,1),data(i,2),data(i,3),sz,'filled','k')
        hold on
    end
end
xlabel('x')
ylabel('y')
zlabel('z')
title(sprintf('Locations where Durbin-Watson test result = 2 +/- %g',thresh))
grid on
longticks
hold off

% plot all values cube
figure(2)
scatter3(data(:,1),data(:,2),data(:,3),sz,data(:,6),'filled')
a = colorbar;
a.Label.String = 'p-value';
a.FontSize = 11;
xlabel('x')
ylabel('y')
zlabel('z')
grid on
longticks

% plot slice at z=0
slice = data(find(data(:,3) == 0),:);

figure(3)
scatter(slice(:,1),slice(:,2),sz,slice(:,6),'filled')
a = colorbar;
a.Label.String = 'p-value';
a.FontSize = 11;
xlabel('x')
ylabel('y')
grid on
longticks