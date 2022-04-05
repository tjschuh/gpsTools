function dwplot(data,thresh)
% ADD HEADER

defval('thresh',0.5)

sz = 50;

% plot ellipsoid
figure(1)
for i=1:length(data)          
    if data(i,4) > 2-thresh & data(i,4) < 2+thresh
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

% plot all values cube
figure(2)
scatter3(data(:,1),data(:,2),data(:,3),sz,data(:,4),'filled')
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
scatter(slice(:,1),slice(:,2),sz,slice(:,4),'filled')
a = colorbar;
a.Label.String = 'p-value';
a.FontSize = 11;
xlabel('x')
ylabel('y')
grid on
longticks