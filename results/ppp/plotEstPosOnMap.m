function plotEstPosOnMap(pos_ecef, pos_error)

% Convert ECEF to latitude, longitude, and altitude
lla = ecef2lla(pos_ecef');

% Create a figure
figure;

% Create a colormap for the errors
colors = [0 1 0; 0 1 1; 0.5 0 0.5; 1 0 0; ]; % green, blue, purple, red

% Create a color scale based on the error
c = discretize(pos_error, [0, 1, 3, 80, inf]);

% Create the geographic axes
ax = geoaxes;
geobasemap(ax, 'satellite');

% Plot the positions with color based on error
hold on
for i = 1:max(c)
    idx = c == i;
    geoscatter(ax, lla(idx,1), lla(idx,2), 10, colors(i,:), 'filled');
end
hold off

% Add a colorbar
colormap(colors)
cb = colorbar;
cb.Ticks = [1/8, 3/8, 5/8, 7/8]; % Adjust these values as needed
cb.TickLabels = {'< 1.0m', '1.0m-3m', '3-80m', '>80m'};

% Enable data cursor mode
dcm_obj = datacursormode(gcf);
set(dcm_obj, 'UpdateFcn', {@myupdatefcn, lla});
end

function txt = myupdatefcn(~, event_obj, lla)
% Custom function to display index when a point is clicked
pos = event_obj.Position;
idx = find(lla(:,1) == pos(1) & lla(:,2) == pos(2), 1);
if ~isempty(idx)
    txt = ['Index: ', num2str(idx)];
else
    txt = '';
end
end