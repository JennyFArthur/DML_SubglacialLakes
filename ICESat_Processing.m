%% Code for processing ICESat GLA12 ice surface elevation data for Arthur et al. (submitted) 
%% 'Subglacial lake activity beneath the East Antarctic Ice Sheet in coastal Dronning Maud Land from ICESat-2 laser altimetry'.
% 
% ICESat files used in this script already contain elevation anomalies calculated using the plane-fitting method described in Moholdt et al. (2010) (Details: https://doi.org/10.1016/j.rse.2010.06.008).

% 1. Reads in individual ICESat-1 point measurements for individual tracks (from mat file).
% 2. Subsets track data to a specified Region of Interest.
% 3. Plots along-track elevation anomalies.

%% Read in ICESat elevation anomaly data (subset to ROI covering lake location)
clear all; close all;
% contains all point measurements used in the plane fitting (99% of all data) and their elevation anomaly (v) with respect to the plane fits. 
mat = matfile('G:\ExampleICESatfile.mat'); %Update filepath accordingly
ICESat1_data = struct();
ICESat1_data.date = mat.GLA(:, 1);
ICESat1_data.lon = mat.GLA(:, 2);
ICESat1_data.lat = mat.GLA(:, 3);
ICESat1_data.east = mat.GLA(:, 4);
ICESat1_data.north = mat.GLA(:, 5);
ICESat1_data.height = mat.GLA(:, 6);
ICESat1_data.v1 = mat.GLA(:, 34);
ICESat1_data.v2 = mat.GLA(:, 35);
ICESat1_data.dhdt = mat.GLA(:, 33);

%% Subset to ROI covering smaller region around potential lake
ICESat1_ROI = shaperead('C:\ExampleShapefilePath.shp'); %Update filepath accordingly
RoiBLake_ROI_X = ICESat1_ROI.X;
RoiBLake_ROI_Y = ICESat1_ROI.Y;
intersectsROI = inpolygon(ICESat1_data.east, ICESat1_data.north, RoiBLake_ROI_X, RoiBLake_ROI_Y);
ICESat1_ROI = struct();
ICESat1_ROI.date = ICESat1_data.date(intersectsROI);
ICESat1_ROI.east = ICESat1_data.east(intersectsROI);
ICESat1_ROI.north = ICESat1_data.north(intersectsROI);
ICESat1_ROI.dhdt = ICESat1_data.dhdt(intersectsROI);
ICESat1_ROI.v1 = ICESat1_data.v1(intersectsROI);
ICESat1_ROI.v2 = ICESat1_data.v2(intersectsROI);
ICESat1_ROI.height = ICESat1_data.height(intersectsROI);
ICESat1_ROI.max_mag_v = ICESat1_ROI.v1 .* (abs(ICESat1_ROI.v1) >= abs(ICESat1_ROI.v2)) + ... % determine maximum magnitude value out of v1 and v2
                        ICESat1_ROI.v2 .* (abs(ICESat1_ROI.v1) < abs(ICESat1_ROI.v2));
date = ICESat1_ROI.date;
east = ICESat1_ROI.east;
north = ICESat1_ROI.north;
dhdt = ICESat1_ROI.dhdt;
height = ICESat1_ROI.height;
maxmag = ICESat1_ROI.max_mag_v;
ICESat1_ROI = horzcat(east, north, date, dhdt, height, maxmag);

% Calculate along-track distance in km for plotting
diff_easting = diff(ICESat1_ROI.east);
diff_northing = diff(ICESat1_ROI.north);
pairwise_distance = sqrt(diff_easting.^2 + diff_northing.^2); %Calculate Euclidean distance between neighboring points
along_track_distance = cumsum([0; pairwise_distance]); % Add a 0 to beginning to represent starting point. Along-track dist in metres.
along_track_distance = along_track_distance/1000; % Convert to km.

% Date conversion for plotting
dates_numeric = str2double(string(ICESat1_ROI.date)); % Convert date strings to datetime objects
dates_datetime = datetime(dates_numeric, 'ConvertFrom', 'yyyymmdd');
dates_formatted = datestr(dates_datetime, 'yyyy/mm/dd'); % Convert datetime objects to strings in the desired format "YYYY/MM/DD"
ICESat1_ROI.date_formatted = dates_formatted; % Add the formatted dates as a new field in the struct

%% Plot ICESAt-1 elevation anomalies for all cycles
unique_dates = unique(ICESat1_data.date);
num_dates = numel(unique_dates);
colors = parula(num_dates);
hold on;
for i = 1:num_dates
    date_indices = find(ICESat1_ROI.date == unique_dates(i));
    scatter(along_track_distance(date_indices), ICESat1_ROI.max_mag_v(date_indices), 20, colors(i, :), 'filled');
end
hold off;
xlabel('Along-track distance (km)', 'FontSize', 16,'FontWeight', 'bold');
ylabel('Elevation anomaly (m)','FontSize', 16,'FontWeight', 'bold');
hAx = gca;
hAx.FontSize = 18; 
hAx.Box = 'on';
line(get(gca, 'XLim'), [0 0], 'Color', 'k');
grid on;
unique_dates = unique(ICESat1_ROI.date);
LegNames = cell(num_dates, 1);
for i = 1:num_dates
    date_str = num2str(unique_dates(i));
    year = date_str(1:4);
    month = date_str(5:6);
    day = date_str(7:8);
    LegNames{i} = [year '/' month '/' day];
end
legend(LegNames, 'FontSize', 16); 
title('ICESat-1 Track 21 (Roi Baudouin)');
