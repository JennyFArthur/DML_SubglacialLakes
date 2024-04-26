%% Code for processing ICESat-2 ATL11 ice surface elevation data for Arthur et al. (submitted) 
%% 'Subglacial lake activity beneath the East Antarctic Ice Sheet in coastal Dronning Maud Land from ICESat-2 laser altimetry'.
% 
% 1. Reads in individual ICESat-2 track (hdf5 file).
% 2. Subsets track data to a specified Region of Interest.
% 3. Calculates elevation anomalies
% 4. Plots along-track elevation anomalies.
% 5. Plots along-track elevation anomalies alongside REMA surface elevations and BedMachine bed elevations.
% 6. Estimates lake volume (based on user-defined lake boundary shapefile,e.g. derived mnanually-derived from REMA strip differencing).


%% Read in ICESat-2 ATL11 elevation data
clear all; close all;
hdf_file = 'C:\ExampleFilePath\ATL11_061610_0319_006_04.h5'; %Update filepath accordingly

%Read in datasets of interest from pt1 group. Include extras as needed.
lat=h5read(hdf_file,'/pt1/latitude');
lon=h5read(hdf_file,'/pt1/longitude');
h_corr=double(h5read(hdf_file,'/pt1/h_corr'))'; %WGS-84 height, corrected for the ATL11 surface shape
h_corr_sigma=double(h5read(hdf_file,'/pt1/h_corr_sigma'))'; %Per-point errors in height estimates. Includes the errors related to the accuracy of the reference surface and the precision of the ICESat-2 range estimates, which are uncorrelated between adjacent reference points.
h_corr_sigma_systematic =double(h5read(hdf_file,'/pt1/h_corr_sigma_systematic'))'; %includes the contribution of uncertainties in measurement geolocation and the satellite's radial orbit errors to the measurement errors
quality_summary=double(h5read(hdf_file,'/pt1/quality_summary'))'; %data quality flag, 1 or 0 (1 = cloudy, rough or low surface-reflectance).
cycle_number=double(h5read(hdf_file,'/pt1/cycle_number'));
rgt=double(h5read(hdf_file,'/pt1/crossing_track_data/rgt')); %Reference Ground Track
rgt = rgt(find(lat)); %Find all indices of 'lat' and extract 'rgt' for only those indices.
ref_pt=double(h5read(hdf_file,'/pt1/ref_pt')); %Segment_id corresponding to center of ATL06 data used for each ATL11 point; 20m spacing.
delta_time=h5read(hdf_file,'/pt1/delta_time'); %Mean number of GPS seconds since 2018-01-01 (ATLAS SDP epoch).
atlas_sdp_gps_epoch=h5read(hdf_file,'/ancillary_data/atlas_sdp_gps_epoch'); %Add to delta_time to get full gps_seconds for each data point. Number of GPS seconds between GPS epoch (1980-01-06T00:00:00.000000Z UTC) and  ATLAS Standard Data Product (SDP) epoch (2018-01-01:T00.00.00.000000 UTC).
fit_quality=h5read(hdf_file,'/pt1/ref_surf/fit_quality'); %0: no problem identified, 1: One or more polynomial coefficients has an error of 10 or larger, 2: One or more surface slope components is greater than 0.02, 3: both 1 and 2
dem_h=h5read(hdf_file,'/pt1/ref_surf/dem_h'); % DEM elevation, derived from ATL06 (derived from REMA DEM filtered to 40-m resolution before interpolation to ICESat-2 segment centers.

%Convert lat,lon to x,y
map_proj = 'spolar';
[E,N] = geod2utm(lon,lat,map_proj); %Reproject to Antarctic Polar Stereographic (EPSG:3031) projection.

%Convert GPS delta_time to UTC time
epoch=datetime('2018-01-01', 'InputFormat', 'yyyy-MM-dd'); % Define start time
utcDateTime = epoch + seconds(delta_time); % Convert GPS delta_time to UTC time
utcDateTime = utcDateTime'; % Transpose to match other arrays
utcDateNum = datenum(utcDateTime); % Convert to datenum
utcDateNum(isinf(utcDateNum)) = NaN; %set any 'Inf' values to NaN
utcDateStr = datestr(utcDateNum,'yyyy-mm-dd-HH-MM-SS');

%Build data matrix, where each h_corr height value (1 per cycle) is a separate attribute
ATL = horzcat(E, N, h_corr, rgt,dem_h,h_corr_sigma,h_corr_sigma_systematic);

%Replace any NoData h_corr values with NaN
h_corr_columns = 3:21;
for b = h_corr_columns
    ATL(ATL(:, b) > 100000, b) = NaN;
end 

% Check maximum per-point and systematic errors
% h_corr_sigma(h_corr_sigma > 100000) = NaN;
% max(h_corr_sigma)
% h_corr_sigma_systematic(h_corr_sigma_systematic > 100000) = NaN;
% max(h_corr_sigma_systematic)

%% Create a new structure from ATL11 hdf-5 file data matrix for plotting
ATL_data = struct();
ATL_data.E = ATL(:, 1); % Assign variables to the structure fields
ATL_data.N = ATL(:, 2);
ATL_data.cycle_3 = ATL(:, 3);
ATL_data.cycle_4 = ATL(:, 4);
ATL_data.cycle_5 = ATL(:, 5);
ATL_data.cycle_6 = ATL(:, 6);
ATL_data.cycle_7 = ATL(:, 7);
ATL_data.cycle_8 = ATL(:, 8);
ATL_data.cycle_9 = ATL(:, 9);
ATL_data.cycle_10 = ATL(:, 10);
ATL_data.cycle_11 = ATL(:, 11);
ATL_data.cycle_12 = ATL(:, 12);
ATL_data.cycle_13 = ATL(:, 13);
ATL_data.cycle_14 = ATL(:, 14);
ATL_data.cycle_15 = ATL(:, 15);
ATL_data.cycle_16 = ATL(:, 16);
ATL_data.cycle_17 = ATL(:, 17);
ATL_data.cycle_18 = ATL(:, 18);
ATL_data.cycle_19 = ATL(:, 19);
ATL_data.RGT = ATL(:, 20);
ATL_data.dem = ATL(:, 21);

%Calculate along-track distance in km for plotting
diff_E = diff(ATL_data.E); %Calculate differences in E and N points
diff_N = diff(ATL_data.N);
pairwise_distance = sqrt(diff_E.^2 + diff_N.^2); %Calculate Euclidean distance between neighboring points
along_track_distance2 = cumsum([0; pairwise_distance]); % Add a 0 to beginning to represent starting point. Along-track dist in metres.
along_track_distance2 = along_track_distance2/1000; % Convert to km.


%% Read in Region Of Interest and subset along-track distance and ice-surface heights to ROI for plotting elevation profile
roi=shaperead('C:\ExampleShapefilePath.shp'); %Update filepath accordingly
roiGeometry=roi.Geometry;
roiX = roi.X;
roiY = roi.Y;
points=[ATL_data.E, ATL_data.N];
intersectsROI = inpolygon(ATL_data.E, ATL_data.N, roiX, roiY); %selects which points intersect with polygon
selectedE = E(intersectsROI); % Filter E based on intersection
selectedN = N(intersectsROI); % Filter N based on intersection
% Re-calculate along-track distance
diff_E = diff(selectedE); %Calculate differences in E and N points
diff_N = diff(selectedN);
pairwise_distance = sqrt(diff_E.^2 + diff_N.^2); %Calculate Euclidean distance between neighboring points
along_track_distance_new = cumsum([0; pairwise_distance]); % Add a 0 to beginning to represent starting point. Along-track dist in metres.
along_track_distance_new = along_track_distance_new/1000; % Convert to km.

roi_cycle3 = ATL_data.cycle_3(intersectsROI);
roi_cycle4 = ATL_data.cycle_4(intersectsROI);
roi_cycle5 = ATL_data.cycle_5(intersectsROI);
roi_cycle6 = ATL_data.cycle_6(intersectsROI);
roi_cycle7 = ATL_data.cycle_7(intersectsROI);
roi_cycle8 = ATL_data.cycle_8(intersectsROI);
roi_cycle9 = ATL_data.cycle_9(intersectsROI);
roi_cycle10 = ATL_data.cycle_10(intersectsROI);
roi_cycle11 = ATL_data.cycle_11(intersectsROI);
roi_cycle12 = ATL_data.cycle_12(intersectsROI);
roi_cycle13 = ATL_data.cycle_13(intersectsROI);
roi_cycle14 = ATL_data.cycle_14(intersectsROI);
roi_cycle15 = ATL_data.cycle_15(intersectsROI);
roi_cycle16 = ATL_data.cycle_16(intersectsROI);
roi_cycle17 = ATL_data.cycle_17(intersectsROI);
roi_cycle18 = ATL_data.cycle_18(intersectsROI);
roi_cycle19 = ATL_data.cycle_19(intersectsROI);
roi_dem = ATL_data.dem(intersectsROI); 
selectedRGT = rgt(intersectsROI);


%% Calculate elevation anomalies for each cycle
h_corr = [roi_cycle3, roi_cycle4, roi_cycle5, roi_cycle6, roi_cycle7, roi_cycle8, roi_cycle9, roi_cycle10, roi_cycle11, roi_cycle12, roi_cycle13, roi_cycle14, roi_cycle15, roi_cycle16, roi_cycle17, roi_cycle18, roi_cycle19];
elevation_anomaly = h_corr - roi_cycle3; %calculating elev anomalies as relative to first cycle
median_elevanom = nanmedian(elevation_anomaly, 2); %median elevation anomaly across all cycles
ROI = inpolygon(selectedE, selectedN, roiX, roiY);
validTime = utcDateTime(ROI,:);
%Calculate running mean of elevation anomalies for smoothing along-track elev anomaly plots
window_size = 10; %number of data points used for calculating running mean
running_mean = zeros(size(elevation_anomaly)); %initialise array
for i = 1:numel(elevation_anomaly)
    start_index = max(1, i - window_size + 1);
    end_index = i;
    running_mean(i) = mean(elevation_anomaly(start_index:end_index), 'omitnan');
end


%% Plot ICESat-2 along-track elevation profile for all cycles
figure;
hold on;
numColors = size(h_corr, 2); %calculate number of columns in h_corr
colors = parula(numColors);  % Generate a custom colormap with distinct colors
%Legend set-up
LegNames = cell(13,1);
for i = 1:size(h_corr, 2) %plot h_corr for each cycle
    scatter(along_track_distance_new, h_corr(:, i), 19, colors(i, :), 'filled');
    dateInd = find(~isnan(utcDateNum(:,i)), 1); %find first non-NaN instance of DateNum in column i  for inclusion in Legend
    LegName = ['Cycle ', num2str(cycle_number(i)), ': ', datestr(utcDateTime(dateInd, i), 'dd mmm yyyy')];
    LegNames{i} = LegName;
end
hold off;
legend(LegNames, 'Location', 'northeast');
leg = legend('show');
leg.FontSize = 18; 
xlabel('Along-track distance (km)','FontSize', 16,'FontWeight', 'bold');
ylabel('Elevation (m)','FontSize', 16,'FontWeight', 'bold');
hAx = gca;
hAx.FontSize = 18; 
hAx.Box = 'on';
hAx.LineWidth = 1;
hAx.TickLength = [0.005, 0.1];
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
grid on;
%print(gcf, 'alongtrack_elev_allcycles_ATL11_108212_0315_005_03.png', '-dpng', '-r300');
title('ICESat-2 RGT 990, 004512 Along-track elevation (Region C)');


%% Plot ICESat-2 along-track elevation profile for all cycles alongside REMA along-track elevation profile
figure;
hold on;
numColors = size(h_corr, 2); %calculate number of columns in h_corr
colors = parula(numColors);  % Generate a custom colormap with distinct colors
%Legend set-up
LegNames = cell(13,1);
for i = 1:size(h_corr, 2) %plot h_corr for each cycle
    scatter(along_track_distance_new, h_corr(:, i), 19, colors(i, :), 'filled');
    dateInd = find(~isnan(utcDateNum(:,i)), 1); %find first non-NaN instance of DateNum in column i  for inclusion in Legend
    LegName = ['Cycle ', num2str(cycle_number(i)), ': ', datestr(utcDateTime(dateInd, i), 'dd mmm yyyy')];
    LegNames{i} = LegName;
end
scatter(along_track_distance_new, roi_dem, 10, 'k', 'filled'); %Plot REMA-extracted elevations provided as an ATL11 attribute
hold off;
LegendText = [LegNames; 'REMA elevation'];
legend(LegendText, 'Location', 'northeast');
leg = legend('show');
leg.FontSize = 14; 
xlabel('Along-track distance (km)','FontSize', 16,'FontWeight', 'bold');
ylabel('Elevation (m)','FontSize', 16,'FontWeight', 'bold');
hAx = gca;
hAx.FontSize = 14; 
hAx.Box = 'on';
hAx.LineWidth = 1;
hAx.TickLength = [0.005, 0.1];
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
grid on;
title('ICESat-2 and REMA Along-track elevations');


%% Plot ICESat-2 along-track elevation profile for all cycles alongside REMA along-track elevations and BedMachine bed elevations
figure;
hold on;
yyaxis left; % Set the primary y-axis
numColors = size(h_corr, 2); %calculate number of columns in h_corr
colors = parula(numColors);  % Generate a custom colormap with distinct colors
LegNames = cell(13,1);
for i = 1:size(h_corr, 2) %plot h_corr for each cycle
    scatter(along_track_distance_new, h_corr(:, i), 19, colors(i, :), 'filled');
    dateInd = find(~isnan(utcDateNum(:,i)), 1); %find first non-NaN instance of DateNum in column i  for inclusion in Legend
    LegName = ['Cycle ', num2str(cycle_number(i)), ': ', datestr(utcDateTime(dateInd, i), 'dd mmm yyyy')];
    LegNames{i} = LegName;
end
scatter(along_track_distance_new, roi_dem, 10, 'k', 'filled');
xlabel('Along-track distance (km)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Ice Surface Elevation (m)', 'FontSize', 16, 'FontWeight', 'bold');
hAx = gca;
hAx.FontSize = 18; 
hAx.Box = 'on';
hAx.LineWidth = 1;
hAx.TickLength = [0.005, 0.1];
set(gca, 'FontSize', 18, 'FontWeight', 'bold');
grid on;
set(gca, 'Ylim', [450, 560]);
set(gca, 'Xlim', [280, 310]);
yyaxis right; % Set the secondary y-axis
BedMachine = 'G:\BedMachine\bedelev.tif';
BedElev = raster2point_extract(selectedE, selectedN, BedMachine, 'nearest');
plot(along_track_distance_new, BedElev, 'k', 'LineWidth', 2);
ylabel('BedMachine Bed Elevation (m)', 'FontSize', 16, 'FontWeight', 'bold');
hAx.YAxis(2).Color = 'k'; % Set secondary y-axis label color to black
%hAx.YAxis(2).TickLabelFormat.Color = 'k'; 
hAx = gca;
hAx.FontSize = 16; 
hAx.Box = 'on';
hAx.LineWidth = 1;
hAx.TickLength = [0.005, 0.1];
set(gca, 'FontSize', 16, 'FontWeight', 'bold');
grid on;
% ylim_right = [min(BedElev), max(BedElev)];
% ylim(ylim_right);
legend(LegNames,'Location', 'northeast');
leg = legend('show');
leg.FontSize = 18; 
hold off;


%% Plot Along-Track Elevation Anomalies
figure;
hold on;
numColors = size(h_corr, 2); %calculate number of columns in h_corr
colors = parula(numColors);  % Generate a custom colormap with distinct colors
%Legend set-up
LegNames = cell(13,1);
handles = [];
for i = 1:size(h_corr, 2) %plot h_corr for each cycle
    h = plot(along_track_distance_new, running_mean(:, i), 'Color', colors(i, :), 'LineWidth', 2);
    handles = [handles; h]; % Store the handle for this line
    dateInd = find(~isnan(utcDateNum(:,i)), 1); %find first non-NaN instance of DateNum in column i  for inclusion in Legend
    LegName = ['Cycle ', num2str(cycle_number(i)), ': ', datestr(utcDateTime(dateInd, i), 'dd mmm yyyy')];
    LegNames{i} = LegName;
end
hold off;
LegendText = [LegNames];
legend(handles, LegendText, 'Location', 'south', 'NumColumns', 3, 'FontSize', 18); 
leg = legend('show');
leg.FontSize = 20; 
xlabel('Along-track distance (km)','FontSize', 20,'FontWeight', 'bold');
ylabel('Elevation Anomaly (m)','FontSize', 20,'FontWeight', 'bold');
line(get(gca, 'XLim'), [0 0], 'Color', 'k');
hAx = gca;
hAx.FontSize = 20; 
hAx.Box = 'on';
hAx.LineWidth = 1;
hAx.TickLength = [0.005, 0.1];
set(gca, 'FontSize', 20, 'FontWeight', 'bold','FontName', 'Arial');
grid on;
title('ICESat-2 ATL11 along-track elevation anomalies Lake M2 065512');

%% Calculate lake volume estimate for each ICESat-2 cycle within input lake boundary shapefile based on median elev change within lake boundary
lakeBoundary=shaperead('C:\ExampleLakeShapefile.shp'); %Update filepath accordingly
lakeBoundaryGeometry=lakeBoundary.Geometry;
lakeBoundaryX = lakeBoundary.X;
lakeBoundaryY = lakeBoundary.Y;
num_cycles = size(h_corr, 2);
intersectsLake = inpolygon(ROIdata(:, 1), ROIdata(:, 2), lakeBoundaryX, lakeBoundaryY); %check which ICESat-2 points along track intersect with REMA-derived lake boundary
num_points = sum(intersectsLake); % Number of ICESAt-2 track points inside the lake boundary
elevation_anomaly_inside_lake = NaN(num_points, num_cycles); % initialise matrix
elevation_anomaly_inside_lake(:, 1:num_cycles) = elevation_anomaly2(intersectsLake, :); %stores elev anomalies inside lake boundary for each cycle
LakeElevAnom = nanmedian(elevation_anomaly_inside_lake, 1); %calculate all median elev anomalies within lake polygon
LakeElevAnom(isnan(LakeElevAnom)) = 0; % replace any NaN values with zero to exclude them from lake volume calculations
%Convert anomalies to km and multiply by REMA-derived lake area in km2, to calculate lake volume for each cycle
LakeVol = (abs(LakeElevAnom)/1000) * 39.4734; % n.b remember to change Lake area accordingly. Calculated within GIS from user-provided lake boundary shapefile.

