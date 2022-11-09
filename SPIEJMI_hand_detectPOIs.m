% detect points of high curvature in an image given the number of critical
% points/max/min

clear; close all; clc;

% load data
SPIEJMI_hand_cleanandbinarizedata

close all;

% input number of POIs
POI_no = 11;

% detect POIs
R_edgecoords = SPIEJMI_hand_findedgecoords(dataR);
R_edgecoords = unique(R_edgecoords,'rows','stable');
% % % % % figure('Name','Binarized dataR'); imshow(dataR,[]); hold on;
tic;
[R_POIs,R_ordered_POIs] = SPIEJMI_hand_min_intangle(R_edgecoords,POI_no);
RLMdetection_toc = toc
T_edgecoords = SPIEJMI_hand_findedgecoords(dataT);
T_edgecoords = unique(T_edgecoords,'rows','stable');
% % % % % figure('Name','Binarized dataT'); imshow(dataT,[]); hold on;
[T_POIs,T_ordered_POIs] = SPIEJMI_hand_min_intangle(T_edgecoords,POI_no);

close all;

POI_indexset = 1:size(T_ordered_POIs,1);
color_arr = rand([POI_no,3]);


iptsetpref('ImshowBorder','tight'); h = figure('Name','Binarized dataR'); %set(h, 'Visible', 'off');
imshow(dataR,[]); hold on;
for i = 1:POI_no
    curr_color = color_arr(i,:);
    plot(R_ordered_POIs(i,1),R_ordered_POIs(i,2),'o','MarkerSize',10,'LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor',curr_color);
    text(R_ordered_POIs(i,1)-20,R_ordered_POIs(i,2)+20,num2str(i),'Color',curr_color,'FontSize',18,'FontWeight','bold');
    pause;
end

iptsetpref('ImshowBorder','tight'); figure('Name','Binarized dataT'); imshow(dataT,[]); hold on;
for i = 1:POI_no
    curr_color = color_arr(i,:);
    plot(T_ordered_POIs(i,1),T_ordered_POIs(i,2),'o','MarkerSize',10,'LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor',curr_color);
    text(T_ordered_POIs(i,1)-20,T_ordered_POIs(i,2)+20,num2str(i),'Color',curr_color,'FontSize',18,'FontWeight','bold');
    pause;
end

% Once the first point in the ordered list of POIs is identified, specify
% the number of sample points to be detected in each subinterval
% Only need to ask once because the number of sample points should be the
% same for the reference and template images