clc;
clear;

dirName ='../data/Dynamic/PlottingAllRoute';
dirName_image ='../data/Dynamic/PlottingAllRouteOnMap'; 
namelist = dir(dirName);
len = length(namelist);
colors = zeros(len,3);

% Plot the map
f = figure;
axis tight;
axis([-122.5, -122.0, 37.3, 37.8]);
hold on;
I = imread([dirName_image '/MTV.png']); 
h = image('XData',[-122.5, -122.0],'YData',[37.8,37.3],'CData',I);%note the latitude (y-axis) is flipped in vertical direction
uistack(h,'bottom'); %move the image to the bottom of current stack

% Loop the data file
for i = 3:len
    file_name = namelist(i).name;
    if strcmp(file_name, '.DS_Store')
        continue;
    end
    Data_input = load(dirName + "/" + file_name);
    figure(f);
    plot(Data_input(:,2),Data_input(:,3), 'LineWidth',3);hold on;    
end