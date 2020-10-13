function [centers_srt,meanRadius,radii_srt] = manualAlign(refimg,nTrays,traySize)
%function for manual well alignment (UI)

%% change log
% written July 2019 by Killian Brennan

figure;
imshow(refimg);
title('select center of upper left well','fontsize',14);
disp('select center of upper left well');
[ulx,uly] = ginput(1); % allow user to select coordinates of upper left well
title('select center of lower right well','fontsize',14);
disp('select center of lower right well');
[lrx,lry] = ginput(1); % allow user to select coordinates of lower right well

x_s = linspace(ulx,lrx,6*nTrays)';
y_s = linspace(uly,lry,16)';

centers_srt = [repelem(x_s,16),repmat(y_s,nTrays*6,1)];

meanRadius = abs(mean([centers_srt(1,1)-centers_srt(17,1),...
    centers_srt(1,2)-centers_srt(2,2)])*0.33); % aproximate radii from well spacing
radii_srt = repelem(meanRadius,nTrays*traySize)';
close all;
end
