function [] = autowell(foldername,varargin)

% No further functionality is to be added to this code that would not
% be usable to all users. These options must be implemented outside of
% 'autowell'.

% autowell requires the following functions to opperate:
% autoAlign
% manualAlign
%

% Function to loop through all of the image files in a folder and calculate
% when the wells froze.  The code automatically identifies the wells and
% then looks at the change in light intensity between images to determine
% when the well changes from liquid to ice. A set threshold intensity can
% be set in order to correctly account for wells that do not freeze
% completely.

% will add Frz_T vector and plots to folder where the DRINCZ
% data is stored

% prefeered folder format: yyyymmdd_IDID but you can choose any


% % I/O
% possible example inputformats:
% autowell('20190228_test3')  (runs autowell with one tray)
% autowell('20190228_test3','nTrays',2)  (runs autowell with two trays)
% autowell('20190228_test3',2)  (runs autowell with two trays)
% autowell('20190228_test3','nTrays',2,'names',{'sample1','sample2'})
% autowell('20190228_test3','manual',1)  (runs autowell directly in manual modus)

% Input
% foldername = folder name for analysis needs to be entered as a string '...'
% further name pair arguments:
% 'peak_thresh' optional input threshold for a well change to be considered
% as a freezing
% event.  This helps ensure that wells that don't freeze are not
% misclassified as frozen and is optional (.6 is the default)
% 'calibration','on' or 'off'(default)
% 'names',{'tray1','tray2','tray3','tray4'}(default)
% 'nTrays', 1(default number of trays, max:4)

% % change log
% updated July 2019 by Killian Brennan:
% detect wrong well alignment and notify user
% manual option for well alignment (manualAlign)
% outputs a metadata logfile for trouble shooting and traceability
% outputs matlab workspace for trouble shooting and advanced reanalysis
% added combined output structure
% removed automatic well alignment to external function (autoAlign)
% added option to skip images during analysis of freezing (runs faster but lower temperature resolution)
% added figure for freezing evaluation to output
% changed exposure normalisazion reference (only over trays)

% updated Mar 2019 by Killian Brennan:
% added spacial search to find missdetected wells

% updated Feb 2019 by Killian Brennan:
% adapted to DRINC-384
% added multitray support
% streamlined data processing and output
% normalizes exposure according to 99.9 precentile of pixles (brightest pixles)

% updated Jul 2018 by Killian Brennan:
% included optional calibration of temperature

% updated Jun 2018 by Killian Brennan:
% decreased well mask radius (0.7) to avoid "waviness"

% updated Apr 2018 by Killian Brennan:
% generalized function for use by all users
% more efficient search algorithmus for determining correct circleSensitivity
% specified EdgeThreshold to solidify circle detection and avoid false detection

% updated Feb 2018 by Killian Brennan:
% advanced user feedback
% error messages for undetected wells
% added mac compatibility
% swiched from ls() to dir() for folder exploration

% addapted from well_compare function by R.O.D Jan 2018

% + magic_args="parameters"
subsampling = .5; % subsampling factor (rescales image to improve processing time)
skip = 0; % number of images skipped for analysis
radiusCorrection = 1; % modulates detected radius before generating mask
radiiRange = [30 40]*subsampling; % select radius size range (in pixels)
EdgeThreshold = 0.01; % edge threshold for well detection
traySize = 96; % number of wells per tray
default_nTrays = 3; % default number of trays
default_peak_thresh=0.8; % Automatically set the ice threshold if it is not is provided
default_calibration = 'off'; % turn calibration on or off
% -

colors = {'r','g','b','c'}; % plotting colors

auto_available = 0; % flag for logfile output

% calibration (optional)
intersect = 1.3;
slope = 0.917;

% + magic_args="close all open figures"
close all

% + magic_args="save parent directory"
parentDir = cd;

% + magic_args="parse imputs"
p = inputParser;
% -

addRequired(p,'foldername');

% set default number of trays
nTrays = default_nTrays;
addOptional(p,'nTrays',default_nTrays,@isnumeric);

% Set name for structure if none is provided
% folder name must be space free and a valid varriable name

% Set the ice threshold if none is provided
addParamValue(p,'peak_thresh',default_peak_thresh,@isnumeric);

% Set the calibration option if none is provided
checkString = @(s) any(strcmp(s,{'off','on'}));
addParamValue(p,'calibration',default_calibration,checkString);

% Name samples
addParamValue(p,'names',{'tray1','tray2','tray3','tray4'});

addParamValue(p,'manual',0,@isnumeric);

parse(p,foldername,varargin{:});
peak_thresh = p.Results.peak_thresh;
calibration  = p.Results.calibration;
nTrays  = p.Results.nTrays;
names = p.Results.names;
manual = p.Results.manual;

nWells = nTrays*traySize; % number of wells being analyzed

% + magic_args="Identify the well location (read in the first image) and retrieve the well change over the run"


cd(foldername) % must be in the correct folder or use the whole path here i.e. cd(strcat('D\DropFreezing\BaselTest\',foldername))
dirstruct = dir('*.jpg'); % find info about .jpg's in folder
pcrnmAll = char(extractfield(dirstruct,'name')'); % extract name information from structure & transform for correct format & convert to char vector
pcrnm = pcrnmAll(1:skip+1:size(pcrnmAll,1),:); % only use every skipth image for analysis
for i=1:size(pcrnm,1)
    if i==1 % analyze first image
        refimg=imresize(imread(pcrnm(i,:)),subsampling);
        pic_size=size(refimg);
        [rr,cc]=meshgrid(1:pic_size(2),1:pic_size(1));
        
        if manual == 0
            auto_available = 1;
            [centers_srt,faults,absents,radii_srt,meanRadius,tries,circleSensitivity] =...
                autoAlign(refimg,nTrays,radiiRange,EdgeThreshold,nWells);
            
            % Creates a figure to ensure all of the wells are correctlyidentified (or not)
            figure;
            imshow(refimg);
            
            for ii = 1:nTrays
                viscircles(centers_srt(1+traySize*(ii-1):traySize*ii,:),...
                    repmat(meanRadius*radiusCorrection,...
                    [length(centers_srt(1+traySize*(ii-1):traySize*ii,:)),1]),...
                    'Color',colors{ii});
                text(centers_srt(1+traySize*(ii-1),1),...
                    centers_srt(1+traySize*(ii-1),2)-2.5*radii_srt(1+traySize*(ii-1)),...
                    names{ii},'Color',colors{ii},'FontSize',14);
            end
            
            hold on
            scatter(faults(:,1),faults(:,2),500,'sy')
            scatter(absents(:,1),absents(:,2),500,'sk')
            %plot([xnn(:,1,minInd) centers_srtOld(:,1)]',[xnn(:,2,minInd) centers_srtOld(:,2)]','-r');
            print([foldername,'_Mask_Auto'],'-dpng','-r300');
            
            % proceed with manual well selection in case too many faults
            if size(faults,1)>=24
                manual = 1;
            end
        end
        
        % manual well selection
        if manual == 1
            disp('ERROR: wells could not be automatically detected');
            disp('please manually select upper left and lower right wells in GUI');
            [centers_srt,meanRadius,radii_srt] = manualAlign(refimg,nTrays,traySize); % overwrite automatic with manual alignment
            
            figure;
            imshow(refimg);
            
            for ii = 1:nTrays
                viscircles(centers_srt(1+traySize*(ii-1):traySize*ii,:),...
                    repmat(meanRadius*radiusCorrection,...
                    [length(centers_srt(1+traySize*(ii-1):traySize*ii,:)),1]),...
                    'Color',colors{ii});
                text(centers_srt(1+traySize*(ii-1),1),...
                    centers_srt(1+traySize*(ii-1),2)-2.5*radii_srt(1+traySize*(ii-1)),...
                    names{ii},'Color',colors{ii},'FontSize',14);
            end
            
            print([foldername,'_Mask_Manual'],'-dpng','-r300');
            
        end
        
      
        wellsCenter(1) = mean(centers_srt(:,1));
        wellsCenter(2) = mean(centers_srt(:,2));
        gridSpacing = abs(mean([centers_srt(1,1)-centers_srt(17,1),...
            centers_srt(1,2)-centers_srt(2,2)]));
        
        
        mask = zeros(pic_size(1),pic_size(2)); % create mask blank
        
        
        disp('Creating a mask for the well locations.'); %Give the user some feedback
        
        h = waitbar(0,'Creating well masks.'); % show progress in waitbar
        for j = 1:length(radii_srt)
            mask = mask + j*(sqrt((rr-centers_srt(j,1)).^2+...
                (cc-centers_srt(j,2)).^2)<=radiusCorrection*meanRadius); %Create a mask of the wells.
            waitbar(j / length(radii_srt));
        end
        close(h); % close waitbar
        
        %refimg = im2uint8(im2double(refimg) / quantile(im2double(refimg(:)),0.999)); % normalize exposure
        refimg = im2uint8(im2double(refimg) /...
            mean(mean(im2double(refimg(round(centers_srt(1,2)):round(centers_srt(end,2)),...
            round(centers_srt(1,1)):round(centers_srt(end,1))))))); % normalize exposure

        for k = 1:length(radii_srt)
            wells(i,k) = mean(refimg(mask == k)); %get the mean pixel value for each well for the first image
        end
        disp('First image processed, looping through the rest of the pictures.') %Give the user some feedback
    else
        if i == 2
            h = waitbar(0,'Analyzing following images.'); % show progress in waitbar
        end
        waitbar(i / length(pcrnm));
        
        img = imresize(imread(pcrnm(i,:)),subsampling);
        %img = im2uint8(im2double(img) / quantile(im2double(img(:)),0.999)); % normalize exposure
        img = im2uint8(im2double(img) /...
            mean(mean(im2double(img(round(centers_srt(1,2)):round(centers_srt(end,2)),...
            round(centers_srt(1,1)):round(centers_srt(end,1))))))); % normalize exposure

        for k = 1:length(radii_srt)
            wells(i,k) = mean(img(mask == k)); %get the mean pixel value for the remaining pictures
        end
    end
end
close(h); % close waitbar
% -

% % Determine when freezing occured and rank the wells when they froze
% normalize the well change and look for the first change where the value is
% larger than the peak_thresh of the max change of the well

Frz_Index = nan(1,nWells);
for iFrz = 1:nWells
    
    % ROD method
    Frz_Index(1,iFrz) = find(abs(diff((wells(:,iFrz)-mean(wells)./std(wells))))./...
        max(abs(diff((wells(:,iFrz)-mean(wells)./std(wells)))))>=peak_thresh,1,'first')+1;
    
    % alternate method
    [maxima Frz_Index_New(1,iFrz)] = max(-diff(wells(:,iFrz))); 
    Frz_Index_New(1,iFrz) = Frz_Index_New(1,iFrz)+1;
end

[~,~,Frz_Rank] = unique(Frz_Index');

Frz_name = pcrnm(Frz_Index,:);
Frz_time = Frz_name(:,1:end-4);

Scan_time = pcrnm(:,1:end-4);

% + magic_args="Calculate Frozen Fraction for each time step"
nFrz = zeros(size(wells,1),1);
for iFF = 1:size(wells,1)
    if iFF == 1
        nFrz(iFF) = sum(Frz_Index == iFF);
    else
        nFrz(iFF) = sum(Frz_Index == iFF)+nFrz(iFF-1);
    end
end


% + magic_args="Find the temperature data"
dirstruc = dir('*.txt');
tname = dirstruc.name;
if isempty(tname)
    error('ERROR: No temperature data available!');
    return;
else
    ramp = readtable('ramp.txt'); % must be ramp.txt
    rawtstmp = datestr(table2array(ramp(:,1)));
    tstmp = datetime(rawtstmp(:,13:20),'Inputformat','HH:mm:ss'); % times corresponding to logged temperatures
    Scan_date = datetime(Scan_time(:,:),'Inputformat','HH_mm_ss'); %convert the time into the correct freezing time
    Frz_date = datetime(Frz_time(:,:),'Inputformat','HH_mm_ss'); %convert the time into the correct freezing time
    
    tdata = table2array(ramp(:,2)); % logged temperatures
    save(strcat(foldername),'tdata'); % save logged temperatures outside workspace
    
    
    for dt = 1:length(Frz_date)
        [Frz_t_uncert(1,dt),Frz_T_index(1,dt)] = min(abs(tstmp-Frz_date(dt))); %extract the index temperature for when each well froze
    end
    
    for dt2=1:size(pcrnm,1)
        [t_uncert(1,dt2),T_index(1,dt2)] = min(abs(tstmp-Scan_date(dt2))); %extract the index temperature for each picture
    end
    
    FrzTall=tdata(Frz_T_index); %Temperature when each well froze
    
    figure;
    
    heatmap(reshape(FrzTall,[16,6*nTrays]),'colormap',redbluecmap); %add max(Frz_Rank) if needed also could use flipud(redbluecmap)
    
    set(gca,'fontsize',12);
    print([foldername,'_frzOrder'],'-dpng','-r300');
% -

end

% plot simple figures
figure;
hold on
for ii = 1:nTrays
    plot(sort(FrzTall(1+traySize*(ii-1):traySize*ii)),...
        linspace(1,0,traySize),'Color',colors{ii});
end

legend(names);
xlabel('Temperature ( ^\circC)');
ylabel('Frozen fraction');
grid on;
set(gca,'fontsize',12);
print(foldername,'-dpng','-r300');

figure; % freezing evaluation figure
hold on
image(-diff(wells)')
scatter(Frz_Index,1:nWells,'xw')
ylabel('well nr')
xticks(1:numel(T_index));
xticklabels(pcrnm);
set(gca, 'TickLabelInterpreter', 'none')
xtickangle(90)
c = colorbar;
c.Label.String = 'diff(intensity)';
for i = 1:nTrays-1
%yline(96*i,'--w','lineWidth',1)
end
ylim([1 nTrays*96])
xlim([min(Frz_Index)-1 max(Frz_Index)+1])

savefig([foldername,'_freezing.fig']);


if strcmp(calibration,'on')
    FrzTall = intersect+FrzTall*slope;
end

% save the data to the current folder
for ii = 1:nTrays
    FrzT = FrzTall(1+traySize*(ii-1):traySize*ii);
    save(strcat(foldername,'_',names{ii}),'FrzT');
end
save(strcat(foldername),'FrzTall');

save('workspace') % save workspace for trouble shooting

% save metadata log file

fid=fopen('logfile.txt','w');
fprintf(fid, [ 'logfile for:                  ' foldername '\n']);
fprintf(fid, [ 'autowell run on:              ' datestr(now) '\n']);
fprintf(fid, [ 'autowell run in directory:    ' cd '\n']);
fprintf(fid, ['\n'  '**********   general info   ***********' '\n' ]);

fprintf(fid, [ 'TFrz output calibration was   ' calibration '\n']);
fprintf(fid, [ 'number of trays:              ' num2str(nTrays) '\n']);
fprintf(fid, [ 'number of wells:              ' num2str(nWells) '\n']);
fprintf(fid, [ 'number of images captured:    ' num2str(size(pcrnm,1)) '\n']);
fprintf(fid, [ 'image subsampling:            ' num2str(subsampling) '\n']);
fprintf(fid, [ 'number of images skipped:     ' num2str(skip) '\n']);


fprintf(fid, ['\n'  '**********   well alignment   ***********' '\n' ]);
if manual == 1 & auto_available == 1
    fprintf(fid, 'automatic well alignment failed, manually aligned by user\n');
elseif manual ==1 & auto_available == 0
    fprintf(fid, 'wells where manually aligned by user\n');
elseif manual == 0
    fprintf(fid, 'wells where automatically aligned\n');
end

if auto_available == 1 % only write this if auto alignment was performed
    fprintf(fid, [ 'number of wells not found:    ' num2str(size(absents,1)) '\n']);
    fprintf(fid, [ 'number of wells missdetected: ' num2str(size(faults,1)) '\n']);
    fprintf(fid, [ 'circle detection sensitivity: ' num2str(circleSensitivity) '\n']);
    fprintf(fid, [ 'tries to reach correct count: ' num2str(tries) '\n']);
end

fprintf(fid, [ 'grid spacing:                 ' num2str(gridSpacing) '\n']);
fprintf(fid, [ 'mean radius:                  ' num2str(meanRadius) '\n']);
fprintf(fid, [ 'wellsCenterX:                 ' num2str(wellsCenter(1)) '\n']);
fprintf(fid, [ 'wellsCenterY:                 ' num2str(wellsCenter(2)) '\n']);
fprintf(fid, ['\n'  '**********   end of logfile   ***********' '\n' ]);

%fprintf(fid, '%f %f r\n', [A B]');
fclose(fid);

% return to parent directory
cd(parentDir);
disp('autowell complete');
end
