function [] = FINCrunner(measName,varargin)
% possible example inputformats:
% FINCrunner('test3')  (runs autowell with 3 trays)
% FINCrunner('test3','nTrays',2)  (runs autowell with two trays)
% FINCrunner('test3','nTrays',2,'names',{'sample1','sample2'})
% adds current date to foldername automatically 

%% I/O
% measName = string of measurement name for folder generation
% further name pair arguments passed on to autowell: (varargin -> variable argument input
% 'peak_thresh' (ice detection threashold) 
% 'calibration','on' or 'off'(default)
% 'names',{'sample1','sample2','sample3','sample4'}(default)
% 'nTrays', 1(default number of trays, max:4) [but only actually use 3]

%% change log
% updated July 2019 by Killian Brennan:
% changed triggering of pictures and temperature from time to temperature intervall

% written Feb 2019 by Killian Brennan

% creates a measurement specific folder according to 'foldername'
% run in global measurement folder
foldername = strcat(datestr(today,'yyyymmdd'),'_',measName);%generates folder name
oldDIR = cd; %saves this as the current directory
mkdir(foldername); %making a new folder of the foldername
cd(foldername); %putting yourself in the current directory of the new folder

close all; % closes all open figures
% close all open waitbars
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%parameters
tempEnd = -30; % end temperature in degC - bath temperature not set temperature - don't go lower than -32
tempStart = 0; % start temperature in degC
ramp = -1; % ramp speed in degC/min
measInt = .2; % measuring intervall in degC - how often the data is recorded (photos and ramp)

%% Initialize bath
COMport = seriallist; %looks at all the connected ports 
if ~strcmp(COMport,'COM5') %COM5 is where the LAUDA is normally connected
    error('Instrument not connected or wrong COM port issue!')
end

delete(instrfindall); % close serial port if it was open - closes any previously open conversations

LAUDA = serial('COM5'); %naming the connection
set(LAUDA,'BaudRate',9600); %frequency of the conversation, just has to be the same, find in LAUDO manual
LAUDA.Terminator = 'CR/LF'; %"carriage return / line feed" -- what type of communication between matlab and lauda- what type of break
fopen(LAUDA); %opens the communication channel, which is on COM5

pause(.1) 
fprintf(LAUDA,'START'); %turns on Lauda, tells Lauda to start. these commands are in Lauda manual
status = fgetl(LAUDA); %gets the infomration from Lauda 
if strcmp(status,'OK') %lauda will respond OK if it turns on
    disp('LAUDA is turned on');
else
    error('LAUDA could not be turned on');
end
pause(.1)
fprintf(LAUDA,strcat('OUT_SP_00_',num2str(tempStart))); %the command for Lauda to go to the start temperature (also in Lauda manual), and we define what the tempstart is (line 36)
status = fgetl(LAUDA); 
if ~strcmp(status,'OK')
    error('LAUDA not responding');
end

%% Initialize camera
% camera interface adapted from R.O.D.'s Feb. 2018 camera_pic.m & cam_loader.m; modified while loop from: https://ch.mathworks.com/matlabcentral/answers/261417-matlab-webcam-inconsistencies
imaqreset; % reset image aquisition
opengl software %Tells matlab to use opengl software for graphics, this should fix the crashing problem

cam_name = 'HD USB Camera';  % assigns the webcam name
clear(cam_name); % clear existing camera
maxAttempts = 10; % sets the number of attempts to obtain the properties
attempt = 1; % initialize the attempt number
hasProperty = false; % initialize the lack of properties
while attempt < maxAttempts  % start the while loop, load cam until properties are there
    cam = webcam(cam_name); % load the camera
    hasProperty = max(ismember(properties(cam), 'Contrast')); % look for the Contrast property, thereby looking for all properties
    if hasProperty
        % We got all the properties so we can quit trying.
        fprintf('We have all camera properties after %d attempts.  Now continuing with program.\n', attempt);
        break;
    end
    fprintf('We do not have all camera properties after %d attempts.  Trying again.\n', attempt);
    attempt = attempt + 1;
    clear(cam_name);
end
if attempt == maxAttempts
    error('Camera not found!')
end
cam.Resolution = '3264x2448'; %depends on the camera
cam.ExposureMode='manual'; %vensures the exposure doesn't change
cam.Exposure=-8; % This may need to be adjusted once the camera is mounted
cam.WhiteBalanceMode='manual';
cam.WhiteBalance=4600; %normally 4600
cam.Saturation = 0; % makes image Monochrome

preview(cam) %starts the camera preview
pause(10);

%% UI stop measurement
%ButtonHandle = uicontrol(cam,'Style', 'PushButton', ...
    %'String', 'Stop Measurement', ...
   % 'Callback', 'delete(gcbf)');

%% Measurement
% wait for temperature to reach 0
disp('Bath temperature initializing') %writes into the command window
pause(.1);
fprintf(LAUDA,strcat('OUT_SP_00_',num2str(tempStart)));
status = fgetl(LAUDA);
if ~strcmp(status,'OK')
    error('LAUDA com error');
end
pause(.1);
fprintf(LAUDA,'IN_PV_00'); %gets the measured bath temperature from LAUDA
tempBath = str2num(fgetl(LAUDA)); %converts that bath temp to a numeric 
tempOriginal = tempBath; % initial starting temperature for waitbar calculation

h = waitbar(0,'Current bath temperature is ','Name','Bath temperature initializing...'); % show progress in waitbar

while round(tempBath,1) ~= tempStart % wait for temp to get within 0.05C of tempStart
    pause(.1);
    fprintf(LAUDA,'IN_PV_00');
    tempBath = str2num(fgetl(LAUDA));
    waitbar(min(max((tempOriginal-tempBath)/(tempOriginal-tempStart),0),1),h,...
        ['Current bath temperature is ',num2str(tempBath),'C'],'Name','Bath temperature initializing...'); % update waitbar
    pause(1);
end
close(h); %close the waitbar

disp('Instrument ready, starting measurement.')
h = waitbar(0,'Current bath temperature is ','Name','Measurement in progress...'); % show progress in waitbar

tempRamp = cell2table(cell((tempStart-tempEnd)*60/(-ramp)/6+10,2),'VariableNames',{'time' 'temp_degC'}); % initiate output "ramp" table where it will write to
timeStart = now*3600*24; % start time in s (now is the current day since something)
i = 1;
breakinner = false; %use this to break a loop from the inside to go to the next loop
while tempBath>tempEnd % measure until tempEnd is reached
    %while round(mod(now*3600*24-timeStart,1),1)~=0 % update set temp and read it every 1s
        pause(.01);
        waitbar(min(max((tempStart-tempBath)/(tempStart-tempEnd),0),1),h,...
            ['Current bath temperature is ',num2str(tempBath),'C'],'Name','Measurement in progress...'); % update waitbar
%        if ~ishandle(ButtonHandle)
%             disp('Measurement stopped by user');
%             breakinner = true;
%             break;
%         end
%         if(breakinner)
%             break;
%         end
    %end
%     if(breakinner)
%         break;
%     end
    
    fixedNow = now; %save the current time of the measurement so that we can later write it to the file
    timeElapsed = fixedNow*3600*24-timeStart;
    tempSet = round(tempStart + timeElapsed/60*ramp,2);
    if tempSet > -55 % avoid exceeding instrument temperature set bounds, has to be done once in the loop
        fprintf(LAUDA,strcat('OUT_SP_00_',num2str(tempSet))); % set lauda temp
    end
    status = fgetl(LAUDA);
    if ~strcmp(status,'OK')
        error('LAUDA com error');
    end
    pause(.01);
    fprintf(LAUDA,'IN_PV_00'); %get the current bath temperature from LAUDA
    tempBath = str2num(fgetl(LAUDA));
    
    if tempBath <= -measInt*(i-1) % save image and record temperature every measInt (temperature triggered), if we're at the right spot to measure then we take a photo and leave the 
        refImg = rgb2gray(snapshot(cam));
        refImg = im2uint8(im2double(refImg) / mean(quantile(im2double(refImg(:)),0.999))); % normalize exposure
        imageTitle = [datestr(fixedNow,'HH_MM_SS') '.jpg'];
        imwrite(refImg,imageTitle); %write the image
        tempRamp.time(i) = {datestr(fixedNow,'dd-mmm-yyyy HH:MM:SS')};
        tempRamp.temp_degC(i) = {num2str(tempBath)};
        writetable(tempRamp(1:i,:),'ramp.txt');
        i = i+1;
    end
end %we reach end when we get below the final temperature
close(h); %close the waitbar

pause(.1);
fprintf(LAUDA,strcat('OUT_SP_00_',num2str(tempStart))); %tell Lauda to warm back up
status = fgetl(LAUDA);
if ~strcmp(status,'OK')
    error('LAUDA com error');
end

% close lauda com port
%fprintf(LAUDA,'STOP');
fclose(LAUDA);
delete(LAUDA);
clear LAUDA

% close camera connection
clear(cam_name);

disp('Measurement complete, LAUDA is warming back up.');

% run autowell
cd(oldDIR);
disp('Running autowell');
%autowell(foldername,varargin{:});

end
%below not doing anything- killian's failed attempt
function [] = setBathTemp(temp)
fprintf(LAUDA,strcat('OUT_SP_00_',num2str(temp)));
status = fgetl(LAUDA);
if ~strcmp(status,'OK')
    error('LAUDA com error');
end
end
