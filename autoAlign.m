function [centers_srt,faults,absents,radii_srt,meanRadius,tries,circleSensitivity] = autoAlign(refimg,nTrays,radiiRange,EdgeThreshold,nWells)
%function for automatic well alignment

%% change log


% extracted from autowell2 function July 2019 by Killian Brennan
% still returns if the wrong number of wells where detected to allow for subsequent manual alignment

% originally addapted from well_compare function by R.O.D Jan 2018

tries = 1;
found = 0;
smaller = 0;
bigger = 0;
searchSize = 20; % ends search after this many attempts
circleSensitivity = 0.95670; % .96 starting value
dSensitivity = .02; % initial search gap size

% converge on correct sensitivity
while tries < searchSize+1 && found == 0
    if tries == 1
        fprintf('Trying circle sensitivity: %.3f \n', circleSensitivity);
    end
    
    [centers,radii]=imfindcircles(refimg,radiiRange,'ObjectPolarity','bright',...
        'Sensitivity',circleSensitivity,'EdgeThreshold',EdgeThreshold); %,'Method','TwoStage'); %finds the wells
    
    fprintf('Number of wells detected: %d \n', length(radii));
    if length(radii) == nWells
        found = 1;
        
    elseif length(radii) < nWells
        smaller = 1;
        if bigger == 1
            dSensitivity = dSensitivity/2;
            bigger = 0;
        end
        circleSensitivity = circleSensitivity + dSensitivity;
        fprintf('Increasing circle sensitivity to: %.3f \n', circleSensitivity);
        
    elseif length(radii) > nWells
        bigger = 1;
        if smaller == 1
            dSensitivity = dSensitivity/2;
            smaller = 0;
        end
        circleSensitivity = circleSensitivity - dSensitivity;
        fprintf('Decreasing circle sensitivity to: %.6f \n', circleSensitivity);
    end
    
    tries = tries+1;
end


% sort the centers so that they are row by row
z = centers(:,1)./(2.*mean(radii)).*(max(centers(:,2))-min(centers(:,2)))+centers(:,2);
[~,srt_ord] = sortrows(z);
centers_srt = centers(srt_ord,:);
radii_srt = radii(srt_ord);
meanRadius = mean(radii);


% spacial search to fix missdetected wells

wellsCenter(1) = mean(centers_srt(:,1));
wellsCenter(2) = mean(centers_srt(:,2));

distFromCenter(:,1) = (centers_srt(:,1) - wellsCenter(1));
distFromCenter(:,2) = (centers_srt(:,2) - wellsCenter(2));

distFromCenterSorted(:,1) = sort((distFromCenter(:,1)));
distFromCenterSorted(:,2) = sort((distFromCenter(:,2)));

centersSize(1) = abs(distFromCenterSorted(end-8,1))+abs(distFromCenterSorted(8,1)); % exclude half of the outer rows (and outliers)
centersSize(2) = abs(distFromCenterSorted(end-nTrays*3,2))+abs(distFromCenterSorted(nTrays*3,2)); % exclude half of the outer rows (and outliers)

gridSpacing = centersSize(1)/nTrays/6;

gridDat(:,:,1) = [repelem(linspace(1,centersSize(1),nTrays*6)',16),repmat(linspace(1,centersSize(2),16)',[nTrays*6,1])];

gridDat(:,1,1) =  gridDat(:,1,1)+mean(wellsCenter(1))-centersSize(1)/2;
gridDat(:,2,1) =  gridDat(:,2,1)+mean(wellsCenter(2))-centersSize(2)/2;

shift = [0,1,-1,2,-2];

for iii = 1:5
    gridDat(:,:,iii) = gridDat(:,:,1);
    
    gridDat(:,1,iii) =  gridDat(:,1,iii)+gridSpacing*shift(iii);
    
    dt = delaunayTriangulation(gridDat(:,:,iii));
    gridIndex(:,iii) = nearestNeighbor(dt, centers_srt);
    
    xnn(:,:,iii) =  gridDat(gridIndex(:,iii),:,iii);
    
    distanceX(:,iii) = centers_srt(:,1)-(xnn(:,1,iii));
    distanceY(:,iii) = centers_srt(:,2)-(xnn(:,2,iii));
    
    quantileDistanceX(iii) = quantile(abs(distanceX(:,iii)),0.8);
    quantileDistanceY(iii) = quantile(abs(distanceY(:,iii)),0.8);
    
    distanceXtrim(:,iii) = distanceX(abs(distanceX(:,iii))<quantileDistanceX(iii),iii);
    distanceYtrim(:,iii) = distanceY(abs(distanceY(:,iii))<quantileDistanceY(iii),iii);
    
    gridDat(:,1,iii) =  gridDat(:,1,iii)+mean(distanceXtrim(:,iii));
    gridDat(:,2,iii) =  gridDat(:,2,iii)+mean(distanceYtrim(:,iii));
    
    distances(:,iii) = sqrt((xnn(:,1,iii)-centers_srt(:,1)).^2 + (xnn(:,2,iii)-centers_srt(:,2)).^2);
        
    sortedDistances(:,iii) = sort(distances(:,iii));
    meanDistance(iii) = mean(sortedDistances(80*nTrays:end-nTrays,iii)); % ignore closest and furthest wells
    
end
[minDist,minInd] = min(meanDistance);
gridDatFinal = gridDat(:,:,minInd);
faults = centers_srt(distances(:,minInd)>mean(radii)*1.2,:);

absents =  gridDatFinal(~ismember(1:96*nTrays,gridIndex(:,minInd)),:);

centers_srtOld = centers_srt;
centers_srt =  gridDatFinal; % overwrite original centers coordinates
end