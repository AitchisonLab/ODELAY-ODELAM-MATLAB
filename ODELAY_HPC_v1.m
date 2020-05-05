function ODELAM_HPC_v1(varargin)
% Needs pars data into ODELAM_HPC_v1;
% when the ODELAM_HPC_v1 function is called on the HPC it needs to know the
% directory path that data is stored and that it will write data.
% ODELAM_HPC_v1 should only return errors to the command line and won't
% report errors that occur during data collection.

inArgs.dirLoc = ' ';
inArgs.well  = 0;
inArgs.Init = ' ';
well = 0;

args = fieldnames(inArgs);

if mod(nargin,2) ~= 0
      error('Expecting an even number of arguments');
end


for idx = 1:2:nargin
    
    switch varargin{idx}
        
        case 'well'
            isWellInput = true;
            well = str2double(varargin{idx+1});
            
        case 'InitExp' % This command should initialize processing the experiment in the directory given
            
            currentDir = varargin{idx+1};
%             try 
                cd(currentDir(:));
                [monData, mP] = InitializeExperiment;
%             catch
%                 error(['Directory Location is not correct ', currentDir(:)'])
%             end
            
        case 'dirLoc' % Find the directory where data is and Experiment files should be;
            initArgs.dirLoc = varargin{idx+1};
            try
                cause = ['line 41', initArgs.dirLoc];
                fprintf(initArgs.dirLoc);
                cd(initArgs.dirLoc);
                cause = 'line 45';
                dirInit = true;
                odIndexFile = dir('*Index_ODELAYData.mat');
                numFiles = size(odIndexFile,1);
                cause = 'line 48';
                if numFiles==0
                    error('ODELAY Experiment Not Initialized. Use InitExp command');
                elseif numFiles>1
                    error('ODELAY Experiment too many Index Files');
                else
                    cause = 'line 55';
                    loadData = false;
                    attempts = 0;
                    
                    while ~loadData & attempts<6
                        try
                            load(fullfile(odIndexFile.folder,odIndexFile.name), 'monData');
                            loadData = True;
                        catch
                            pause(2*randn)
                            attempts = attempts+1;
                        end
                    end 
                end
                
            catch 
                error(['Directory Path not found.',cause]);
            end
            
        otherwise
            error('Value "%s" for variable "%s" is not a known command', varargin{idx+1}, varargin{idx});
    end

end

%% Check if everything is working
isMonData = exist('monData')==1;
isWellInput = well<=monData.numWells & well ~= 0;

    if isMonData & isWellInput 
        % check well variable to make sure it is not bigger than the
        % initialized well
        analyzedIm = WellProcess(well, monData);
    elseif well ~= 0
        error('Could not process well because well variable or monData were out of range');
    end

                
end

%% Main Function for Tracking And Analysis
function analyzedIm = WellProcess(well, monData)
%% Process two images and provide tracking data
% Needs files analyzed.  Last list of files, new list of files and the
% number of images to analyze.  

% Load Well Data
oD = LoadDataState(well, monData);

analyzedIm = monData.AnalyzedImages(well,:)==1;
numImAnalyzed = sum(analyzedIm);
% Check how many images are analyzed use this vector to 
% Analyze Data until all current images are found

cd(fullfile(monData.expDir, monData.wellDir{well}));
currentFiles = dir('*.mat');
[~,inds] = sort([currentFiles.datenum]);
currentFiles = currentFiles(inds);

fileInd = 1:numel(currentFiles);
analyzIndex = fileInd(~analyzedIm(1:numel(currentFiles)));

oD.WellData.FileList = currentFiles;


% NewLoadImage +
% LoadOldImage +
% ThresholdOldImage +
% ThresholdNewImage +
% PhaseCorrelate Old+New Evaluate SampleDrift+
% BlobAnalysis+
% Object Track+
% EnterData into ObjectNext and ObjectTrack Data+
% Estimate Growth curves +
% Save Low Bit Depth Image for display +
% Update well analysis +
% Shut down workers once caught up. +
% 
%% Start Processing Data Here
for aI = analyzIndex
            
    imageFile = currentFiles(aI).name;
    % load New Image
    [RawImage, StitchMeta, RawStack,  fluorIm, fT1] = ODELAM_MatReadv5(imageFile,...
                                                                       monData.magnification,...
                                                                       monData.background);
    
    %% Generate a thumbnail of the stitched image for use in the GUI later
    [yD, xD] = size(RawImage);
    oD.WellOut.MovieFrame(aI).cdata = imresize(im2uint8(imadjust(RawImage,[0,max(double(RawImage(:)))/2^16], [])), [105, round(105*xD/yD)]);

    % Use an adaptive background filter to correct for spacial intisity
    % differences
    [Sx,Sy] = imgradientxy(RawStack(:,:,5));
    [sImage, ~] = imgradient(Sx,Sy);

    [Gx,Gy] = imgradientxy(RawImage);
    [AnImage, ~] = imgradient(Gx,Gy);
    clear RawImage Sx Sy Gx Gy 

    %% Binarize Image
    % Threshold the image and record variables that describe the Threshold
    % in terms fo the filtered image histogram
    OffsetMult = 1.6;%ImageVars.OffsetMult;
    ImageHist = imhist(AnImage, double(intmax('uint16')));

    [thrshVal] = Image_Segment(AnImage, ImageHist);

    % fix this bump up an put in Image Segment method and Otsu Method
    [~,maxind] = max(ImageHist(2:end-1));
    bwMaxoffset = abs(thrshVal-maxind)+1;
    bwThresh = round(OffsetMult*bwMaxoffset+maxind);
    bwMinoffset = bwThresh;...-iMinMax(1);
    bwThreshRatio = bwMinoffset;.../(iMinMax(2)-iMinMax(1));
   
    bwImage = AnImage>bwThresh;
    % Binarize the image and store it for later use in linking objects over
    % time
    
    % add these values to ImageVars
    dialateObject = strel('disk',2,0);
    bwImage = imdilate(bwImage, dialateObject);
    bwImage = imfill(bwImage, 'holes');
    dialateObject = strel('disk',2,0);
    bwImage = imerode(bwImage, dialateObject);
    bwImage = bwareaopen(bwImage,9,8);
    
    
    % Try filling in edges of bwImage.  This will remove all objects within
    % edge of image as the same object and set that object to object #1
    bwImage(1,1:end) = true;
    bwImage(end,1:end) = true;
    bwImage(1:end,1) = true;
    bwImage(1:end,end) = true;
    
    
    % Generate a matrix that follows the individual linking lists 
    oD.WellOut.BinaryData.bwThresh(aI,1)      = bwThresh;
    oD.WellOut.BinaryData.bwThreshRatio(aI,1) = bwThreshRatio;
    oD.WellOut.BinaryData.bwMaxoffset(aI,1)   = bwMaxoffset;
    oD.WellOut.BinaryData.bwMinoffset(aI,1)   = bwMinoffset;
    oD.WellOut.BinaryData.bwArea(aI,1)        = double(sum(bwImage(:)));
   
    %% Blob Analysis
    CC = bwconncomp(bwImage, 8);
    CC.NumObjects = CC.NumObjects-1;
    CC.PixelIdxList(1) = [];
    bwLabel = labelmatrix(CC);
    ImageStats = regionprops(CC,...
                             'Centroid',...
                             'FilledArea',...
                             'PixelIdxList',...
                             'PixelList');

    fieldList = fieldnames(ImageStats);
    numObj = size(ImageStats,1);
    % Organize Image stats fields into a reasonable order
    for n = 1:size(fieldList,1)
    numcol = 0;

        switch fieldList{n}
            case 'Centroid'

                xyCent = zeros(numObj,2);
                for m = 1:numObj
                    xyCent(m,:) = ImageStats(m,1).(fieldList{n});
                end
                oD.WellOut.RegionStats(aI).(fieldList{n}) = xyCent;
                numcol = numcol+2;

            case 'FilledArea'

                FilledArea = zeros(numObj,1, 'uint32');
                for m = 1:numObj
                    FilledArea(m,1) = ImageStats(m,1).(fieldList{n});
                end
                oD.WellOut.RegionStats(aI).(fieldList{n}) = FilledArea;
                numcol = numcol+1;

            case 'Eccentricity'

                Eccent = zeros(numObj,1);
                for m = 1:numObj
                    Eccent(m,1) = ImageStats(m,1).(fieldList{n});
                end
                oD.WellOut.RegionStats(aI).(fieldList{n}) = Eccent;
                numcol = numcol+1;

            case 'MajorAxisLength'

                majAxis = zeros(numObj,2);
                for m = 1:numObj
                    majAxis(m) = ImageStats(m,1).(fieldList{n});
                end
                oD.WellOut.RegionStats(aI).(fieldList{n}) = majAxis;
                numcol = numcol+1;

            case 'MinorAxisLength'
                minAxis = zeros(numObj,2);
                for m = 1:numObj
                    minAxis(m) = ImageStats(m,1).(fieldList{n});
                end
                oD.WellOut.RegionStats(aI).(fieldList{n}) = minAxis;
                numcol = numcol+1;
                
        end
    end
    
    %% Extract Fluorescence data from Fluoresences image
    if ~isempty(fluorIm)
       fStates = fieldnames(fluorIm);
       if aI==1
       
           for fI = 1:numel(fStates)
               oD.WellOut.fData.(fStates{fI}).sumIntesity = NaN(5000,size(monData.AnalyzedImages,2));
               oD.WellOut.fData.(fStates{fI}).trackIntesity = NaN(5000,size(monData.AnalyzedImages,2));
               for n = 1:numObj
         
                 vals = fluorIm.(fStates{fI}).stitchImage(CC.PixelIdxList{n});
                 oD.WellOut.fData.(fStates{fI}).sumIntesity(n,aI) = sum(double(vals(:)));
               end
           end
       else
          for fI = 1:numel(fStates)
               for n = 1:numObj
                 vals = fluorIm.(fStates{fI}).stitchImage(CC.PixelIdxList{n});
                 oD.WellOut.fData.(fStates{fI}).sumIntesity(n,aI) = sum(double(vals(:)));
               end
           end
       end
    end

    if aI==1
        
        oD.WellOut.numObj(aI)     = numObj;
        oD.WellOut.prevBwImage    = bwImage;
        oD.WellOut.prevLblImage   = bwLabel;
        oD.WellOut.prevFFT        = fT1;
        oD.WellOut.prevbw         = sImage>bwThresh;
        oD.WellOut.initDateNum    = monData.StartODELAY;
        oD.WellOut.TimePoints(aI) = StitchMeta.ImageDateNum-oD.WellOut.initDateNum;
%         oD.WellOut.prevImageStats = ImageStats;
       fNames = fieldnames(StitchMeta);
        for n = 1:numel(fNames)
            oD.WellOut.StitchMeta(aI).(fNames{n}) = StitchMeta.(fNames{n});
        end

        oD.WellOut.ObjectNext = NaN(size(oD.WellOut.ObjectNext));
        oD.WellOut.ObjectTrack = NaN(size(oD.WellOut.ObjectNext));
        oD.WellOut.ObjectCentX = NaN(size(oD.WellOut.ObjectNext));
        oD.WellOut.ObjectCentY = NaN(size(oD.WellOut.ObjectNext));
        
        
        for n = 1:numObj
            oD.WellOut.ObjectTrack(n,aI) = n;
            oD.WellOut.ObjectArea(n,aI) = FilledArea(n);
            oD.WellOut.ObjectCentX(n,aI) = xyCent(n,1);
            oD.WellOut.ObjectCentY(n,aI) = xyCent(n,2);
        end

   else
        timeVec     = oD.WellOut.TimePoints;
        objectNext  = oD.WellOut.ObjectNext;
        objectTrack = oD.WellOut.ObjectTrack;
        objectArea  = oD.WellOut.ObjectArea;
        objectCentX = oD.WellOut.ObjectCentX;
        objectCentY = oD.WellOut.ObjectCentY;
        prevNumObj  = oD.WellOut.numObj(aI-1);
        xyDisp      = oD.WellOut.xyDisp;
    
    %% Centroid Association
         % Figure out what the image shift is from the previous Images
         imDim = size(sImage);
         bw1 = sImage>bwThresh;
         bw2 = oD.WellOut.prevbw;
         % Use FFT phase corelation to determin the offet
         fT2 = oD.WellOut.prevFFT;
         
         fT = fT1.*conj(fT2);
         fTabs = fT./abs(fT);
         fmag1 = ifft2(fTabs);
         fmag1(1) = 0;%the first index of fmag is always 1 so ignor it.
         [~, ind] = max(fmag1(:)); 
         [r,c] = ind2sub(imDim, ind);
         row = min([imDim(1)-r, r]);
         col = min([imDim(2)-c, c]);
         xyDisp(aI,:) = [row,col];
         
         % Once the offset is found the direction is not known and could be
         % sementric in 4 directions.  So test all for and look for maximum
         % overlap of the images.
         v1 = [0,0; 1,0; 0,1; 1,1];
         v2 = double(~v1);  % Positive direction in one image is the negative direction for the other.
         rowVec = 1:imDim(1);
         colVec = 1:imDim(2);
         rv1 = zeros(1,imDim(1));
         rv2 = zeros(1,imDim(1));
         cv1 = zeros(1,imDim(2));
         cv2 = zeros(1,imDim(2));
         cond = zeros(1,4);
   
        % Calculate difference between centroids in two sequential images
        % rows are the index of centroids in the first image
        % columns are the index of centroids in the second image or the next image
        % in this way the ObjectNext is the row index of the next image timepoint

         c = 1;
         for n = 1:4
             sw1 = zeros(imDim(1) + row , imDim(2) + col);
             rv1 = rowVec + v1(n,1)*row;
             cv1 = colVec + v1(n,2)*col;
             rv2 = rowVec + v2(n,1)*row;
             cv2 = colVec + v2(n,2)*col;
             sw1(rv1, cv1) = bw1;
             sw1(rv2, cv2) = sw1(rv2, cv2).*bw2;
             cond(c) = sum(sw1(:));
             c = c+1;
         end
         
         [~, ind] = max(cond);        
%       this gives the overlap vector for aligning the images
%       Set image Dimensions so they are identical.
        imageDimA = size(bwImage);
        imageDimB = size(oD.WellOut.prevBwImage);
        maxdim = max([imageDimA; imageDimB],[],1);
        
        %To do include translation from images earlier.
        alDim = floor((maxdim-imageDimA)./2)+1;
        auDim = maxdim-ceil((maxdim-imageDimA)./2);
        
        blDim = floor((maxdim-imageDimB)./2)+1;
        buDim = maxdim-ceil((maxdim-imageDimB)./2);
        
        arV = (alDim(1):auDim(1)) + v1(ind,1)*row;
        acV = (alDim(2):auDim(2)) + v1(ind,2)*col;
        brV = (blDim(1):buDim(1)) + v2(ind,1)*row;
        bcV = (blDim(2):buDim(2)) + v2(ind,2)*col;
               
        A = zeros(maxdim + [row, col]);
        B = zeros(maxdim + [row, col]);
        aLbl = zeros(maxdim + [row, col]);
        bLbl = zeros(maxdim + [row, col]);
        
        A(arV,acV)  =  bwImage;
        B(brV,bcV)  =  oD.WellOut.prevBwImage;
        aLbl(arV,acV)  =  bwLabel;
        bLbl(brV,bcV)  =  oD.WellOut.prevLblImage;
        
% Multiply black and white Images together.  This makes a mask
% where colonies overlap.

        M = A.*B;
        ALbl = aLbl.*M;%Current Labeled Image
        BLbl = bLbl.*M;%Prev Labeled Image
        CCM = bwconncomp(M, 8);
        for m = 1:CCM.NumObjects
            
            aVals = ALbl(CCM.PixelIdxList{m});
            bVals = BLbl(CCM.PixelIdxList{m});
            aDiff = diff(aVals);
            bDiff = diff(bVals);
            if ~isempty(aVals) && ~isempty(aVals) &&... % Check to make sure that aVals and bVals are not empty.
               sum(aDiff(:)) == 0 && sum(bDiff) == 0 &&... % make sure that all values are identical or the diff is zero
               bVals(1) ~=0 && aVals(1)~=0 % also make sure the values don't point to nothing.
           
               objectNext(bVals(1),aI-1) = aVals(1);
            end
        end
       
% Place objects that were linked in the Object Next list into an easier to
% address Object Track List.  
% Check for douplacates in the objectNext column array and remove them
    
    objectNext(objectNext==0) = NaN;
    rNum = ~isnan(objectNext(:,aI-1));
    if sum(rNum)>0
        objCnt = accumarray(objectNext(rNum,aI-1),1);
        idVec = 1:numel(objCnt);
        id = idVec(objCnt>1);
        for n = id
            inds = objectNext(:,aI-1) == n;
            objectNext(inds,aI-1) = NaN;
        end
    end
    
    % Place objects that were linked in the Object Next list into an easier to
    % address Object Track List.  
    if aI == 2
        objectTrack(:,1)= 1:size(objectTrack,1);
    end

    nextInds = 1:size(objectNext,1);
    nextList = nextInds(~isnan(objectNext(:,aI-1)));
    sizeNext = size(objectNext,1);
    sizeTrack = size(objectTrack,1);
    if sizeNext > sizeTrack
        objectTrack = [objectTrack;NaN(sizeNext-sizeTrack,size(objectTrack,2))];
    end
    
    if ~isempty(nextList)
        for n = 1:nextList(end)
            objTemp = objectTrack(n,aI-1);

            if n<size(objectTrack,1) &&...
               ~isnan(objTemp) &&...
               objTemp~=0 &&...
               ~isnan(objectNext(objTemp,aI-1)) &&...
               objectNext(objTemp,aI-1)<=sizeTrack
               objectTrack(n,aI) = objectNext(objTemp,aI-1);
            end

        end
    end
    
    if aI ==2
        
       for m = 1:size(oD.WellOut.RegionStats(1).FilledArea,1)
         if ~isnan(objectTrack(m,1))
             ind = objectTrack(m,1); 
             if sum(objectTrack(:,1)==ind)==1

              objectArea(m,1)  = oD.WellOut.RegionStats(1).FilledArea(ind,1);
              objectCentX(m,1) = oD.WellOut.RegionStats(1).Centroid(ind,1);
              objectCentY(m,1) = oD.WellOut.RegionStats(1).Centroid(ind,2);
              if isfield(oD.WellOut, 'fData')
                for fI = 1:numel(fStates) 
                  oD.WellOut.fData.(fStates{fI}).trackIntesity(m,1) = ...
                        oD.WellOut.fData.(fStates{fI}).sumIntesity(ind,1);
                end
              end

              end
         end
        end
    end  
    
    if size(objectArea,1)<size(objectTrack,1);
          objectArea  = [objectArea; NaN(size(objectTrack,1)-size(objectArea,1),size(objectTrack,2))];
          objectCentX = [objectCentX;NaN(size(objectTrack,1)-size(objectArea,1),size(objectTrack,2))];
          objectCentY = [objectCentY;NaN(size(objectTrack,1)-size(objectArea,1),size(objectTrack,2))];
    end
          
    %Remove Tracks that are fewer than 5 timepoints long.  
    for m = 1:size(objectTrack,1)
         if ~isnan(objectTrack(m,aI)) && objectTrack(m,aI)~=0
             ind = objectTrack(m,aI); 
            
             if sum(objectTrack(:,aI)==ind)==1
                  objectArea(m,aI)  = oD.WellOut.RegionStats(aI).FilledArea(ind,1);
                  objectCentX(m,aI) = oD.WellOut.RegionStats(aI).Centroid(ind,1);
                  objectCentY(m,aI) = oD.WellOut.RegionStats(aI).Centroid(ind,2);
                  if isfield(oD.WellOut, 'fData')
                      for fI = 1:numel(fStates) 
                       oD.WellOut.fData.(fStates{fI}).trackIntesity(m,aI) = ...
                            oD.WellOut.fData.(fStates{fI}).sumIntesity(ind,aI);
                      end
                  end
             else
                 objectTrack(objectTrack(:,aI) == ind, aI) = NaN;
             end
         else
             objectTrack(m,aI) = NaN;
         end
    end


    %Generate Timepoints for this Data-Set
%     if oD.WellOut.initDateNum > StitchMeta.ImageDateNum 
%         oD.WellOut.initDateNum = StitchMeta.ImageDateNum;
%     end
    timeVec(aI) = StitchMeta.ImageDateNum-oD.WellOut.initDateNum;

    %Time is stored as a value using the MATLAB function "datenum".  The time 
    %difference is calculated and then converted to minutes.  
    DateTime = datevec(timeVec);

    % Calculate the Time points in minutes
    ClockPoints = DateTime(:,3).*24*60 ...
                 +DateTime(:,4).*60 ...
                 +DateTime(:,5)...
                 +DateTime(:,6)./60;

    oD.WellOut.ClockPoints = ClockPoints;
    oD.WellOut.timeVec     = timeVec;
    oD.WellOut.TimePoints  = timeVec;
    oD.WellOut.ObjectTrack = objectTrack;
    oD.WellOut.ObjectArea  = objectArea;
    oD.WellOut.ObjectCentX = objectCentX;
    oD.WellOut.ObjectCentY = objectCentY;
    oD.WellOut.ObjectNext  = objectNext;
    oD.WellOut.StitchMeta(aI)  = StitchMeta;

    oD.WellOut.prevBwImage  = bwImage;
    oD.WellOut.prevLblImage = bwLabel;
    oD.WellOut.prevFFT      = fT1;
    oD.WellOut.prevbw       = bw1;
    oD.WellOut.numObj(aI)   = numObj;
    oD.WellOut.xyDisp       = xyDisp;
    end
%     fprintf(['Well: ', oD.WellData.WellID, ' Image: ',num2str(aI),'\n']);
    analyzedIm(aI) = true;
end
%  assignin('base', 'xyDisp',xyDisp);
%  assignin('base','oD', oD);

%% FitData
%Check to see if the object area shrinks between the first and second
%timepoint if it does the focus may be bad and use the second time point to
%start fitting.
if isempty(analyzIndex)
   analyzedIm(1:numel(currentFiles)) = 1; 
else
    if size(objectArea,2)>1; %Catch the problem if ObjectArea is filled
        dArea = diff(log2(objectArea(:,1:2)),1,2);
        if dArea<-1; tp1 = 2;else tp1 = 1;end
    else tp1 = 1;
    end
objectArea(objectArea==0) = NaN;
diffAreaStd = abs(diff(log2(objectArea),1,2));
dbInds = diffAreaStd>0.6;
bgSteps = cumsum(dbInds,2)~=0;
objectArea(bgSteps)= NaN;

trackLength = sum(~isnan(objectArea),2);
trackVec = 1:size(objectArea,1);
trackInds = trackVec(trackLength>10);

if size(objectArea,1)>size(oD.Tracks2.ObjectInfo.FitDataGompDT,1)
    oD.Tracks2.ObjectInfo.FitDataGompDT = NaN(size(objectArea,1),17);
    oD.Tracks2.ObjectInfo.FitDataLinear = NaN(size(objectArea,1),5); 
end
    
%Set up fitting with gompertsfit functions.  
if ~isempty(trackInds)
    for n = trackInds
            inds = ~isnan(objectArea(n,:));
            tData    = ClockPoints(inds)';
            areaData = objectArea(n,inds);
            initTvTLag = oD.Tracks2.ObjectInfo.FitDataGompDT(n,1:4);
            oD.Tracks2.ObjectInfo.FitDataLinear(n,1:end) = linearFitParam(tData, areaData);
            oD.Tracks2.ObjectInfo.FitDataGompDT(n,1:end-1) = gompertzFitDT(tData, areaData, initTvTLag);
        
    end
end

%% FitData Column Key 
%Col     1   2   3   4     5     6     7     8     9       10        11        12       13
%Header 'a' 'b' 'c' 'd' 'fval' 'Tlag' 'Td' 'Tex' 'ATex' 'Aplateau' 'TdFlag' 'TexFlag' 'TVmax'
%     a:        Gompertz function value, a for y=a+b*exp(-exp(-ct+d))
%     b:        Gompertz function value, b for y=a+b*exp(-exp(-ct+d))
%     c:        Vmax Td=(2*log(2))/(3*Vmax); or Tlag
%     d:        Tlag or DT
%     fval:     Sum of squares (modeled versus experimental) from fmincom
%     Tlag:     Extracted lag time
%     Td:       Extracted doubling time
%     Tex:      Time at which modeled curve exits exponential phase
%               Defined by second x-intercept of the third derivative  
%     ATex:     Area of the modeled microcolony at exit from exponential phase
%     Aplateau: Area of the modeled microcolony once carrying capacity is reached  
%     TdFlag:   1 if maximum growth, the point at which doubling time is
%               calculated, occurs after the last experimental data point.
%               0 if doubling time is calculated at a time point within the
%               period monitored.
%     TexFlag:  1 if exit from exponential phase occurs after the last
%               experimental data point, 0 if within period monitored. 
%     TVmax:    Time at which maximal growth rate is achieved and doubling
%               time is calculated.
%     For the above columns 'NaN' is entered for groups excluded from
%     modeling due to paucity of experimental data or lack of growth. 
% 
oD.Tracks2.ObjectInfo.TimePoints  = ClockPoints;
oD.Tracks2.ObjectInfo.ObjectArea  = objectArea;
oD.Tracks2.ObjectInfo.ObjectTrack = objectTrack;
if isfield(oD.WellOut, 'fData')
    oD.Tracks2.ObjectInfo.fData = oD.WellOut.fData;
end
% %% Package Data for Saving

%% If Data is done processing SaveDataState
    SaveDataState(oD, monData, well)
end

end %function

function [k_segment] = Image_Segment(AnImage, counts)
%% Image Segmentation
%==========================================================================
%  Image Segment breaks an image down into descreat parts and uses the
%  maximum pixel value in the individual reagion. Then from each region the
%  max value is accumulated and the pixel value with the highest number of
%  votes becomes the threshold value. 
%==========================================================================
      % Determin Image Bit Depth if not already known
      classname = class(AnImage);
      bitdepth = intmax('uint16');
   
      % Segment image by creating vectors from linear segments of the image
      [nrow, ncol] = size(AnImage);
      
      dcol = 1:25:ncol;
      drow = 1:25:nrow;
      
      max_pix_hist = zeros(1,bitdepth);
      
      for c = 1:numel(dcol)-1
          for r = 1:numel(drow)-1
              imseg = AnImage(drow(r):drow(r+1),dcol(c):dcol(c+1));
              max_pix = floor(max(imseg(:))+1);
              if max_pix>bitdepth; max_pix = bitdepth; end
              max_pix_hist(1,max_pix) = max_pix_hist(1,max_pix)+1;
          end
      end  
          
      %Find max value of histogram segments
      [~, k_segment] = max(max_pix_hist(2:bitdepth-1));
      k_segment = k_segment+1;
    
end

function [RawImage, RawMeta, RawStack, fluorIm,fTrans] = ODELAM_MatReadv5(fileName, magnification, background)
%% ODELAY_MetaMorph_Read
%  This function reads in sequential images written and 
%  uses phase correlation to stitch them into a single image.  
% Author: Thurston Herricks
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% email:  Thurston.Herricks@systemsbiology.org


numTry = 0;
success = false;

while numTry<=10 & success == false
    try
        imageData = load(fileName);
        success = true;
    catch
        pause(3);
        numTry = numTry+1;
    end
end

CorrSwitch=1;
PixSize = 6.45/magnification;


numImage = size(imageData.rawImage,3);
xDim = size(imageData.rawImage,2);
yDim = size(imageData.rawImage,1);
% Get stage position values for each image in the group
stageXYZ = imageData.xyzTime(:,1:3);

RawStack = zeros(max(yDim), max(xDim), numImage, 'uint16');

% Load Background Image if it exists
if isempty(background)
    background = zeros(yDim,xDim,'uint16');
end
%Read in images into RawStack

missingImage = false;

for imNum = 1:numImage 
       RawStack(:,:,imNum) = imageData.rawImage(:,:,imNum)-background;
end

%Correct Stage positions and set them to a minimum coordiante axis.  
minstageX = min(stageXYZ(1,1));
maxstageX = max(stageXYZ(:,1));
minstageY = min(stageXYZ(1,2));
maxstageY = max(stageXYZ(:,2));
stageXYrel = abs([stageXYZ(:,2)-minstageY,...
                  stageXYZ(:,1)-minstageX]);

pixXYrel = stageXYrel./PixSize;
overlapXY = NaN(numImage);

% Determin relative displacement between images
distRows = numImage*(numImage-1)/2; %number of comparisons 
distXY = zeros(distRows,10);
tempInds = false(distRows,1);
% 
cnt = 1;
for col = 1:numImage-1
    for row = col+1:numImage
        vecXY = pixXYrel(row,:)-pixXYrel(col,:);
        distXY(cnt,3) = sqrt(sum(vecXY.^2));
        distXY(cnt,1:2)= [col,row];
        distXY(cnt,4:5) = (vecXY);
        distXY(cnt,6) = atan2(vecXY(1),vecXY(2));
        
        if vecXY(1)==0 & vecXY(2)<xDim & (prod(vecXY) == 0)
            tempInds(cnt) = true;
            overlapXY(row,col) =  1;
            overlapXY(col,row) =  3;         
        elseif vecXY(2)==0 & vecXY(1)<yDim & (prod(vecXY) == 0)
            tempInds(cnt) = true;
            overlapXY(row,col) =  2;
            overlapXY(col,row) =  4;  
        end
        cnt = cnt+1;
    end
end


imOverlap = distXY(tempInds,:);
minDist = imOverlap;
 
%  TODO figure out how to use overlapXY and image order to determine image
%  comparison order imCompOrd
   
% constrain by calculate upper right hand corner for each image and see which is smallest    
quadDisp = [0,xDim-1;
            yDim-1,xDim-1;
            yDim-1,0;
            0,0];

ind = zeros(size(minDist,1),1);
magt = zeros(size(minDist,1),1);
fTN = zeros(yDim,xDim,numImage,'double');
for n = 1:numImage
    fTN(:,:,n) = fft2(double(RawStack(:,:,n)));
end
fTrans = fTN(:,:,5);
for n = 1:size(minDist,1)

     %get the FFT of the two images
     fT1 = fTN(:,:,minDist(n,1));
     fT2 = fTN(:,:,minDist(n,2));

     %Calculate the cross-correlation of the images
     fT = fT1.*conj(fT2);
     fTabs = fT./abs(fT);
     fmag1 = ifft2(fTabs);
     fmag1(1) = 0;%the first index of fmag is always 1 so ignor it.

     [magt(n), ind(n)] = max(fmag1(:)); 
end
[row,col] = ind2sub([yDim,xDim], ind);
minDist(:,7:8) = [row,col];
dXYPhase = [row,col];

% calculate the displacement vector diffPhaseXY which is the XY
% displacement from the stage corrdinates.  The smallest displacement
% in this case is the correction since larger displacements are
% probably incorrect.

magDT = zeros(size(minDist,1),4);
diffPhaseXY = zeros(size(minDist,1),2);
angDT = zeros(size(minDist,1),4);
magMin = zeros(size(minDist,1),1);
angMin = zeros(size(minDist,1),1);
angCheck = zeros(size(minDist,1),1);
magCheck = zeros(size(minDist,1),1);

for n = 1:size(minDist,1)
        D = minDist(n,4:5);
    for j = 1:4

        T = minDist(n,7:8)-quadDisp(j,:);
        magDT(n,j) = sqrt(sum((D-T).^2));
        angDT(n,j) = acos(dot(T,D)/(sqrt(sum(D.^2))*sqrt(sum(T.^2))));

    end
    [magCheck(n), magMin(n)]  = min(magDT(n,:));
    [angCheck(n), angMin(n)]  = min(angDT(n,:));
    T = minDist(n,7:8) - quadDisp(magMin(n),:);
    TD(n,:) = T-D;
end
%round the angles between the vectors so that the numbers that are
%close can be calculated.
minDist(:,9:10) = TD;


% Find the Regions with the same displacement vectors
sameVec = false(size(minDist,1),size(minDist,1));
 for m = 1:size(minDist,1)
         sameVec(m,:) = minDist(:,6) == minDist(m,6);
 end
 redVec = sum(rref(double(sameVec)),2)>0;
 baseVec = sameVec(redVec,:);
 
%round the angles between the vectors so that the numbers that are
%close can be calculated.
angFlag = zeros(size(angCheck));
magFlag = zeros(size(magCheck));
angProbs = zeros(size(angCheck));
magProbs = zeros(size(magCheck));
for n = 1:numel(angCheck); angFlag(n) = sum(abs(angCheck(sameVec(n,:))'-angCheck(n))>0.01); end
for n = 1:numel(magCheck); magFlag(n) = sum(abs(magCheck(sameVec(n,:))'-magCheck(n))>4); end


%This means there is a bad vector as all should be identical
magProbs(baseVec(1,:)) = magFlag(baseVec(1,:))~=min(magFlag(baseVec(1,:)));
magProbs(baseVec(2,:)) = magFlag(baseVec(2,:))~=min(magFlag(baseVec(2,:)));
angProbs(baseVec(1,:)) = angFlag(baseVec(1,:))~=min(angFlag(baseVec(1,:)));
angProbs(baseVec(2,:)) = angFlag(baseVec(2,:))~=min(angFlag(baseVec(2,:)));
numProbs = sum(magProbs|angProbs);

if numProbs~=0
    
    electAng = mean(angCheck(~angProbs)); 
    electMag = mean(magCheck(~magProbs)); 
    vecList = 1:size(minDist,1);
    fixInds = magProbs|angProbs;
    fixList = vecList(fixInds);
       
    for n = fixList
         tempV = round(minDist(n,4:5));
         
         sameInd = sameVec(n,:);
         sameInd(fixList)=0;
%          for m = vecList(~fixInds)
%             if tempV(1) == round(minDist(m,4)) && tempV(2) == round(minDist(m,5))
%                 sameInd(m) = true;
%                 % Need check for quadrant correction
%             end 
%          end
         TD(n,:) = mean(TD(sameInd,:),1);
    end
    
   if sum(isnan(TD(:)))>0
       TD = zeros(size(minDist(:,4:5)));
   end
%         
end
%    
% Find vectors paths to each image;
imagePath = NaN(numImage);
numSteps = zeros(numImage,1);

imagePath(1,1) = 0;
im = 2;
prevImage = 0;
cnt = 0;

for imNum =  2:numImage
    revOrder = zeros(1,numImage);
    prevImage = imNum;
    img = imNum;
    cnt = 0;
    while prevImage~=1 & cnt < numImage;
        cnt = cnt+1;
        row = find(minDist(:,2)==prevImage,1,'first');
        prevImage = minDist(row,1);
        revOrder(cnt) = row; % this is the reverse order of the path and will be flipped
    end
    imagePath(imNum,1:cnt) = fliplr(revOrder(1:cnt));
    numSteps(imNum) = cnt;
end

% correct the phaseXY from displacements in the individual images to
% the cumulative displacement from correcting each image
phaseXY = zeros(numImage,2);    
cumPhaseXY = zeros(numImage,2);
for imNum = 2:numImage
        cumPhaseXY(imNum,:) = sum(TD(imagePath(imNum,1:numSteps(imNum)),:),1);    
end
phaseXY = round(pixXYrel+cumPhaseXY);
% Finnally zero out the corrdinate system for assembling the images.
phaseXYcor  = [phaseXY(:,1)-min(phaseXY(:,1)),...
               phaseXY(:,2)-min(phaseXY(:,2))];

% ToDo:  Check displacements and make sure they average out across all directions to other images.
            
dXY = diffPhaseXY;          
      
% Convert to pixels
% phaseXYcor=[106,    0;
%              53, 1183;
%               0, 2366;
%             987,   42;
%             934, 1224;
%             881, 2408;
%            1872,   84;
%            1819, 1266;
%            1766, 2450];
    
    
imPix = phaseXYcor;%  imPix = stageXYrotcorr;
% Determin size of stitched image
stitchDim = max(imPix)+[yDim, xDim];
    
% Create array for stitched image
stitchImage   = zeros(stitchDim, 'uint16');
stitchDevisor = zeros(stitchDim,'uint16');
   
% Create an array of indicies that represent pixels in the original
yindp = 1:yDim;
xindp = 1:xDim;

corindy = zeros(numImage,yDim);
corindx = zeros(numImage,xDim);

% Generate a stitchDevisor...this is a matrix that gives the number of
% times an individual pixel is overlapped so that each pixes is averaged
% together appropriately.
for n = 1:numImage
        corindy(n,:) = imPix(n,1)+yindp;
        corindx(n,:) = imPix(n,2)+xindp;
        stitchDevisor(corindy(n,:), corindx(n,:)) = stitchDevisor(corindy(n,:), corindx(n,:))+1;                 
end

imagevec = 1:numImage;  % used to trouble shoot which images are added
for j = 1:numel(imagevec)
    n = imagevec(j);
    imagedata = RawStack(:,:,n)./stitchDevisor(corindy(n,:), corindx(n,:));
    stitchImage(corindy(n,:), corindx(n,:))= imagedata+stitchImage(corindy(n,:), corindx(n,:));    
    
end 

if isfield(imageData, 'fluorImage')
    
    fluorExist = true;
    fluorList = fieldnames(imageData.fluorImage);
    numfluorIm = numel(fluorList);   
    for n = 1:numfluorIm
       fluorIm.(fluorList{n}).stitchImage = zeros(stitchDim, 'uint16');  
    end
    
    for n = 1:numfluorIm  
       for j = 1:numel(imagevec)
            m = imagevec(j);
            imagedata = imageData.fluorImage.(fluorList{n}).rawImage(:,:,m)./stitchDevisor(corindy(m,:), corindx(m,:));
            fluorIm.(fluorList{n}).stitchImage(corindy(m,:), corindx(m,:))= imagedata+fluorIm.(fluorList{n}).stitchImage(corindy(m,:), corindx(m,:));    
       end
    end
else
    fluorIm = [];
end

OverlapCount = 0;
if missingImage == true
    OverlapCount = -1;
end


RawImage = stitchImage;
RawMeta(1).StagePos = stageXYZ;
% RawMeta(1).ImagePhaseDxy = dXY;
RawMeta(1).ImagePhasePos = phaseXYcor;
RawMeta(1).zFocusVals = imageData.zRecord;
RawMeta(1).ImageDateNum = imageData.xyzTime(5,4);
RawMeta(1).ImageDim = size(stitchImage);
RawMeta(1).dxyMax = [0,0];
RawMeta(1).OverlapCount = OverlapCount;

end

function prmsall = linearFitParam(tdata, iareas)

        xVals = [tdata',ones(numel(tdata),1)];

        [x, stdx, mse] = lscov(xVals, log2(iareas'));
        prmsall(1,1:2) = x';
        prmsall(1,3:4) = stdx';
        prmsall(1,5) = mse;
        
end

function prmsAll = gompertzFitDT(tData, iareas, initVals)
%% gompertzFitParam calculates an initial estimate of the Parameterized 
%% Gompertz Function 
% using a coarse grid optimization and then attempts to find 
% a constrained minimum of the Gompertz function at this initial estimate 
% using fmincon MATLAB function.
% Notes:
%==========================================================================
%% Author: Alexander V. Ratushny
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% Email: alexander.ratushny@systemsbiology.org
% or
% Seattle Biomed
% 307 Westlake Ave N
% Seattle, WA 98109 USA
% Email: alexander.ratushny@seattlebiomed.org
%==========================================================================
%% Author: Thurston Herricks 
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% Emails: 
% Thurston.Herricks@systemsbiology.org
%==========================================================================
% October 2010; Last revision: 2016/03/01
%
    numtimepoints = sum(~isnan(iareas));

    idata = log2(iareas(1:numtimepoints));
    tdata =  tData(1:numtimepoints);
    Nsteps=40;
    areaMax=max(iareas);  
    factor=1.05;
    smthArea = smoothdata(idata,'movmean',5);
    x = tdata([1,numtimepoints]);
    y = smthArea([1,numtimepoints]);
    m = diff(y)/diff(x);
    yVals = m*tdata + y(1)-m*x(1);
    diffVals = smthArea-yVals;
    cumVals  = cumsum(diffVals);
    

    [~, lagInd] = min(diffVals);
    [~, texInd] = max(diffVals);
    [~ ,vmxInd] = min(cumVals);
    
%     numPos = sum(cumVals(vmxInd:end)>0);
   
   if  lagInd < vmxInd && vmxInd < texInd
        estTlag  = tdata(lagInd);
        estDT = tdata(vmxInd)-tdata(lagInd);
        estB  = idata(vmxInd); 
   elseif lagInd<vmxInd && (texInd<lagInd | texInd<vmxInd)
       estTlag  = tdata(lagInd);
       estB  = idata(vmxInd); 
       estDT = tdata(vmxInd)-tdata(lagInd);

   elseif lagInd<texInd && (vmxInd<lagInd | vmxInd<texInd)
       estTlag  = tdata(lagInd);
       estB     = idata(texInd); 
       estDT    = (tdata(texInd)-tdata(lagInd))./2;
       
   else
    % Use course grid optimization function findPrmsGompF to find 
    % a local minima based on the 
    
        vecDT=linspace(1,2*tdata(end),Nsteps);
        bmin = 0;
        bmax = 10;
        vecTlag = linspace(1,tdata(end),Nsteps);
        vecB = linspace(bmin,bmax,Nsteps);
        [estB, estTlag,estDT,  fssq1] = findPrmsGompBDt(vecB,  vecTlag, vecDT, tdata, idata);
    
    end
    meanArea = mean(idata(1:5));
    stdArea  = std(idata(1:5),0,2);
    Klag = log((3+sqrt(5))/2);
    prms0 = [meanArea estB  estTlag estDT];

    meanArea = mean(idata(1:5));
    stdArea  = std(idata(1:5),0,2);
    Klag = log((3+sqrt(5))/2);
    
    aLow = meanArea-3*stdArea;
    aUp  = meanArea+3*stdArea;
    dTLow = 1;
    dTUp  = max(tdata);
    bLow  = 0;
    bUp   = 10;
    lagLow = 0;
    lagUp = max(tdata);
    
    lb    = [aLow   bLow   lagLow dTLow];
    ub    = [aUp    bUp    lagUp  dTUp];
    gompFunc = @(prms)gompParBDt(prms, tdata, idata);
    options = optimoptions('fmincon','Algorithm','interior-point',...
                           'Display','off');
    Aeq = [0, 1, 0, 1]; B = [tdata(end)];
    warning('off','all');
    [prms, fssq, exitflag] = fmincon(gompFunc,prms0,Aeq,B,[],[],lb,ub,[],options);
    warning('on','all');
    a    = prms(1);
    b    = prms(2);
    Tlag = prms(3);
    dT   = prms(4);
    
    Klag = log((3+sqrt(5))/2);
    Kex  = log((3-sqrt(5))/2);
    c = Klag/dT;
    d = Tlag*c+Klag;
    
    TVmax = d/c;
    Tex = (d-Kex)/c;
    Vmax = b*c*exp(-1);
    Td=1/Vmax;

    if(TVmax>tdata(end))
        TdFlag=1;
    else
        TdFlag=0;
    end
    
    if(Tex>tdata(end))
        TexFlag = 1;
    else
        TexFlag = 0;
    end

    Tplat = 0;
    
    ATex     = gompBDt(prms, Tex);
    Aplateau = gompBDt(prms,1e50);
    prmsAll=[prms fssq Tlag Td Tex ATex Aplateau TdFlag TexFlag TVmax Tplat exitflag fssq/numtimepoints]; 
                
end
 
function [estB, estTlag, dTmin, fmin ] = findPrmsGompBDt(vecB, vecTlag, vecDt, tdata, idata )
%% findPrmsGompF calculates an initial estimate of the Gompertz function using a coarse grid optimization.
%
% Author: Alexander V. Ratushny
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% Email: alexander.ratushny@systemsbiology.org
% or
% Seattle Biomed
% 307 Westlake Ave N
% Seattle, WA 98109 USA
% Email: alexander.ratushny@seattlebiomed.org
%% Author Thurston Herricks
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% Email: Thurston.Herricks@systemsbiology.org
% October 2010; Last revision: 2016/03/01

    flag=0;
    a = nanmean(idata(1:5));
    K = log((3+sqrt(5))/2);
    numTPs = numel(vecTlag);
    tVec = 1:numTPs;
    
   for B = vecB
       for tp = tVec(1:end-1)
           tlag = vecTlag(tp);
          
           vecDt = vecTlag(tp+1:end)-vecTlag(tp);
             
           for dT = vecDt
                yn=a+B.*exp(-exp((K/dT).*(dT+tlag-tdata)));
                ifunc = nansum((idata-yn).^2);
                if(flag==0||(flag==1&&ifunc<fmin))
                    fmin    = ifunc;
                    dTmin   = dT;
                    estB    = B;
                    estTlag = tlag;
                    flag=1;
                end
            end
        end
    end

end

function f = gompParBDt(x, tdata, idata)
%% gompFparam calculates least squares distance of the paramerized gompertz function.
%
% from the logarithm of Gompertz-approximated growth curve values data through time.
%
%% Author: Alexander V. Ratushny
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% Email: alexander.ratushny@systemsbiology.org
% or
% Seattle Biomed
% 307 Westlake Ave N
% Seattle, WA 98109 USA
% Email: alexander.ratushny@seattlebiomed.org
%% Author: Thurston Herricks
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% Email: thurston.herricks@systemsbiology.org
% October 2010; Last revision: 2016/03/01
%%
Klag = log((3+sqrt(5))/2);

a    = x(1);
b    = x(2);
tlag = x(3);
dT   = x(4);

yn=a + b.*exp(-exp((Klag/dT).*(dT+tlag-tdata)));
f=nansum((yn-idata).^2);
         
 end

function yn = gompBDt(x, tdata)
%% gompFparam calculates least squares distance of the paramerized gompertz function.
%
% from the logarithm of Gompertz-approximated growth curve values data through time.
%
%% Author: Alexander V. Ratushny
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% Email: alexander.ratushny@systemsbiology.org
% or
% Seattle Biomed
% 307 Westlake Ave N
% Seattle, WA 98109 USA
% Email: alexander.ratushny@seattlebiomed.org
%% Author: Thurston Herricks
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% Email: thurston.herricks@systemsbiology.org
% October 2010; Last revision: 2016/03/01
%%

Klag = log((3+sqrt(5))/2);

a    = x(1);
b    = x(2);
tlag = x(3);
dT   = x(4);

yn= a + b.*exp(-exp((Klag/dT).*(dT+tlag-tdata)));
         
end 

function oD = LoadDataState(well, monData)
%% Use ODELAYLoadData to load well data + previous image data from 
%Use uigetfile to load data from datasets without using the comand line

cd(monData.saveDir);
fileList = dir('*ODELAYData.mat');
ExData = load(fullfile(fileList(well).folder,fileList(well).name));
ExFields = fieldnames(ExData);
    
for ctr = 1:size(ExFields,1)
    switch ExFields{ctr}

        case 'WellDataTemp'

            WellDataFields = fieldnames(ExData.WellDataTemp);

            for fieldNum = 1:size(WellDataFields,1);
                oD.WellData.(WellDataFields{fieldNum}) = ...
                    ExData.WellDataTemp.(WellDataFields{fieldNum,:});

            end
           

        case 'WellOutTemp'
            WellOutFields = fieldnames(ExData.WellOutTemp);
             for fieldNum = 1:size(WellOutFields,1);
                    oD.WellOut.(WellOutFields{fieldNum}) = ...
                        ExData.WellOutTemp.(WellOutFields{fieldNum,:});
             end


        case 'Tracks2Temp'
            Tracks2Fields = fieldnames(ExData.Tracks2Temp);
             for fieldNum = 1:size(Tracks2Fields,1);
                    oD.Tracks2.(Tracks2Fields{fieldNum}) = ...
                        ExData.Tracks2Temp.(Tracks2Fields{fieldNum,:});
             end


    end
end
cd(monData.expDir);
end

function Background = ODELAY_MatCollectBackground(monData, mP)
%% ODELAY_Collect_Background
% function randomly sellects a number of images and averages them together
% to create a background image for background subtraction

%% Author: Thurston Herricks 
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% Emails: 
% Thurston.Herricks@systemsbiology.org
% Background = evalin('base', 'Background');
numWells = 5; %size(monData.wellDir,1);
imageTot = numWells*mP.numTiles;

% Load and average the first image from the first XX wells
for well = 1:numWells
    expDir   = monData.expDir;
    wellDir  = monData.wellDir{well};
    imageFile = [monData.wellDir{well},'_1.mat'];
    fileName = fullfile(expDir,wellDir, imageFile);
    
    numTry = 0;
    success = false;

    while numTry<=10 & success == false
        try
            imageData = load(fileName);
            success = true;
        catch
            pause(3);
            numTry = numTry+1;
        end
    end

    for cnt = 1:size(imageData.rawImage,3);
            RawImage = imageData.rawImage(:,:,cnt);
            if exist('accume','var')==0;
                accume = zeros(size(RawImage), 'double');
            end
            accume(:) = accume(:)+double(RawImage(:))./imageTot;    
    end

end

accume = accume-min(accume(:));

Background = uint16(floor(accume));

end

function SaveDataState(oD, monData, well)
 
%% Save well data for next round.
%Determin Variables in the Caller Workspace

stInd = regexp(monData.expDir, filesep);
[path,ExperimentName, ~] = fileparts(monData.expDir);

WellDataTemp = oD.WellData;
WellOutTemp  = oD.WellOut;
Tracks2Temp  = oD.Tracks2;

cd(monData.saveDir);

saveWellName =  [ExperimentName,'_Well_',num2str(well,'%0.2d'),'_ODELAYData'];
save(saveWellName, 'WellDataTemp','WellOutTemp','Tracks2Temp');

cd(monData.expDir);

end

function ODELAYSaveData

WellDataExist = false;
WellOutExist  = false;
Tracks2Exist  = false;
Tracks3Exist  = false;
%Determin Variables in the Caller Workspace
callerVariables  = evalin('caller','who');

if ~isdir('ODELAY Well Data')
    mkdir('ODELAY Well Data');
    cd('ODELAY Well Data');
    saveDir = cd;
else isdir('ODELAY Well Data')
    cd('ODELAY Well Data');
    saveDir = cd;
end

for n = 1:size(callerVariables,1)
    switch callerVariables{n,:}
         case 'monData'
            monData         = evalin('caller', 'monData');
            monData.saveDir = saveDir;
        case 'ExperimentName'
            ExperimentName       = evalin('caller', 'ExperimentName');
        case 'Experiment_Name'
            ExperimentName       = evalin('caller', 'Experiment_Name');
        case 'CentroidAssociationVersion'
            CentroidAssociationVersion   = evalin('caller','CentroidAssociationVersion');
        case 'WellProcessVersion'
            WellProcessVersion   = evalin('caller','WellProcessVersion');
        case 'DataDirectory'
            DataDirectory    = evalin('caller','DataDirectory');
            DataDirctoryExist = true;
        case 'ImageDirectory'
            ImageDirectory   = evalin('caller','ImageDirectory');
        case 'CurrentCD'
            CurrentCD        = evalin('caller', 'CurrentCD');
        case 'ImageVars'
            ImageVars        = evalin('caller','ImageVars');
        case 'TimeDilute'
            TimeDilute       = evalin('caller','TimeDilute');
        case 'numwells'
            numwells         = evalin('caller','numwells');
        case 'saveWellFile'
            saveDir     = evalin('caller','saveWellFile');
        case 'WellData'
            WellData         = evalin('caller','WellData');
            WellDataExist = true;
        case 'WellOut'
            WellOut          = evalin('caller','WellOut');
            WellOutExist = true;
        case 'Tracks2'
            Tracks2          = evalin('caller','Tracks2');
            Tracks2Exist = true;  
        case 'Tracks3'
            Tracks3          = evalin('caller','Tracks3');
            Tracks3Exist = true;
         case 'tiffImageFiles'
            tiffImageFiles   = evalin('caller','tiffImageFiles');
            tiffImageFilesExist= true;
    end
end


 
if WellDataExist&WellOutExist&Tracks2Exist&Tracks3Exist 
    for well = 1:numwells
       WellDataTemp = WellData(well);
       WellOutTemp  = WellOut(well);
       Tracks2Temp  = Tracks2(well);
       Tracks3Temp  = Tracks3(well);
       saveWellName =  [ExperimentName,'_Well_',num2str(well,'%0.2d'),'_ODELAYData'];
       save(saveWellName, 'WellDataTemp','WellOutTemp','Tracks2Temp', 'Tracks3Temp');

    end
elseif WellDataExist&WellOutExist&Tracks2Exist&~Tracks3Exist 
    for well = 1:numwells
       WellDataTemp = WellData(well);
       WellOutTemp  = WellOut(well);
       Tracks2Temp  = Tracks2(well);
       saveWellName =  [ExperimentName,'_Well_',num2str(well,'%0.2d'),'_ODELAYData'];
       save(saveWellName, 'WellDataTemp','WellOutTemp','Tracks2Temp');

    end
elseif WellDataExist&WellOutExist&~Tracks2Exist&~Tracks3Exist
     for well = 1:numwells
       WellDataTemp = WellData(well);
       WellOutTemp  = WellOut(well);
       saveWellName =  [ExperimentName,'_Well_',num2str(well,'%0.2d'),'_ODELAYData'];
       save(saveWellName, 'WellDataTemp','WellOutTemp');

     end
elseif WellDataExist&~WellOutExist&~Tracks2Exist&~Tracks3Exist
     for well = 1:numwells
       WellDataTemp = WellData(well);
       saveWellName =  [ExperimentName,'_Well_',num2str(well,'%0.2d'),'_ODELAYData'];
       save(saveWellName, 'WellDataTemp');

    end
end

cd(fileparts(cd));

saveVariables = {'CurrentCD';...
                 'ExperimentName';...
                 'DataDirectory';...
                 'ImageDirectory';...
                 'ImageVars';...
                 'monData';...
                 'numwells';...
                 'saveWellFile';...
                 'tiffImageFiles'};

save([ExperimentName,'_Index_ODELAYData'],saveVariables{1},'-v7.3'); 
for cnt = 2:numel(saveVariables)
     if exist(saveVariables{cnt},'var')
        save([ExperimentName,'_Index_ODELAYData'],saveVariables{cnt},'-append');
     end
end

end

function [Experiment_Name, ImageDirectory, DataDirectory, ImageVars] = ODELAY_importExpDescription(fileName)

%% Import the data
[~, ~, data] = xlsread(fileName, 'Sheet1', '');
data(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),data)) = {''};

for cntr = 1:size(data,1)

    switch data{cntr,1}
        case 'Experiment Name'
            Experiment_Name = data{cntr,2};
        case 'ImageDirectory'
            ImageDirectory = data{cntr,2};
        case 'DataDirectory'
            DataDirectory = data{cntr,2};
        case 'Dilution Time'
             ImageVars.TimeDilute =...
             datenum([data{cntr,2},' ',data{cntr,3}],'mm/dd/yyyy HH:MMPM'); 
        case 'Microscope ID'
            ImageVars.Microscope_ID = data{cntr,2};
        case 'CorrSwitch'
            ImageVars.CorrSwitch = data{cntr,2};
        case  'TiledRows'
            ImageVars.TiledRows = data{cntr,2};
        case  'TiledCols'
            ImageVars.TiledCols= data{cntr,2};
        case  'binningRate'
            ImageVars.binningRate= data{cntr,2};
        case  'Threshold_Method'
            ImageVars.Threshold_Method = data{cntr,2};
        case  'numwells'
            ImageVars.NumWells = data{cntr,2};
        case 'numImages'
            ImageVars.NumImages = data{cntr,2};
        case  'Magnification'
            ImageVars.Magnification = data{cntr,2};
        case  'PixSize'
            ImageVars.PixSize = data{cntr,2};
        case  'CameraRotaionAngle'
            ImageVars.CameraRotaionAngle= data{cntr,2};
        case 'collect_background'
            ImageVars.collect_background = data{cntr,2};
        case  'OffsetMult';
            ImageVars.OffsetMult = data{cntr,2};
        case  'image_per_well';
            ImageVars.image_per_well = data{cntr,2};
        case  'RotAngle';
            ImageVars.RotAngle = data{cntr,2};
        case  'phaseXYcor'
            ImageVars.phaseXYcor = [data{cntr:cntr+8,2};data{cntr:cntr+8,3}]';
        case 'StrainID Headers'
            ImageVars.StrainIDHeaders = data(cntr,3:10);
        case  'StrainID'
            n=1;
            row = cntr;
            tot = size(data,1);
            while row<tot; 
                if isempty(data{row,3})
                    row = row-1;
                    break
                end
                row = row+1;
            end
            ImageVars.StrainID = data(cntr:row,3:10);
    end
end
end

function [monData, mP] = InitializeExperiment
% global monData mP
monData = struct([]);
% First Check for Stage Data file and load it.
ODELAY_StageData = dir('ODELAY_StageData.mat');

if ~isempty(ODELAY_StageData)
    monData(1).stageFile = ODELAY_StageData.name;
    monData(1).expDir = ODELAY_StageData.folder;
else
    error('No StageFile In this Directory');
end
       %% add in Error function 
cd(monData(1).expDir); % go to experiment Directory
fileName = fullfile(monData(1).expDir, monData(1).stageFile);
load(fileName);

wellIdx = sort(mP.wellIdx);
for well = 1:mP.numWells
    if isdir(mP.wellID{wellIdx(well)})
        monData.wellDir(well,1) = mP.wellID(wellIdx(well));
    end
end

monData.wellDir = sort(monData.wellDir);
monData.AnalyzedImages = false(mP.numWells, mP.totIter);

% Generate ODELAY ImageVars variables.
monData.IndexFile = dir('*Index_ODELAYData.mat');
if isempty(monData.IndexFile)
    
    % Check for Experiment Description  
    ExpDiscriptionFile = dir('*ODELAYExpDisc.xlsx');

    ImageVars.DirectoryList = monData.wellDir;
    stagefileInfo = dir('ODELAY_StageData.mat');
    TimeDilute       = stagefileInfo.datenum;
    ImageVars.TimeDilute =  mP.StartODELAY;
    ImageVars.Microscope_ID = 'Microscope';
    ImageVars.CorrSwitch = true;
    ImageVars.TiledRows = 3;
    ImageVars.TiledCols = 3;
    ImageVars.binningRate = 104;
    ImageVars.NumWells = mP.numWells;
    ImageVars.Magnification = 10;...mP.mag;
    ImageVars.PixSize = 6.5;...mP.pixSize;
    ImageVars.CameraRotaionAngle = 0;
    ImageVars.collect_background = false;
    ImageVars.OffsetMult = 1.6;
    ImageVars.RotAngle = 0;
    ImageVars.phaseXYcor = [0,2;
                            1653,3;
                            3131,0;
                            3,1558;
                            1564,1556;
                            3130,1557;
                            1,3126;
                            1566,3120;
                            3133,123];

    ImageVars.StrainIDHeaders = {{'Plate Well'},{'ODELAY Well'},{'Plot Name'},{'Temperature'},{'Base Media'},{'Misc2'}};
    ImageVars.StrainID = cell(mP.numWells,10);
    for n = 1:mP.numWells; ImageVars.StrainID(n,1) =  monData.wellDir(n); end;
    ImageVars.DirectoryList = monData.wellDir;
    ImageVars.ImageDirectory = monData.expDir;
    ImageVars.NumWells = mP.numWells;
    if ~isempty(ExpDiscriptionFile)
        [~, ~, ~,ImageVarsTemp] = ODELAY_importExpDescription(ExpDiscriptionFile(1,:)); 
        ImageVars.StrainID = ImageVarsTemp.StrainID;
    end

    
    if mP.iterNum>1 && ~isfield(mP,'background') && ~isfield(monData, 'background')
        monData.background = ODELAY_MatCollectBackground(monData, mP);
    end
    
    ImageVars.Background = monData.background; 
    % Figure out the number of wells in the experiment
    cd(monData.expDir);
    dirList = dir;
    
    dirList(1:2) = [];
    isDirList = [dirList.isdir]';
    dirList(~isDirList)=[];
    WellData = struct([]);
    numwells =  mP.numWells;
       
    stInd = regexp(monData.expDir, filesep);
    
    [path,ExperimentName,ext]    = fileparts(monData.expDir);
    CentroidAssociationVersion   = 'ODELAY_RemotProcess_v2';
    WellProcessVersion   = 'ODELAY_RemotProcess_v2';
    DataDirectory     = fullfile(path,ExperimentName,'ODELAY Well Data');
    DataDirctoryExist = true;
    ImageDirectory    = monData.expDir;
    CurrentCD         = monData.expDir;
    stagefileInfo     = dir('ODELAY_StageData.mat');
    TimeDilute        =  mP.StartODELAY;
   
    for well = 1:numwells
                cntr = strcmp({dirList(:).name;},monData.wellDir(well,1));
                if sum(well)~=0
                    WellData(well,1).WellID = dirList(cntr).name;
                    WellData(well,1).WellIndex = well;
                end
    end
   
    numTimePoints = mP.totIter;
    maxObj = 5000;
     for well = 1:mP.numWells
        WellOut(well,1).WellID = WellData(well,1).WellID;
        WellOut(well,1).WellIndex = WellData(well,1).WellIndex;
        WellOut(well,1).BinaryData.bwThreshCorr  = zeros(numTimePoints,1);
        WellOut(well,1).BinaryData.bwThresh      = zeros(numTimePoints,1);
        WellOut(well,1).BinaryData.bwThreshRatio = zeros(numTimePoints,1);
        WellOut(well,1).BinaryData.bwMaxoffset   = zeros(numTimePoints,1);
        WellOut(well,1).BinaryData.bwMinoffset   = zeros(numTimePoints,1);
        WellOut(well,1).BinaryData.bwArea        = zeros(numTimePoints,1);
        
        WellOut(well,1).xyDisp        = NaN(numTimePoints,2);
        WellOut(well,1).ClockPoints = NaN(numTimePoints,1);
        WellOut(well,1).TimePoints = NaN(numTimePoints,1);
        
        WellOut(well,1).ImageIndex  = zeros(numTimePoints,1);
        WellOut(well,1).ObjectNext  = NaN(maxObj, numTimePoints);
        WellOut(well,1).ObjectTrack  = NaN(maxObj, numTimePoints);
        WellOut(well,1).ObjectArea  = NaN(maxObj, numTimePoints);
        WellOut(well,1).ObjectCentX = NaN(maxObj, numTimePoints);
        WellOut(well,1).ObjectCentY = NaN(maxObj, numTimePoints);
        WellOut(well,1).numObj = NaN(numTimePoints,1);
        WellOut(well,1).trackLength = zeros(maxObj,1);
        Tracks2(well,1).ObjectInfo.WellID        = WellData(well,1).WellID;
        Tracks2(well,1).ObjectInfo.FitDataGomp   = NaN(maxObj,17);
        Tracks2(well,1).ObjectInfo.FitDataGompDT   = NaN(maxObj,17);
        Tracks2(well,1).ObjectInfo.FitDataLinear = NaN(maxObj,5);
        Tracks2(well,1).ObjectInfo.TimePoints    = NaN(numTimePoints, 1);
        Tracks2(well,1).ObjectInfo.TimePointsl   = false(numTimePoints,1);
        Tracks2(well,1).ObjectInfo.ImageIndex    = NaN(numTimePoints, 1);
        Tracks2(well,1).ObjectInfo.ObjectArea    = NaN(maxObj, numTimePoints);
        Tracks2(well,1).ObjectInfo.flagIndex     = false(maxObj,1);
        Tracks2(well,1).ObjectInfo.ObjectTrack   = NaN(maxObj, numTimePoints);
        Tracks2(well,1).ObjectInfo.ObjectTrack(:,1) = (1:maxObj);
     end
    % Add book keeping variable to latest image in each well stack that has
    % been processed.  
    monData.numWells = numwells;
    monData.StartODELAY = mP.StartODELAY;
    monData.AnalizedImages = zeros(numwells, numTimePoints); 
    monData.IndexFile = dir('*Index_ODELAYData.mat');
    monData.magnification = ImageVars.Magnification;
    monData.pixSize = ImageVars.PixSize;
    ODELAYSaveData % Save the Data using ODELAYSaveData
    
else

    load(fullfile(monData.IndexFile.folder,monData.IndexFile.name), 'monData');
end


end









