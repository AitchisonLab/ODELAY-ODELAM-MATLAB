function ODELAY_MicroscopeControl
clear global
clear all
global mmc mP hCont hImageDisp hNav 

 
 %% Initialize Microscope
ActivateMicroscope;

%% Initialize Microscope Properites 
InitializeMicroscopeProperties;  %***Generates mP Structure to store Figure and Dataset Props***
    
%% Initialize Microscope Control Gui
hCont      = InitializeControlGui;           % ***Generates Figure with dependencies on mP and mcc***   
hImageDisp = InitializeImageDisplayGui;      %  ***Separate File***




mP.hVideoTimer = timer('BusyMode','drop',...
                       'ExecutionMode','fixedSpacing',...
                       'StartFcn',@VideoInit,...
                       'StopFcn', @VideoStop,...
                       'TimerFcn', @VideoUpdate,...
                       'Period',0.05);
                  
mP.hODELAYTimer = timer('BusyMode','queue',...
                        'ExecutionMode','fixedRate',...
                        'StartFcn',@ODELAYStart,...
                        'StopFcn', @ODELAYStop,...
                        'TimerFcn',@ODELAYScanPlate,...
                        'ErrorFcn',@ODELAYError,...
                        'TasksToExecute',mP.totIter,...
                        'Period',mP.iterPeriod);
                        
%% Microscope Controls GUI Handles List 
%               hControlPanel: 1
%             hResetFunctions: 173.0128
%              hStagePosPanel: 0.0129
%                hFocusContAx: 2.0129
%                hStageContAx: 3.0129
%                  hNavContAx: 4.0129
%                 FocusContIm: 5.0129
%                hStageContIm: 6.0134
%                  hNavContIm: 7.0134
%                      htextX: 8.0134
%                  hStagePosX: 9.0129
%                      htextY: 10.0129
%                  hStagePosY: 11.0129
%                      htextZ: 12.0129
%                  hStagePosZ: 13.0129
%                   htextWell: 14.0129
%               hStagePosWell: 15.0129
%                htextOriginX: 16.0129
%               hStageOriginX: 17.0129
%                htextOriginY: 174.0128
%               hStageOriginY: 175.0128
%                htextOriginZ: 176.0128
%               hStageOriginZ: 177.0128
%             hStageOriginSet: 178.0128
%              hStageGoOrigin: 179.0128
%                  hODELAYDir: 180.0128
%                  hRunODELAY: 181.0128
%                 hIllumPanel: 182.0128
%                  hTransLamp: 184.0128
%              hTransLampEdit: 185.0128
%                   hTransApp: 186.0128
%               hTransAppEdit: 187.0128
%                   hFluorApp: 188.0128
%               hFluorAppEdit: 189.0128
%               hTransShutter: 190.0128
%               hFluorShutter: 191.0128
%                  hCubePanel: 192.0128
%                    hCubeSel: [194.0128 197.0128 200.0128 203.0128 206.0128 209.0128 212.0128]
%                  hCubeFocus: [195.0128 198.0128 201.0128 204.0128 207.0128 210.0128 213.0128]
%                    hCubeExp: [196.0128 199.0128 202.0128 205.0128 208.0128 211.0128 214.0128]
%                hCameraPanel: 215.0128
%                 hCameraSnap: 217.0128
%                hCameraFocus: 218.0128
%             hAutofocusPanel: 219.0128
%            hAutofocusButton: 221.0128
%            hAutofocusSelect: 222.0128
%             hAutofocusRange: 223.0128
%         hAutofocusRangeText: 224.0128
%             hAutofocusSteps: 225.0128
%         hAutofocusStepsText: 226.0128
%         hAutofocusTargetInc: 227.0128
%     hAutofocusTargetIncText: 228.0128

%% Microscope Handle CallBacks
set(hCont.hControlPanel,  'CloseRequestFcn',@CloseMicroscopeControl);
set(hCont.hStageContIm,   'ButtonDownFcn',@guiMoveStage);
set(hCont.FocusContIm ,   'ButtonDownFcn',@guiFocusCont);
set(hCont.hNavContAx,     'ButtonDownFcn',@guiPlateNav);
set(hCont.hTransLamp,     'Callback',@guiTransLamp)
set(hCont.hTransLampEdit, 'Callback',@guiTransLampEdit) 
set(hCont.hTransApp,      'Callback',@guiTransApp)
set(hCont.hFluorApp,      'Callback', @guiFluorApp) 
set(hCont.hTransShutter,  'Callback', @guiTransShutter)
set(hCont.hFluorShutter,  'Callback', @guiFluorShutter)

set(hCont.hCubeSel,       'Callback',@guiCubeSelect);
set(hCont.hRecordSel,     'Callback',@guiRecSelect);

set(hCont.hCameraSnap,       'Callback', @guiCameraSnap);
set(hCont.hCameraFocus,      'Callback', @guiCameraFocus); 
set(hCont.hResetFunctions,   'Callback', @ResetFunctions);
set(hCont.hAutofocusButton,  'Callback',@guiAutofocus);
set(hCont.hAutofocusButton2, 'Callback',@guiAutofocus2)
set(hImageDisp.hImagePanel,  'CloseRequestFcn',@CloseMicroscopeControl);
set(hCont.hStageOriginSet,   'Callback',@guiSetOrigin);
set(hCont.hStageGoOrigin,    'Callback',@guiGoOrigin);

set(hCont.hStagePosX, 'String',num2str(mmc.getXPosition(mP.xyDrive),'%6.2f')); %hCont.hStagePosX
set(hCont.hStagePosY, 'String',num2str(mmc.getYPosition(mP.xyDrive),'%6.2f'));% hCont.hStagePosY 
set(hCont.hStagePosZ, 'String',num2str(mmc.getPosition(mP.zDrive),'%6.2f'));% hCont.hStagePosZ 

set(hCont.hSlideSelect, 'CallBack',@guiSlideSelect);
set(hCont.hRunODELAY, 'CallBack', @guiRunODELAY);

set(hCont.hAutofocusButton2, 'Callback',@guiAutofocus2)
set(hCont.hPosList,   'Callback', {@guiUpDateListPos,0})
set(hCont.hLoadList,  'Callback',  @guiLoadPosList)
set(hCont.hUpdatePos, 'Callback',  @guiUpdatePos)
set(hCont.hMovePrev,  'Callback', {@guiUpDateListPos,1})
set(hCont.hMovNext,   'Callback', {@guiUpDateListPos,2})

set(hCont.hNumSlices,   'Callback',{@guiSliceUpdate});
set(hCont.hSlicesThick, 'Callback',{@guiSliceUpdate});   
                               
hCont.hTotalThickness   
drawnow;
end

function ResetFunctions(hFig, events)
global mmc mP hCont hImageDisp hNav
%%
%               hControlPanel: 1
%             hResetFunctions: 173.0128
%              hStagePosPanel: 0.0129
%                hFocusContAx: 2.0129
%                hStageContAx: 3.0129
%                  hNavContAx: 4.0129
%                 FocusContIm: 5.0129
%                hStageContIm: 6.0134
%                  hNavContIm: 7.0134
%                      htextX: 8.0134
%                  hStagePosX: 9.0129
%                      htextY: 10.0129
%                  hStagePosY: 11.0129
%                      htextZ: 12.0129
%                  hStagePosZ: 13.0129
%                   htextWell: 14.0129
%               hStagePosWell: 15.0129
%                htextOriginX: 16.0129
%               hStageOriginX: 17.0129
%                htextOriginY: 174.0128
%               hStageOriginY: 175.0128
%                htextOriginZ: 176.0128
%               hStageOriginZ: 177.0128

%                  hODELAYDir: 180.0128

%                 hIllumPanel: 182.0128
%                  hTransLamp: 184.0128
%              hTransLampEdit: 185.0128
%                   hTransApp: 186.0128
%               hTransAppEdit: 187.0128
%                   hFluorApp: 188.0128
%               hFluorAppEdit: 189.0128
%               hTransShutter: 190.0128
%               hFluorShutter: 191.0128
%                  hCubePanel: 192.0128
%                    hCubeSel: [194.0128 197.0128 200.0128 203.0128 206.0128 209.0128 212.0128]
%                  hCubeFocus: [195.0128 198.0128 201.0128 204.0128 207.0128 210.0128 213.0128]
%                    hCubeExp: [196.0128 199.0128 202.0128 205.0128 208.0128 211.0128 214.0128]
%                hCameraPanel: 215.0128
%                 hCameraSnap: 217.0128
%                hCameraFocus: 218.0128
%             hAutofocusPanel: 219.0128
%            hAutofocusButton: 221.0128
%            hAutofocusSelect: 222.0128
%             hAutofocusRange: 223.0128
%         hAutofocusRangeText: 224.0128
%             hAutofocusSteps: 225.0128
%         hAutofocusStepsText: 226.0128
%         hAutofocusTargetInc: 227.0128
%     hAutofocusTargetIncText: 228.0128
% lampVal = mmc.getProperty('Transmitted Light','Level');
% 
%  mmc.defineConfig('ImageMode', 'BrightField', 'Transmitted Light','Level',lampVal);
% mP         = InitializeMicroscopeProperties;

%   set(hCont.hTransLamp,   'value',15*2000/20,...
%                              'min',0,'max',2000,...
%                              'SliderStep',[1/2000 10/2000])
% mP.tilePos    = [mP.tileOrder(:,1).*(mP.sensorSize(2)*mP.pixSize/mP.mag*(1-mP.overlap)),...
%                  mP.tileOrder(:,2).*(mP.sensorSize(1)*mP.pixSize/mP.mag*(1-mP.overlap))];
% 
% scanWell(1,1);
% scanWell(4,1);
% scanWell(3,1);
% scanWell(13,1);
%  correctZPos(1);
% InitializeMicroscopeProperties; 
% scanWell(20, 1);
% fluorShutter = mmc.getProperty('IntensiLightShutter','State');
% diaShutter = mmc.getProperty('TIDiaLamp','State');
% restartFlag=false;
%     if strcmp(get(mP.hVideoTimer, 'Running'),'on'); stop(mP.hVideoTimer); restartFlag = true; end;
%      mmc.snapImage();
%     tempImage = mmc.getImage();
%     firstImage = reshape(typecast(tempImage,'uint16'),1392,1040)';
%     
% focusMethods =  get(hCont.hAutofocusSelect,'String');
% AutoFunc = focusMethods(get(hCont.hAutofocusSelect,'Value'),:);
% zFocusPos = [-20:0.5:20]+mmc.getPosition(mP.zDrive);
% zFocus  = NaN(numel(zFocusPos),size(focusMethods,1)+1);
% zFocus(:,1) = zFocusPos';
% zFocus(:,2) = NaN(numel(zFocusPos),1);
% nSteps = numel(zFocusPos);
% zPos = mmc.getPosition(mP.zDrive);
% zRecord = NaN(2*nSteps+21,3);
% 
% for step = 1:nSteps
%     mmc.setPosition(mP.zDrive, zFocus(step,1))
%     mmc.waitForSystem(); pause(0.05);
%     mmc.snapImage();
%     tempImage = mmc.getImage();
%     cameraImage = reshape(typecast(tempImage,'uint16'),1392,1040)';
%     zFocus(step,1) = mmc.getPosition(mP.zDrive);
%     for n = [1,2,4,7:25]
%         fprintf([focusMethods(n,:), '\n']);
%         zFocus(step,n+1) = fmeasure(cameraImage, focusMethods(n,:), []);
%     end
%     [yout,x]=imhist(cameraImage, double(2^12));
%     set(hImageDisp.hPlot,'XData',x,'YData',yout);
%     set(hImageDisp.hImDisplay,'CData',cameraImage) 
% %     zRecord(recInd,1:2) = zFocus(step,1:2);
% %     recInd = recInd+1;
%     drawnow;
% end
% 
% assignin('base', 'zFocus',zFocus);
% assignin('base', 'focusMethods',focusMethods);
% assignin('base', 'firstImage',firstImage);

% test = get(mP.hODELAYTimer,'TasksExecuted')
% 
% load('ODELAY_StageData.mat'); 
% test = get(mP.hODELAYTimer,'TasksExecuted');
test = 1;
% scanWell(1, 1)

% % 
% set(mP.hODELAYTimer, 'TasksToExecute',280);
% test = get(mP.hODELAYTimer,'TasksToExecute');
% mP.iterPeriod=30*60;
% set(mP.hODELAYTimer,'Period',mP.iterPeriod);
% set(hCont.hNavContAx,     'ButtonDownFcn',@guiPlateNav);
%
mmc.defineConfig('ImageMode', 'BrightField', 'Core','AutoShutter','0');
mmc.defineConfig('ImageMode', 'BrightField', 'IntensiLightShutter','State', '0');
drawnow;
end

%% Callback Functions

function guiMoveStage(hFig, events)
global mmc mP hCont 

curPoint = get(hCont.hStageContAx,'CurrentPoint');
    XLim = get(hCont.hStageContAx,'XLim');
    YLim = get(hCont.hStageContAx,'YLim');
  axOrigin = [diff(XLim)/2, diff(YLim)/2];
  dispXY = 100.*(axOrigin-curPoint(1,1:2))./axOrigin;
  xPos = mmc.getXPosition(mP.xyDrive);
  yPos = mmc.getYPosition(mP.xyDrive);
  
  newPos = [xPos,yPos] - dispXY;
  
  mmc.setXYPosition(mP.xyDrive,newPos(1), newPos(2));

end

function guiFocusCont(hFig, events)
global mmc mP hCont 

curPoint = get(hCont.hFocusContAx,'CurrentPoint');
    YLim = get(hCont.hFocusContAx,'YLim');
  axOrigin = diff(YLim)/2;
  dispZ = 10.*(axOrigin-curPoint(1,2))./axOrigin;
  
  zPos = mmc.getPosition(mP.zDrive);
  newZPos = zPos+dispZ;
  
  mmc.setPosition(mP.zDrive,newZPos);

end

function guiPlateNav(hFig, events)
global mmc mP hCont hImageDisp
    % Get pointer Positions within Axes
    pointerpos = get(hCont.hNavContAx, 'CurrentPoint');
    x = pointerpos(1,1);
    y = pointerpos(1,2);

    cCenters = mP.cCenters(:,[1,2])+mP.cCenters(:,[3,4])./2;
    pointDist =  sqrt((cCenters(:,1)-x).^2 +(cCenters(:,2)-y).^2);
    [val,ind] = min(pointDist);         

    indx = false(numel(hCont.hNavCircles),1);
    indx(ind) = true;

    
%     cColor = get(hCont.hNavCircles(indx),'FaceColor');
    set(hCont.hNavCircles( indx), 'FaceColor',[1 0 0]);
    set(hCont.hNavCircles(~indx), 'FaceColor','none');
    set(hCont.hPosList, 'Value',ind);
    
    xNavPos = mP.stageXYPos(indx,1);
    yNavPos = mP.stageXYPos(indx,2);

    xAbsPos = xNavPos;...+mP.XYZOrigin(1);
    yAbsPos = yNavPos;...+mP.XYZOrigin(2);
    [xAbsPos,yAbsPos]
    set(hCont.hStagePosWell,'String',mP.wellID{indx});
    mmc.setXYPosition( mP.xyDrive,xAbsPos,yAbsPos);
    
end

function guiSlideSelect(hFig,events)
global mP hCont
    listpos = get(hFig, 'Value');
    delete(hCont.hNavCircles)
    
    GenerateSpotPositions(listpos)
    
    switch listpos
        case 1
            singCondIm = imread('Drug 1 chamber Slide Gasket MakerDesign v1.jpg');
            set(hCont.hNavContIm, 'CData',singCondIm);
            
            pixum = 690/(4500*13);
            imDim = size(singCondIm);
            xMax =  imDim(2)/pixum;
            yMax =  imDim(1)/pixum;
            xOrigin = 99/pixum;
            yOrigin = 110/pixum;
            set(hCont.hNavContAx,'XTick',[],'YTick',[],'PickableParts','all','HitTest','on','YDir','reverse',...
                                  'XLim',[0,xMax]-xOrigin, 'YLim',[0,yMax]-yOrigin);

            cols = 0:4500:4500*13;
            rows = 0:4500:4500*7;
            spacing =  4500;
            xoffset = -2250;
            yoffset = -2250;
            radius = spacing;
            
        case 2
            fiveCondIm = imread('Drug 5 chamber Slide Gasket MakerDesign v3a.jpg');
            set(hCont.hNavContIm, 'CData',fiveCondIm);
            
            pixum = 690/(4500*13);
            imDim = size(fiveCondIm);
            xMax =  imDim(2)/pixum;
            yMax =  imDim(1)/pixum;
            xOrigin = 99/pixum;
            yOrigin = 110/pixum;
            set(hCont.hNavContAx,'XTick',[],'YTick',[],'PickableParts','all','HitTest','on','YDir','reverse',...
                                  'XLim',[0,xMax]-xOrigin, 'YLim',[0,yMax]-yOrigin);

            cols = 0:4500:4500*13;
            rows = 0:4500:4500*7;
            spacing =  4500;
            xoffset = -2250;
            yoffset = -2250;
            radius = spacing;
            
         case 3
            twentyfourCondIm = imread('24_Well_Plate.jpg');
            twentyfour = zeros(594,889,3,'uint8');
            twentyfour(:,:,1) = imresize(twentyfourCondIm,[594,889]);
            twentyfour(:,:,2) = imresize(twentyfourCondIm,[594,889]);
            twentyfour(:,:,3) = imresize(twentyfourCondIm,[594,889]);
            set(hCont.hNavContIm, 'CData',twentyfour);
            
            pixum = 670/(19500*5);
            imDim = size(twentyfour);
            xMax =  imDim(2)/pixum;
            yMax =  imDim(1)/pixum;
            xOrigin = 134/pixum;
            yOrigin = 130/pixum;
            set(hCont.hNavContAx,'XTick',[],'YTick',[],'PickableParts','all','HitTest','on','YDir','reverse',...
                                  'XLim',[0,xMax]-xOrigin, 'YLim',[0,yMax]-yOrigin);
                                           
            cols = 0:19500:19500*5;
            rows = 0:19500:19500*3;
            spacing =  19500;
            xoffset = -16280/2;
            yoffset = -16280/2;
            radius =  16280;
            
    end
    

hCont.hNavCircles = gobjects(numel(rows)*numel(cols),1);

cnt = 1;
for row = rows
    for col = cols

        cCenters(cnt,:) = [col+xoffset, row+yoffset, radius,radius];
        hCont.hNavCircles(cnt) = rectangle('Position',cCenters(cnt,:),...
                                               'EdgeColor',[1 0 0],...
                                               'Curvature',[1,1],...
                                               'Parent',hCont.hNavContAx,...
                                               'HitTest','off');
                                           
       

        cnt = cnt+1;
    end
  
end

mP.cCenters = cCenters;
    
end

function guiCubeSelect(hFig,events)
global mmc mP hCont 
    restartFlag=false;
    if strcmp(get(mP.hVideoTimer, 'Running'),'on'); stop(mP.hVideoTimer); restartFlag = true; end;
    if mmc.isSequenceRunning();   mmc.stopSequenceAcquisition();  end
    
% Figure out which Cube button is depressed
    cubeInd = hCont.hCubeSel == hFig;
    numCubes = numel(hCont.hCubeSel);
    cubePos = 0:numCubes-1;
    for n = 1:numCubes
        set(hCont.hCubeSel(n), 'value', cubeInd(n))
    end

%One Slot is used for Brightfield...this could be any slot need to get ID
%from mP structure!!!!!!!!!

if cubePos(cubeInd)== cubePos(numCubes)
  mmc.setConfig('ImageMode','BrightField');
else
%   mmc.setProperty('IntensiLightShutter','State',1);
%   mmc.setProperty('TIEpiShutter','State',0);
%   set(hCont.hFluorShutter, 'Value',0);
%   set(hCont.hTransShutter, 'Value',1);
%   mmc.setProperty('TIFilterBlock1','State', cubePos(cubeInd));
  mmc.setConfig('ImageMode', mP.stateNames{cubeInd});
end
% lampLevel = mmc.getProperty('Transmitted Light','Level');
% set(hCont.hTransLamp,'Value',str2double(lampLevel));

    mP.CubeEngaged = cubeInd;
    if restartFlag; start(mP.hVideoTimer); end
end

function guiRecSelect(hFig,events)
global mmc mP hCont 
   
% Figure out which Cube button is depressed
    cubeInd = hCont.hRecordSel == hFig;
    val = get(hFig,'value');
    if val
        set(hFig,'BackgroundColor',[0 1 0]);
        mP.RecChan(cubeInd)  = val;
    else
        set(hFig,'BackgroundColor',[1 0 0]);
        mP.RecChan(cubeInd)  = val;
    end
    
  
end

function updateGUIXYZ
global mmc mP hCont

% fluorShutter = mmc.getProperty('IntensiLightShutter','State')=='1';
% transShutter = mmc.getProperty('TIEpiShutter','State')=='1';

xAbsPos = mmc.getXPosition(mP.xyDrive);
yAbsPos = mmc.getYPosition(mP.xyDrive);
zAbsPos = mmc.getPosition(mP.zDrive);

xPos = xAbsPos-mP.XYZOrigin(1);
yPos = yAbsPos-mP.XYZOrigin(2);
zPos = zAbsPos-mP.XYZOrigin(3);

set(hCont.hStagePosX,'String',num2str(xPos,'%6.1f')); %hCont.hStagePosX
set(hCont.hStagePosY,'String',num2str(yPos,'%6.1f'));% hCont.hStagePosY 
set(hCont.hStagePosZ, 'String',num2str(zPos,'%6.1f'));% hCont.hStagePosZ

% set(hCont.hTransShutter,'Value',transShutter);
% set(hCont.hFluorShutter,'Value',fluorShutter);

end

function guiTransLamp(hFig,events)
global mmc mP hCont 
    val = get(hCont.hTransLamp, 'Value');
    if val<0
        val = 0;
    elseif val>2000
        val = 2000;
    end
    
try
powerSup = serial('COM4',...
       'BaudRate',9600,...
       'DataBits', 8,...
       'Parity','none',...
       'StopBits',1);
fopen(powerSup);
vSet = num2str(val/100,'%05.2f');
fprintf(powerSup, '%s', ['VSET1:',vSet]);
pause(0.1);
fprintf(powerSup, '%s','VOUT1?'); 
out = fread(powerSup,5); 
fclose(powerSup)
delete(powerSup)
clear powerSup
set(hCont.hTransLampEdit, 'String',char(out(1:5)'));
catch
end
%     if mmc.getProperty('Transmitted Light','State')=='0'; mmc.setProperty('Transmitted Light','State',1);end;
    
%     mmc.setProperty('Transmitted Light','Level', round(val));
%     set(hCont.hTransLamp, 'Value',val);
 

end

function guiTransLampEdit(hFig,events)
global mmc mP hCont 
    strVal = get(hCont.hTransLampEdit, 'String');
    val = round(str2double(strVal));
    if val<0; val = 0; elseif val>255; val = 255; end
    mmc.setProperty('Transmitted Light','Level', round(val));
    set(hCont.hTransLamp, 'Value',val);
    set(hCont.hTransLampEdit, 'String',num2str(round(val),'%i'));

end

function guiTransApp(hFig, events)
global mmc mP hCont
    % Set min and max values from mP and Microscope!!!!
    val = get(hCont.hTransApp, 'Value');
    minVal = 1; maxVal = 46;
    if val<minVal; val=minval; elseif val>maxVal; val = maxVal; end
    set(hCont.hTransApp,'Value',val);
    set(hCont.hTransAppEdit,'String',num2str(round(val),'%i'));
%     mmc.setProperty('TL-FieldDiaphragm', 'Position', round(val));
end

function guiFluorApp(hFig, events)
global mmc mP hCont
    % Set min and max values from mP and Microscope!!!!
    val = get(hCont.hFluorApp, 'Value');
    minVal = 1; maxVal = 12;
    if val<minVal; val=minval; elseif val>maxVal; val = maxVal; end
    set(hCont.hFluorApp,'Value',round(val));
    set(hCont.hFluorAppEdit,'String',num2str(round(val),'%i'));
%     mmc.setProperty('IL-FieldDiaphragm', 'Position', round(val));
end

function guiTransShutter(hFig, events)
global mmc mP hCont
    val = get(hCont.hTransShutter, 'Value');
    if val~=1; val = 0; else val = 1; end
    set(hCont.hTransShutter,'Value',val);
    mmc.setProperty('TIDiaShutter', 'State', val);
end

function guiFluorShutter(hFig, events)
global mmc mP hCont
    val = get(hCont.hFluorShutter, 'Value');
    if val~=1; val = 0; else val = 1; end
    set(hCont.hFluorShutter,'Value',val);
    mmc.setProperty(mP.fluorShutter, 'State', val);
end

function XYPosUpdate(hTimer, events)
global mmc mP hCont

    set(hCont.hStagePosX,'String', num2str(mmc.getXPosition(mP.xyDrive),'%6.2f')); 
    set(hCont.hStagePosY,'String', num2str(mmc.getYPosition(mP.xyDrive),'%6.2f'));
    set(hCont.hStagePosZ,'String', num2str(mmc.getPosition(mP.zDrive),  '%6.2f'));

end

function guiCameraSnap(hFig, events)
global mmc mP hCont hImageDisp

    if strcmp(get(mP.hVideoTimer, 'Running'),'Running'); stop(mP.hVideoTimer); end;
    if mmc.isSequenceRunning();   mmc.stopSequenceAcquisition();  end
    [~,ind] = find(mP.CubeEngaged==1);
    expTime = str2double(get(hCont.hCubeExp(ind),'String'));
    mmc.setExposure(expTime);
    mmc.snapImage
    cameraImage = reshape(typecast(mmc.getImage,'uint16'),1392,1040)';
    set(hImageDisp.hImDisplay,'CData',cameraImage);
    drawnow;
    
end

function guiCameraFocus(hFig, events)
global mP mmc hCont hImageDis 
    Val = get(hFig,'Value');
   maxVal  = get(hFig,'Max'); 
   minVal  = get(hFig,'Min'); 
   
%    [~,ind] = find(mP.CubeEngaged==1);
%    expTime = str2double(get(hCont.hCubeFocus(ind),'String'));  
%     mmc.setExposure(expTime);

    if Val==maxVal
        start(mP.hVideoTimer)
    elseif Val==minVal;
        stop(mP.hVideoTimer)
    end
end

function guiSetOrigin(hFig, events)
global mmc mP hCont hImageDisp

xAbsPos = mmc.getXPosition(mP.xyDrive);
yAbsPos = mmc.getYPosition(mP.xyDrive);
zAbsPos = mmc.getPosition(mP.zDrive);

mP.XYZOrigin = [xAbsPos, yAbsPos, zAbsPos];

rowVec = 1:numel(mP.rowName);
colVec = 1:numel(mP.colName);
cnt = 1;
for row = rowVec;
    for col = colVec;
       mP.stageXYPos(cnt,:) = [mP.xgrid(col)+mP.XYZOrigin(1), mP.ygrid(row)+mP.XYZOrigin(2)];
       mP.stageZPos(cnt,mP.iterNum)  = mP.XYZOrigin(3); %Calculated stage positions before focus
       cnt=cnt+1;
    end
end
mP.twoPhaseFocus = true;


set(hCont.hStageOriginX,'String',num2str(mP.XYZOrigin(1), '%7.1f'));
set(hCont.hStageOriginY,'String',num2str(mP.XYZOrigin(2), '%7.1f'));
set(hCont.hStageOriginZ,'String',num2str(mP.XYZOrigin(3), '%7.1f'));

xPos = xAbsPos-mP.XYZOrigin(1);
yPos = yAbsPos-mP.XYZOrigin(2);
zPos = zAbsPos-mP.XYZOrigin(3);

set(hCont.hStagePosX,'String',num2str(xPos, '%7.1f'));
set(hCont.hStagePosY,'String',num2str(yPos, '%7.1f'));
set(hCont.hStagePosZ,'String',num2str(zPos, '%7.1f'));

end

function guiGoOrigin(hFig,events)
global mmc mP hCont hImageDisp

mmc.setXYPosition(mP.xyDrive, mP.XYZOrigin(1), mP.XYZOrigin(2));
mmc.setPosition(mP.zDrive, mP.XYZOrigin(3));

end

function guiRunODELAY(hFig,events)
global mP hCont

listpos = get(hCont.hSlideSelect, 'Value');

GenerateSpotPositions(listpos);
if ~isfield(mP, 'StartODELAY')
    mP.StartODELAY = now;
end
if strcmp(get(mP.hVideoTimer, 'Running'),'on'); stop(mP.hVideoTimer); end
if strcmp(get(mP.hODELAYTimer, 'Running'),'on')
    stop(mP.hODELAYTimer); 
elseif strcmp(get(mP.hODELAYTimer, 'Running'),'off')
    start(mP.hODELAYTimer); 
end

end

function guiAutofocus(hFig, events)
global hCont
    focusMethods =  get(hCont.hAutofocusSelect,'String');
    AutoFunc = focusMethods(get(hCont.hAutofocusSelect,'Value'),:);
    zRange = str2double(get(hCont.hAutofocusRange,'String'));
    nSteps = str2double(get(hCont.hAutofocusSteps,'String'));
    targetIncrement = str2double(get(hCont.hAutofocusTargetInc,'String'));
    [zFocus, ~] = Autofocus(AutoFunc, zRange, nSteps, targetIncrement);

end

function guiAutofocus2(hFig, events)
global hCont
    focusMethods =  get(hCont.hAutofocusSelect2,'String');
    AutoFunc = focusMethods(get(hCont.hAutofocusSelect2,'Value'),:);
    zRange = str2double(get(hCont.hAutofocusRange2,'String'));
    nSteps = str2double(get(hCont.hAutofocusSteps2,'String'));
    targetIncrement = str2double(get(hCont.hAutofocusTargetInc2,'String'));
    [zFocus, ~] = Autofocus(AutoFunc, zRange, nSteps, targetIncrement);
end

function guiUpDateListPos(hFig,events, eventcase)
global mmc mP hCont
    listpos = get(hCont.hPosList, 'Value');
    indx = false(numel(hCont.hNavCircles),1);
    switch eventcase
        case 0
          indx(listpos) = true;  
        
        case 1
            if listpos<=1;
                 listpos = mP.numWells;
            else
                listpos = listpos-1;
            end
            
        case 2
             if listpos>=mP.numWells;
                listpos = mP.numWells;
             else
                 listpos = listpos+1;
               
            end
    end
    indx(listpos) = true;
    set(hCont.hNavCircles( indx), 'FaceColor',[1 0 0]);
    set(hCont.hNavCircles(~indx), 'FaceColor','none');
    set(hCont.hPosList, 'Value',listpos);

    xNavPos = mP.stageXYPos(indx,1);
    yNavPos = mP.stageXYPos(indx,2);

    xAbsPos = xNavPos+mP.XYZOrigin(1);
    yAbsPos = yNavPos+mP.XYZOrigin(2);
    
    set(hCont.hStagePosWell,'String',mP.wellID{indx});
    mmc.setXYPosition(mP.xyDrive,xAbsPos,yAbsPos);        
    
end

function guiSliceUpdate(hFig, events)
global hCont mP

numSlices = str2double(get(hCont.hNumSlices,   'String'));
if mod(numSlices,2)==1 
    mP.numSlice = abs(numSlices);
end

mP.sliceThick = str2double(get(hCont.hSlicesThick,'String'));   

if mP.numSlice >1 
   mP.stackVec = -floor(mP.numSlice/2)*mP.sliceThick:mP.sliceThick:floor(mP.numSlice/2)*mP.sliceThick;
else
    mP.stackVec = 0;
end

set(hCont.hNumSlices,   'String', num2str(mP.numSlice));
% set(hCont
set(hCont.hTotalThickness, 'String',num2str(2*mP.stackVec(end)));   
end

%% To Do Include in guiRunODELAY...need to account for 
function GenerateSpotPositions(cond)

global mP mmc hCont

mP.numTiles = 9;

magInd = mmc.getProperty(mP.objectiveTurret,'State');

mP.mag = mP.objectivLens(str2double(magInd));
set(hCont.hLensMenu, 'Value',str2double(magInd));

mP.sensorSize = [1040, 1392];
mP.overlap = 0.15;
mP.pixSize = 6.45;
mP.tileOrder = [ 1, -1; 
                 0, -1;
                -1, -1; 
                 1,  0; 
                 0,  0;
                -1,  0;
                 1,  1;
                 0,  1;
                -1,  1];

mP.tilePos    = [mP.tileOrder(:,1).*(mP.sensorSize(2)*mP.pixSize/mP.mag*(1-mP.overlap)),...
                 mP.tileOrder(:,2).*(mP.sensorSize(1)*mP.pixSize/mP.mag*(1-mP.overlap))];

    switch cond
        case 1
            mP.rowName = {'E';'F';'G';'H';'I';'J';'K';'L'};
            mP.colName = {'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';};

            mP.xgrid = -[0:4500:4500*13];
            mP.ygrid =   0:4500:4500*7;

            mP.numWells = numel(mP.rowName)*numel(mP.colName);

            mP.stageXYPos = zeros(mP.numWells, 2);
            mP.stageZPos  = zeros(mP.numWells, mP.numTimePoints+1); %Calculated stage positions before focus
            mP.zFocusPos  = zeros(mP.numWells, mP.numTimePoints+1); %Stage Positions recorded after focus
            mP.wellID =cell(mP.numWells,1);;
            mP.lastImaged = zeros(mP.numWells,1);
%             mP.XYZOrigin(1,:) = [29500,-15500,8500.0]; 
            rowVec = 1:numel(mP.rowName);
            colVec = 1:numel(mP.colName);
            cnt = 1;
            for row = rowVec;
                for col = colVec;
                   mP.stageXYPos(cnt,:)   = [mP.xgrid(col)+mP.XYZOrigin(1), mP.ygrid(row)+mP.XYZOrigin(2)];
                   mP.stageZPos(cnt,1)  = mP.XYZOrigin(3); %Calculated stage positions before focus
                   mP.zFocusPos(cnt,1)  = mP.XYZOrigin(3); %Stage Positions recorded after focus
                   mP.wellID(cnt) = {[mP.rowName{row} mP.colName{col}]};
                   cnt=cnt+1;
                end
            end

            mP.rowVec = 1:2:8;
            mP.colVec = 2:13;   
            
        case 2
            
            mP.rowName = {'E';'F';'G';'H';'I';'J';'K';'L'};
            mP.colName = {'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';};

            mP.xgrid = -[0:4500:4500*13];
            mP.ygrid = 0:4500:4500*7;

            mP.numWells = numel(mP.rowName)*numel(mP.colName);

            mP.stageXYPos = zeros(mP.numWells, 2);
            mP.stageZPos  = zeros(mP.numWells, mP.numTimePoints+1); %Calculated stage positions before focus
            mP.zFocusPos  = zeros(mP.numWells, mP.numTimePoints+1); %Stage Positions recorded after focus
            mP.wellID =cell(mP.numWells,1);;
            mP.lastImaged = zeros(mP.numWells,1);
%             mP.XYZOrigin(1,:) = [29500,-15500,8500.0]; 
            rowVec = 1:numel(mP.rowName);
            colVec = 1:numel(mP.colName);
            cnt = 1;
            for row = rowVec;
                for col = colVec;
                   mP.stageXYPos(cnt,:)   = [mP.xgrid(col)+mP.XYZOrigin(1), mP.ygrid(row)+mP.XYZOrigin(2)];
                   mP.stageZPos(cnt,1)  = mP.XYZOrigin(3); %Calculated stage positions before focus
                   mP.zFocusPos(cnt,1)  = mP.XYZOrigin(3); %Stage Positions recorded after focus
                   mP.wellID(cnt) = {[mP.rowName{row} mP.colName{col}]};
                   cnt=cnt+1;
                end
            end

            mP.rowVec = 1:2:8;
            mP.colVec = [1,2,4,5,7,8,10,11,13,14];     
        case 3
          
            mP.XYZOrigin(1,1:2) = [49000,-30000]; 
            set(hCont.hStageOriginX,'String',num2str(mP.XYZOrigin(1), '%7.1f'));
            set(hCont.hStageOriginY,'String',num2str(mP.XYZOrigin(2), '%7.1f'));
            set(hCont.hStageOriginZ,'String',num2str(mP.XYZOrigin(3), '%7.1f'));
            mP.rowName = {'A';'B';'C';'D'};
            mP.colName = {'01';'02';'03';'04';'05';'06'};

            mP.xgrid = -[0:19500:19500*5];
            mP.ygrid =  [0:19500:19500*3];

            mP.numWells = numel(mP.rowName)*numel(mP.colName);

            mP.stageXYPos = zeros(mP.numWells, 2);
            mP.stageZPos  = zeros(mP.numWells, mP.numTimePoints+1); %Calculated stage positions before focus
            mP.zFocusPos  = zeros(mP.numWells, mP.numTimePoints+1); %Stage Positions recorded after focus
            mP.wellID =cell(mP.numWells,1);
            mP.lastImaged = zeros(mP.numWells,1);
            
            rowVec = 1:numel(mP.rowName);
            colVec = 1:numel(mP.colName);
            cnt = 1;
            for row = rowVec;
                for col = colVec;
                   mP.stageXYPos(cnt,:)   = [mP.xgrid(col)+mP.XYZOrigin(1), mP.ygrid(row)+mP.XYZOrigin(2)];
                   mP.stageZPos(cnt,1)  = mP.XYZOrigin(3); %Calculated stage positions before focus
                   mP.zFocusPos(cnt,1)  = mP.XYZOrigin(3); %Stage Positions recorded after focus
                   mP.wellID(cnt) = {[mP.rowName{row} mP.colName{col}]};
                   cnt=cnt+1;
                end
            end
            
            mP.rowVec = 1:2:4;
            mP.colVec = 1:6;  
    end
    % Set 
    mP.wellIdx = zeros(2*numel(mP.rowVec)*numel(mP.colVec),1);
    mP.numWells = numel(mP.wellIdx);
    for row = mP.rowVec
       sInd = (row-1)*numel(mP.colVec)+1; 
       eInd = sInd+2*numel(mP.colVec)-1;
       mP.wellIdx(sInd:eInd,1) =  [mP.colVec+numel(mP.colName)*(row-1), fliplr(mP.colVec+numel(mP.colName)*(row))];
       mP.rowIdx(sInd:eInd,1) =   [repmat(row, 1,numel(mP.colVec)), repmat(row+1, 1,numel(mP.colVec))];
       mP.colIdx(sInd:eInd,1) =   [mP.colVec, fliplr(mP.colVec)];
    end
    
    set(hCont.hPosList,'String',mP.wellID); 
end

%% Video Timer and Image Callbacks
function VideoInit(hTimer, events)
global mmc mP hCont hImageDisp

    % Check if Sequence is running if so stop the Aquisition
    if mmc.isSequenceRunning();   mmc.stopSequenceAcquisition();  end
%     figure(hImageDisp.hImagePanel);
%     figure(hCont.hControlPanel);
    cubInd = str2num(mmc.getProperty('TIFilterBlock1','State'))+1;
    [~,ind] = find(mP.CubeEngaged==1);
    expTime = str2double(get(hCont.hCubeFocus(ind),'String'));
    mmc.setExposure(expTime);

end

function VideoStop(hTimer, events)
global mmc mP
    mmc.stopSequenceAcquisition();
%     assignin('base','mP',mP);
end
                       
function VideoUpdate(hTimer, events)
global mP mmc hCont hImageDisp

%     bufferDelay = mmc.getBufferIntervalMs;
%     imageRemaining = mmc.getRemainingImageCount();
%     numImages = mmc.getBufferFreeCapacity;
%     if mP.numImages ~= numImages mmc.snapImage();
        mmc.snapImage();
        tempImage = mmc.getImage();
        cameraImage = reshape(typecast(tempImage,'uint16'),1392,1040)';
        [yout,x]=imhist(cameraImage, double(2^12));
        set(hImageDisp.hImDisplay,'CData',cameraImage);
        set(hImageDisp.hPlot,'XData',x,'YData',yout);
        updateGUIXYZ;
%          assignin('base','cameraImage',cameraImage);
%     end
%     mP.numImages = numImages;
   
% mP.step = mP.step+1;

end

function setExposureTime(cubeInd)
global mP mmc hCont    




end

%% Tiff/Image Save Function

%% Directory Parse function

%% AutoFocus version 0.1 Simple iterative focus but might not be sufficent
function [zFocusPos, zRecord] = Autofocus(AutoFunc, zRange, nSteps, targetIncrement)
global mmc mP hCont hImageDisp
    %Stop any other timers currently running to keep them from
    %interfearing.  
    restartFlag=false;
    if strcmp(get(mP.hVideoTimer, 'Running'),'on'); stop(mP.hVideoTimer); restartFlag = true; end

% Calculate Z-Movements
    if isempty(AutoFunc)
        focusMethods =  get(hCont.hAutofocusSelect,'String');
        AutoFunc = focusMethods(get(hCont.hAutofocusSelect,'Value'),:);
    end
    if strcmp(mP.AutoFocusFnc, 'Test')
        
        zFocus = [8571.025;
                  8574.500;
                  8578.025;
                  8581.525;
                  8585.025;
                  8579.775;
                  8583.225;
                  8582.350;
                  8584.100;
                  8582.800;
                  8583.650;
                  8582.975;
                  8583.450;
                  8583.125;
                  8583.325];
        
        zRecord = NaN(2*nSteps+21,3);
           
        for pos = 1:numel(zFocus)
%             fprintf([num2str(zFocus(pos)), '\n']);
            mmc.setPosition(mP.zDrive, zFocus(pos));
            mmc.waitForSystem(); mmc.logMessage('line test'); pause(0.05);
            mmc.snapImage();
            tempImage = mmc.getImage();
            cameraImage = reshape(typecast(tempImage,'uint16'),1392,1040)';
            Temp1 = mmc.getPosition(mP.zDrive);
            Temp2 = fmeasure(cameraImage, AutoFunc, []);
            [yout,x]=imhist(cameraImage, double(2^12));
            set(hImageDisp.hPlot,'XData',x,'YData',yout);
            set(hImageDisp.hImDisplay,'CData',cameraImage);
        end
        
        zFocusPos = Temp1;
    else
    
    zPos = mmc.getPosition(mP.zDrive);
    zRecord = NaN(2*nSteps+21,3);
    upZ = zPos+zRange;
    lowZ = zPos-zRange; %range is +- range to move above and below
   
    zFocus = zeros(nSteps+3,3);
   
    zFocus(1:nSteps,1) = linspace(lowZ,upZ,nSteps)';
    zIncrement = abs(zFocus(2,1)-zFocus(1,1));
    recInd = 1;
    zPhase = 1;
    
    for step = 1:nSteps
        zFocusPos = zFocus(step,1);
        mmc.setPosition(mP.zDrive, zFocusPos);
        mmc.waitForSystem(); mmc.logMessage('line 689'); pause(0.05);
        mmc.snapImage();
        tempImage = mmc.getImage();
        cameraImage = reshape(typecast(tempImage,'uint16'),1392,1040)';
        zFocus(step,1) = mmc.getPosition(mP.zDrive);
        zFocus(step,2) = fmeasure(cameraImage, AutoFunc, []);
        [yout,x]=imhist(cameraImage, double(2^12));
        set(hImageDisp.hPlot,'XData',x,'YData',yout);
        set(hImageDisp.hImDisplay,'CData',cameraImage) 
        zRecord(recInd,1:2) = zFocus(step,1:2);
        zRecord(recInd,3) = zPhase;
        recInd = recInd+1;
        drawnow;
    end
    zPhase = 2;
    [~,maxInds] = max(zFocus(1:nSteps,2));
    ind = maxInds(1);

    if ind ==1 | ind==nSteps %Check to see if focus is at limits
        upZ  = zFocus(ind)+zRange;
        lowZ = zFocus(ind)-zRange;
        zFocus(1:nSteps,1) = linspace(lowZ,upZ,nSteps)';
        for step = 1:nSteps
            zFocusPos = zFocus(step,1);
            mmc.setPosition(mP.zDrive, zFocusPos);
            mmc.waitForSystem(); mmc.logMessage('line 714');pause(0.05);
            mmc.snapImage();
            tempImage = mmc.getImage();
            cameraImage = reshape(typecast(tempImage,'uint16'),1392,1040)';
            zFocus(step,1) = mmc.getPosition(mP.zDrive);
            zFocus(step,2) = fmeasure(cameraImage, AutoFunc, []);
            [yout,x]=imhist(cameraImage, double(2^12));

            set(hImageDisp.hPlot,'XData',x,'YData',yout);
            set(hImageDisp.hImDisplay,'CData',cameraImage);
           zRecord(recInd,1:2) = zFocus(step,1:2);
           zRecord(recInd,3) =  zPhase;
            recInd = recInd+1;
            drawnow;
        end
    end

    [~,maxInds] = max(zFocus(:,2));
    maxInd = maxInds(1);
    zFocus(nSteps+2,:) = zFocus(maxInd,:);

% Fine Tune Focus by sampleing between previous maximum points to find a
% local maximum any spot that is low due to flickering of lamp will be
% sorted out. 

numIter = 1;
cmdFocus = zeros(size(zFocus,1),20);
 while (zIncrement>targetIncrement) & (numIter <= 20) 
   % Calculate space between already measured top three maxima
  if maxInd == 1;
       zFocus(nSteps+1,1) =  zFocus(maxInd,1) - zIncrement;
       stepVec = nSteps+1;
       zPhase = 4.1;
       
  elseif maxInd  == nSteps ...| maxInd == nSteps+3;
       zFocus(nSteps+1,1) =  zFocus(nSteps-1,1);
       zFocus(nSteps+3,1) =  zFocus(maxInd,1) + zIncrement;
       stepVec = nSteps+3;
       zPhase = 4.2;
      
  else
       zIncrement = zIncrement/2;
       zFocus(nSteps+1,1) = zFocus(maxInd,1)-zIncrement;
       zFocus(nSteps+3,1) = zFocus(maxInd,1)+zIncrement;
       stepVec = [nSteps+1,nSteps+3];
       zPhase = zIncrement;
  end
    zFocusPos = zFocus(nSteps+1,1)-10;
    mmc.setPosition(mP.zDrive, zFocusPos); % This extra move takes out slack in the z-drive;
    mmc.waitForSystem();pause(0.05);
    cmdFocus(:,numIter) = zFocus(:,1)- zPos;      
        
  for step = stepVec
        
        zFocusPos = zFocus(step,1);
        mmc.setPosition(mP.zDrive, zFocusPos);
        mmc.waitForSystem();mmc.logMessage('line 770');pause(0.05);
        mmc.snapImage;
        tempImage = mmc.getImage;
        cameraImage = reshape(typecast(tempImage,'uint16'),1392,1040)';
        zFocus(step,1) = mmc.getPosition(mP.zDrive);
        zFocus(step,2) = fmeasure(cameraImage, 'LAPV', []);
        [yout,x]=imhist(cameraImage, double(2^12));
        set(hImageDisp.hPlot,'XData',x,'YData',yout);
        set(hImageDisp.hImDisplay,'CData',cameraImage);
        zRecord(recInd,1:2) = zFocus(step,1:2);
        zRecord(recInd,3) = zPhase;
        recInd = recInd+1;
        drawnow;
   end
    [~,maxInds] = max(zFocus(:,2));
%     [~,zOrd] = sort(zFocus(:,1));
%     [~,maxInds] = max(diff(zFocus(zOrd,2))./diff(zFocus(zOrd,1)));
    maxInd = maxInds(1);
    zFocus(nSteps+2,:) = zFocus(maxInd,:);
%     [~, zOrd] = sort(zRecord(:,1));
%     [~, maxDiffInd] = max(diff(zRecord(zOrd,2))./diff(zRecord(zOrd,1)));
%     zFocus(nSteps+2,:) = zRecord(maxDiffInd,1:2);
   % The previous algorithm assums that points 2 and 3 have focus positions
   % above and below the max fmeasure value.  If the max fmeasure values is
   % at the highest or lowest focus position, change that position.
   % Otherwise reduce the increment and look above and below the highest scored
   % focus positions.

   numIter = numIter+1;
   
 end
zFocusPos = zFocus(maxInd,1);
[~,zOrd] = sort(zFocus(:,1));
[~,maxInd] = max(diff(zFocus(zOrd,2))./diff(zFocus(zOrd,1)));

%Set position at max zFocus found.
mmc.setPosition(mP.zDrive,zFocusPos);
mmc.waitForSystem();
% mmc.logMessage('line 779 end autofocus');
    end
    
    if restartFlag; start(mP.hVideoTimer); end

% if ~isdeployed;  
assignin('base','mP',mP); 
assignin('base','zRecord',zRecord);
assignin('base','zFocus',zFocus);

updateGUIXYZ;
end

%% Run ODELAY experiment
function ODELAYStart(hTimer, events)
    global mmc mP hCont
%Check for a StageData File.  If one exists then use it. 

A = exist('ODELAY_StageData.mat');
if A==2
    load('ODELAY_StageData.mat'); 
    mP.getFileFlag = false; 
    fprintf('Loaded ODELAY Stage File');
else 
    mP.filePath = uigetdir(cd,'Set ODELAY Directory');
    cd(mP.filePath);
    for well = 1:mP.numWells
                if ~isdir(mP.wellID{mP.wellIdx(well)})
                    mkdir(mP.wellID{mP.wellIdx(well)});
                end
            end
    fileName = 'ODELAY_StageData';
    save(fullfile(mP.filePath, fileName),'mP');
end


set(hCont.hRunODELAY, 'BackgroundColor', [0 1 0], 'String','Stop');
% if mP.getFileFlag
%         mP.filePath = uigetdir('D:\','Set ODELAY Directory');
%         cd(mP.filePath);
%             for well = 1:mP.numWells
%                 if ~isdir(mP.wellID{mP.wellIdx(well)})
%                     mkdir(mP.wellID{mP.wellIdx(well)});
%                 end
%             end
% end
% fileName = 'ODELAY_StageData';
% save(fullfile(mP.filePath, fileName),'mP');
end

function ODELAYStop(hTimer, events)
global mP hCont

%     set(hCont.hRunODELAY, 'BackgroundColor', [0 1 0], 'String','ODELAY!!!');
%     fileName = ['ODELAY_StageData'];
%     save(fullfile(mP.filePath, fileName),'mP');
%     mP.getFileFlag = false;
% %     if mP.ErrorOccured == true && mP.restartCnt<5
% %         mP.ErrorOccured = false;
%         start(mP.hODELAYTimer);
%     end

end

function ODELAYError(hTimer, events)
global  mmc mP 

 if strcmp(get(mP.hVideoTimer, 'Running'),'Running'); stop(mP.hVideoTimer); end
 
%         if strcmp(get(mP.hODELAYTimer, 'Running'),'Running'); stop(mP.hODELAYTimer); end;
%            mP.restartCnt = mP.restartCnt+1;
%            mP.ErrorOccured = true;
%            mmc.unloadAllDevices;
%            mmc.reset;
%            clear mmc
%            pause(5)
%            ActivateMicroscope;
%            pause(5)
%            mP.getFileFlag = false;
end

function ODELAYScanPlate(hTimer,events)
  global mmc mP hCont hNav hImageDisp
  
%   load('ODELAY_StageData.mat'); 
  skipval = false;
  iterNum = mP.iterNum;  %mP.iterNum must be set to 1 initially in Activate Microscope
  if iterNum > mP.totIter
      stop(mP.hODELAYTimer);
  end
  fileName = ['MMLogFile_',num2str(iterNum),'.txt'];
  logFile = fullfile(['D:\'],fileName);
  mmc.setPrimaryLogFile(logFile)
  
  lastImaged = mP.lastImaged;
  mmc.setConfig('ImageMode','BrightField');
%   mmc.setProperty('TIDiaShutter', 'State', 1);
  startWell  = 1; 
  if mP.wellNumber ~= 1   % This may mean that the iterations stopped prematurely and need to be restarted. 
      startWell  = mP.wellNumber;      
  end
  
  for wellNumber = startWell:mP.numWells
      mP.wellNumber  = wellNumber;
      if ~get(hCont.hRunODELAY, 'Value'); stop(mP.hODELAYTimer); skipval = true; break;end
      
      % Main File for recording Data
      scanWell(mP.wellIdx(wellNumber), iterNum);
      
      % Also record stage and program state position
      fileName = 'ODELAY_Monitor';
      mP.lastImaged(mP.wellIdx(wellNumber), 1) = now;
      lastImaged = mP.lastImaged;
      save(fullfile(mP.filePath,fileName), 'lastImaged');
      fileName = 'ODELAY_StageData';
      save(fullfile(mP.filePath, fileName),'mP');
  end
  
  if ~skipval % At the end of the plate cycle do these 
      mP.twoPhaseFocus = false;  % After 1 round of twoPhaseFocus set it to false.
      mP.wellNumber  = 1; % Reset wellNumber to exit system.
      correctZPos(mP.iterNum); % Set new focus Positions based on last zFocus Poistions.
      mP.iterNum = mP.iterNum+1;
      fileName = 'ODELAY_StageData';
      save(fullfile(mP.filePath, fileName),'mP');
      pause(3)
      
  else
      fileName = 'ODELAY_StageData';
      save(fullfile(mP.filePath, fileName),'mP');
  end
  mmc.setProperty(mP.transShutter, 'State', 0);
%   if mod(mP.iterNum,24)==0
%      mmc.unloadAllDevices;
%      mmc.reset;
%      clear mmc
%      pause(1)
%      ActivateMicroscope;
%      pause(1)
%   end
  
 end

%% Scan/Tile Images
function scanWell(wellNum, iterNum)
global mmc mP hCont hNav hImageDisp
   
   
   if strcmp(get(mP.hVideoTimer, 'Running'),'on'); stop(mP.hVideoTimer); restartFlag = true; end;
    
%    currConfig = mmc.getCurrentConfig('ImageMode');
%    if ~strcmpi(currConfig, 'BrightField');
%         mmc.setConfig('ImageMode', 'BrightField');
%         expTime = str2double(get(hCont.hCubeFocus(7),'String'));
%         mmc.setExposure(expTime);
%    end
%     
    
    rawImage   = zeros(mP.sensorSize(1), mP.sensorSize(2), mP.numTiles, 'uint16');
    xyzTime    = zeros(mP.numTiles, 5);
        
    xPos = mP.stageXYPos(wellNum, 1);
    yPos = mP.stageXYPos(wellNum, 2);
    zPos = mP.stageZPos(wellNum,iterNum);
    fprintf(['New Position x=',num2str(xPos),' y=' num2str(yPos), ' z=',num2str(zPos),'\n']);
    mmc.setXYPosition(mP.xyDrive,xPos, yPos);pause(0.1);
    mmc.waitForDevice(mP.xyDrive); 
    mmc.setPosition(mP.zDrive,zPos);pause(0.1);
    mmc.waitForSystem();
    
    indx = false(size(mP.wellID));
    indx(wellNum)=1;
    
    set(hCont.hNavCircles(indx), 'FaceColor',[1 0 0]);
    set(hCont.hNavCircles(~indx), 'FaceColor','none');
    
  % Find Focus using Autofocus
    focusMethods =  get(hCont.hAutofocusSelect,'String');
    AutoFunc = focusMethods(get(hCont.hAutofocusSelect,'Value'),:);
    zRange = str2double(get(hCont.hAutofocusRange,'String'));
    nSteps = str2double(get(hCont.hAutofocusSteps,'String'));
    targetIncrement = str2double(get(hCont.hAutofocusTargetInc,'String'));
    
    focusMethods2 =  get(hCont.hAutofocusSelect2,'String');
    AutoFunc2 = focusMethods(get(hCont.hAutofocusSelect2,'Value'),:);
    zRange2 = str2double(get(hCont.hAutofocusRange2,'String'));
    nSteps2 = str2double(get(hCont.hAutofocusSteps2,'String'));
    targetIncrement2 = str2double(get(hCont.hAutofocusTargetInc2,'String'));
    
    if mP.twoPhaseFocus 
      
        % course then fine focus to ger a good initial focus 
        [zFocus, zRecord] = Autofocus(AutoFunc, zRange, nSteps, targetIncrement);
        newZPos = zFocus-2.0;
        mmc.setPosition(mP.zDrive, newZPos);pause(0.1);
        mmc.waitForSystem(); 
        [zFocus, zRecord] = Autofocus(AutoFunc2 , zRange2, nSteps2, targetIncrement2);
        
    else
        [zFocus, zRecord] = Autofocus(AutoFunc2, zRange2, nSteps2, targetIncrement2); 
    end
    
    if  ~mP.twoPhaseFocus && abs(zFocus - zPos)>7
        mmc.setPosition(mP.zDrive, zPos);
        mmc.waitForSystem(); pause(0.1);
        [zFocus, zRecord] = Autofocus(AutoFunc, zRange, nSteps, targetIncrement);
        newZPos = zFocus-2.0;
        mmc.setPosition(mP.zDrive, newZPos);
        mmc.waitForSystem(); pause(0.1);
        [zFocus, zRecord] = Autofocus(AutoFunc2, zRange2, nSteps2, targetIncrement2);
    end
        
        
    % Record Focus Position
    mP.zFocusPos(wellNum,iterNum) = zFocus;
    % Calculate Tile Postions 
    currentX = mmc.getXPosition(mP.xyDrive);
    currentY = mmc.getYPosition(mP.xyDrive);
    xTilePos = mP.tilePos+currentX;
    yTilePos = mP.tilePos+currentY;

     zStackPos = mP.stackVec + zFocus;
    
    turretPos = 1:mP.nTurretPos-1;
    recordPos = turretPos(mP.RecChan(1:end-1)==1);
    
    pos = 0;
    for sliceNum = 1:numel(mP.stackVec)
        mmc.setPosition(mP.zDrive, zStackPos(sliceNum));
        for tileNum = 1:mP.numTiles
            pos = pos+1;
            mmc.setXYPosition(mP.xyDrive,xTilePos(pos,1), yTilePos(pos,2));
            mmc.waitForDevice(mP.xyDrive);
            pause(0.1);
            mmc.snapImage();
            tempImage = mmc.getImage;
            rawImage(:,:,pos) = reshape(typecast(tempImage,'uint16'),1392,1040)';
            x = mmc.getXPosition(mP.xyDrive);
            y = mmc.getYPosition(mP.xyDrive);
            z = mmc.getPosition(mP.zDrive);
            xyzTime(pos,:) = [x, y, z, now, sliceNum];
            [yout,x]=imhist( rawImage(:,:,pos), double(2^12));
            set(hImageDisp.hPlot,'XData',x,'YData',yout);
            set(hImageDisp.hImDisplay,'CData',rawImage(:,:,pos));
            drawnow;
            if ~isempty(recordPos)
                mmc.setProperty(mP.transShutter, 'State', 0);
                for turretPos = recordPos
                    mmc.setConfig('ImageMode', mP.stateNames{turretPos});
                    expTime = str2double(get(hCont.hCubeExp(turretPos),'String'));
                    mmc.setExposure(expTime);
                    mmc.waitForDevice('TIFilterBlock1');
                    mmc.snapImage();
                    tempImage = mmc.getImage;
                    fluorImage.(mP.stateNames{turretPos}).rawImage(:,:,pos) = reshape(typecast(tempImage,'uint16'),1392,1040)';
                    [yout,x]=imhist( rawImage(:,:,pos), double(2^12));
                    set(hImageDisp.hPlot,'XData',x,'YData',yout);
                    set(hImageDisp.hImDisplay,'CData',fluorImage.(mP.stateNames{turretPos}).rawImage(:,:,pos));
                    drawnow;
                end
                mmc.setConfig('ImageMode', 'BrightField');
                expTime = str2double(get(hCont.hCubeExp(7),'String'));
                mmc.setExposure(expTime);
            end
            
        end
    end
        
  
    fileName = [mP.wellID{wellNum},'_',num2str(iterNum, '%i')];
    
    numTry = 0;
    success = false;
    while numTry<=100 & success == false
        try
            if isempty(recordPos)
               save(fullfile(mP.filePath, mP.wellID{wellNum}, fileName),...
                     'rawImage', 'xyzTime', 'zFocus', 'wellNum', 'iterNum', 'zRecord');
            else
                save(fullfile(mP.filePath, mP.wellID{wellNum}, fileName),...
                     'rawImage', 'fluorImage', 'xyzTime', 'zFocus', 'wellNum', 'iterNum', 'zRecord');
            end
            success = true;
        catch
            pause(1);
            numTry = numTry+1;
        end
    end

end

%% Correct Z positions 
function correctZPos(iterNum)
global mP mcc hCont hImageDisp
% mP.xgrid 
% mP.ygrid 
% mP.stageXYPos numWells x 2 of xy positions
% mP.stageZPos numWells by numTimePoints zPositions over time

listpos = get(hCont.hSlideSelect, 'Value');
%     switch listpos
%         case 1 % single chamber slide
%            for well = 1:mP.numWells;
%                r = mP.rowIdx(well);
%                c = mP.colIdx(well);
%                vq(r,c-1) = mP.zFocusPos(mP.wellIdx(well),iterNum);            
%            end
% 
%             vq_conv = zeros(size(vq)+2);
%             vq_conv(2:end-1,2:end-1) = vq;
% 
%             vq_conv(1,2:end-1)   = 2.*vq(1,1:end)-vq(2,1:end);
%             vq_conv(end,2:end-1) = 2.*vq(end,1:end)-vq(end-1,1:end);
%             vq_conv(2:end-1,1)   = 2.*vq(1:end,1)-vq(1:end,2);
%             vq_conv(2:end-1,end) = 2.*vq(1:end,end)-vq(1:end,end-1);
% 
%             vq_conv(1,1)     = mean([vq_conv(1,1),vq_conv(1,2),vq_conv(2,2),vq_conv(2,1)]);
%             vq_conv(1,end)   = mean([vq_conv(1,end),vq_conv(end-1,end-1),vq_conv(2,end),vq_conv(1,end-1)]);
%             vq_conv(end,end) = mean([vq_conv(end,end),vq_conv(end-1,end),vq_conv(end-1,end-1),vq_conv(end,end-1)]);
%             vq_conv(end,1)   = mean([vq_conv(end,1),vq_conv(end-1,1),vq_conv(end,2),vq_conv(end-1,2)]);
% 
%             vq_test = conv2(vq_conv,[1,1,1;1,1,1;1,1,1]./9,'valid');
% 
%             dDist     = vq-vq_test;
%             rmsdDist  = sum(dDist(:).^2)*1.0;
%             vq_prime  = vq;
%             vq_prime(dDist>rmsdDist) = vq_test(dDist>rmsdDist); 
% 
%             vq_conv = zeros(size(vq)+2);
%             vq_conv(2:end-1,2:end-1) = vq_prime;
% 
%             vq_conv(1,2:end-1)   = 2.*vq_prime(1,1:end)-vq_prime(2,1:end);
%             vq_conv(end,2:end-1) = 2.*vq_prime(end,1:end)-vq_prime(end-1,1:end);
%             vq_conv(2:end-1,1)   = 2.*vq_prime(1:end,1)-vq_prime(1:end,2);
%             vq_conv(2:end-1,end) = 2.*vq_prime(1:end,end)-vq_prime(1:end,end-1);
% 
%             vq_conv(1,1)     = mean([vq_conv(1,1),vq_conv(1,2),vq_conv(2,2),vq_conv(2,1)]);
%             vq_conv(1,end)   = mean([vq_conv(1,end),vq_conv(end-1,end-1),vq_conv(2,end),vq_conv(1,end-1)]);
%             vq_conv(end,end) = mean([vq_conv(end,end),vq_conv(end-1,end),vq_conv(end-1,end-1),vq_conv(end,end-1)]);
%             vq_conv(end,1)   = mean([vq_conv(end,1),vq_conv(end-1,1),vq_conv(end,2),vq_conv(end-1,2)]);
% 
%             vq_zeroed = conv2(vq_conv,ones(3)./3^2,'valid');
% 
%             %once the next Z Position is calculated, uptdate the next mP.stageZPos
%             %array for the next scan.  
%             timePoint = iterNum+1;
%             for  well = 1:mP.numWells
%                  r = mP.rowIdx(well);
%                  c = mP.colIdx(well);
%                  mP.stageZPos(mP.wellIdx(well),timePoint) = vq_zeroed(r,c-1);
%             end
% 
%         case 2 % 5 chamber slide
%             % this will set the next focuse point to the last focus point.
%             % however the Autofocus range should be large enough so that
            % focus does not run away.
            timePoint = iterNum+1;
            for well = 1:mP.numWells;
                   mP.stageZPos(mP.wellIdx(well),timePoint) = mP.zFocusPos(mP.wellIdx(well),iterNum);
            end
%         end
        
end

function  ODELAYLoadState
global mP
    temp = load('ODELAY_StageData.mat', 'mP');
    fieldNames = fieldnames(temp.mP);
    tempInds = ~strcmp(fieldNames,'hVideoTimer')&~strcmp(fieldNames,'hODELAYTimer');
    
    for n = 1:numel(fieldNames)
        if tempInds(n)
            mP.(fieldNames{n}) = temp.mP.(fieldNames{n});
        end
    end
    delete(temp);
end

%% Close Figure Function
function CloseMicroscopeControl(hFig, events)
global mmc mP hCont hImageDisp hNav
        if strcmp(get(mP.hVideoTimer , 'Running'),'Running'); stop(mP.hVideoTimer ); end
        if strcmp(get(mP.hODELAYTimer, 'Running'),'Running'); stop(mP.hODELAYTimer); end
       
        delete(hImageDisp.hImagePanel);
        delete(hCont.hControlPanel);
        delete(mP.hVideoTimer);
        delete(mP.hODELAYTimer);
 
        if exist('mmc','var') 
            mmc.unloadAllDevices;
            mmc.reset;
        end
        
        clear -global mmc mP hCont hImageDisp hNav

        
end

