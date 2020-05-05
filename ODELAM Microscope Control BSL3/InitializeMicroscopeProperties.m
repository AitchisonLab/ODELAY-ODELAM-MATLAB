function InitializeMicroscopeProperties
global mP
%% To do:
% set up file locator
% read data from uM config file and config file for initializing ODELAY or
% other experiments

% Device,TIScope,NikonTI,TIScope
% Device,TIDiaShutter,NikonTI,TIDiaShutter
% Device,TINosePiece,NikonTI,TINosePiece
% Device,TICondenserCassette,NikonTI,TICondenserCassette
% Device,TIFilterBlock1,NikonTI,TIFilterBlock1
% Device,TILightPath,NikonTI,TILightPath
% Device,TIZDrive,NikonTI,TIZDrive
% Device,TIXYDrive,NikonTI,TIXYDrive
% Device,Camera-1,PVCAM,Camera-1
mP.xyDrive         = 'TIXYDrive';
mP.zDrive          = 'TIZDrive';
mP.transShutter    = 'TIDiaShutter';
mP.fluorShutter    = 'TIEpiShutter';
mP.objectiveTurret = 'TINosePiece';
mP.filterTurret    = 'TIFilterBlock1';
mP.lightPath       = 'TILightPath';
mP.AutoFocusFnc    = 'No Test';


%% AutoFocus Props
mP.zRange = 70;
mP.numSteps = 11;
mP.targetIncrement = 0.3;

mP.zRange2 = 10;
mP.numSteps2 = 5;
mP.targetIncrement2 = 0.3;
mP.twoPhaseFocus = true;

%% ODELAY Props
mP.rowName = {'E';'F';'G';'H';'I';'J';'K';'L'};
mP.colName = {'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';};

mP.xgrid = -[0:4500:4500*13];
mP.ygrid = 0:4500:4500*7;

mP.numTimePoints =  432;
mP.numWells = numel(mP.rowName)*numel(mP.colName);
mP.stackVec = 0;
mP.stageXYPos = zeros(mP.numWells, 2);
mP.stageZPos  = zeros(mP.numWells, mP.numTimePoints+1); %Calculated stage positions before focus
mP.zFocusPos  = zeros(mP.numWells, mP.numTimePoints+1); %Stage Positions recorded after focus
mP.wellID = cell(mP.numWells,1);

mP.lastImaged = zeros(mP.numWells,1);
mP.XYZOrigin(1,:) = [29500,-15500,8500.0]; 
 
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

mP.iterNum = 1;
mP.wellNumber = 1;

mP.totIter = mP.numTimePoints;
mP.iterPeriod = 1800;
mP.getFileFlag = true;
mP.ErrorOccured = false;
mP.restartCnt = 0;

mP.numTiles   = 9;
mP.mag = 20;
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
%% Microscope Conditions              
%   mP.XYZOrigin = zeros(totIter, 3);
      
  
%% Figure Properties ***Move to  mP = InitializeMicroscopeProperties; ***
        mP.screenSize = get(0,'ScreenSize');
        mP.posControl = [1080, 400, 800, 550];
        mP.posDisp    = [50,50,800,800];
        mP.ColorList  = [525,565,605,655,705,525, 425];
        mP.ColorNames = {'Cy5';'DAPI';'GFPHQ';'TxRed';'Blue';'AT-AQ';'BF';};
        mP.stateNames = {'Cy5';'DAPI';'GFPHQ';'TxRed';'Blue';'AT-AQ';'BF';};
        
        mP.objectivLens = [20,20,20];
        mP.objTurretPos = [1,2,3];
        mP.LenseMenu     = ['20x|20x|20x'];
        
        mP.ExpTime    = [500;500;500;500;500;500;5];
        mP.nTurretPos = 7;
        mP.numImages  = 9;
        mP.CubeEngaged = [0,0,0,0,0,0,1];
        mP.RecChan     = [0,0,0,0,0,0,1];
      
    % Need to check Screen Size
       numbuttons = numel(mP.stateNames); 
       left = 0.125;
       spacing = 0.045;
       width = (1-left-spacing*(numbuttons+1))/numbuttons;
       height = 0.5;
       bottom = 0.135;
for n = 1:mP.nTurretPos
%     [left bottom width height]
     mP.CubePos(n,1:4)   = [left+(width+spacing)*(n-1), bottom, width, height];
     mP.ExpPos(n,1:4)    = [left+(width+spacing)*(n-1), bottom+0.55, width, 0.1];
     mP.FocusPos(n,1:4)  = [left+(width+spacing)*(n-1), bottom+0.7125, width, 0.1];
     mP.RecPos(n,1:4)    = [left+(width+spacing)*(n-1), bottom-0.125, width, 0.1];
     mP.CubeColor(n,:)   = Wavelength_to_RGB(mP.ColorList(n))./255;
     
end
mP.RecColor = zeros(7,3);
mP.RecColor(mP.RecChan==1,1:3)   = [0 1 0];
mP.RecColor(mP.RecChan~=1,1:3)   = repmat([1 0 0],sum(mP.RecChan~=1),1);

mP.InitImage = zeros(1080, 1392, 'uint16');
mP.colorMap = repmat(linspace(0,1,256)',1,3);
l = 0.15;
s = l;
w = (1-4*s)/3;

mP.EditXYZPos  = [0.40, 0.85, 0.125, 0.1;
                  0.55, 0.85, 0.125, 0.1;
                  0.70, 0.85, 0.125, 0.1;
                  0.85, 0.85, 0.125, 0.1;
                  0.40, 0.65, 0.125, 0.1;
                  0.55, 0.65, 0.125, 0.1;
                  0.70, 0.65, 0.125, 0.1];

mP.sethOriginPos = [0.85, 0.65, 0.125, 0.1];
mP.sethGoOrigin  = [0.05, 0.40, 0.2,   0.1];
mP.sethODELAYDir = [0.05, 0.25, 0.2,   0.1];
mP.sethRunODELAY = [0.05, 0.10, 0.2,   0.1];

mP.AutofocusMenu = ['ACMO|BREN|CONT|CURV|DCTE|DCTR|GDER|GLVA|GLLV|GLVN|GRAE|GRAT|GRAS|HELM|HISE|HISR|LAPE|LAPM|LAPV|LAPD|SFIL|SFRQ|TENG|TENV|VOLA'];

mP.posConfigMenu = [1080,400,400,500];
mP.posSlideSel   = [10,10,100,50];        

end