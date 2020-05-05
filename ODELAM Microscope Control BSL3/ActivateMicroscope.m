
function ActivateMicroscope 
global mmc
% Start MicroManager and Import Core Java objects
import mmcorej.*;
mmc = CMMCore;

mmc.loadSystemConfiguration('C:\Program Files\Micro-Manager-1.4\ODELAY_BF.cfg');

% mmc.enableDebugLog(true);
% mmc.setPrimaryLogFile('D:\MATLAB\MMLogFile.txt',True);
% SidePort 0 = eyepiece; 1 = 50% Left 50% eyes; 2 = 100% Left; 3 = 100% right  
% IL-Turret 0 = 525; 1 = 565; 2 = 605; 3 = 655; 4 = 705; 5 = GFP/ATL/BF; 
stateNames = {'Cy5';'DAPI';'GFPHQ';'TxRed';'Blue';'AT-AQ';'BF';};
% Label,TIFilterBlock1,5,6-AT-AQ
% Label,TIFilterBlock1,4,5-Blue
% Label,TIFilterBlock1,3,4-TxRed
% Label,TIFilterBlock1,2,3-GFPHQ
% Label,TIFilterBlock1,1,2-DAPI
% Label,TIFilterBlock1,0,1-Cy5

% mP.xyDrive         = 'TIXYDrive';
% mP.zDrive          = 'TIZDrive';
% mP.transShutter    = 'TIDiaShutter';
% mP.fluorShutter    = 'IntensiLightShutter';
% mP.objectiveTurret = 'TINosePiece';
% mP.filterTurret    = 'TIFilterBlock1';
% mP.lightPath       = 'TILightPath';
% mP.AutoFocusFnc    = 'No Test';


%% Define hard coaded Imageing modes

 %% Brightfield Image Conditions
mmc.defineConfig('ImageMode', 'BrightField', 'TIFilterBlock1','State','2');
mmc.defineConfig('ImageMode', 'BrightField', 'TIDiaShutter','State','1');
mmc.defineConfig('ImageMode', 'BrightField', 'Core','AutoShutter','0');
% mmc.defineConfig('ImageMode', 'BrightField', 'TIEpiShutter','State', '0');
% mmc.defineConfig('ImageMode', 'BrightField', 'TIDiaLamp','State','1');
% mmc.defineConfig('ImageMode', 'BrightField', 'TIDiaLamp','Intensity','4');
% %  
for m = 0:5    
    mmc.defineConfig('ImageMode', stateNames{m+1}, 'TIFilterBlock1','State',num2str(m));
    mmc.defineConfig('ImageMode', stateNames{m+1}, 'TIDiaShutter','State','0');
%     mmc.defineConfig('ImageMode', stateNames{m+1}, 'TIEpiShutter','State', '1');
    mmc.defineConfig('ImageMode', stateNames{m+1}, 'Core','AutoShutter','0');
%     mmc.defineConfig('ImageMode', stateNames{m+1}, 'TIDiaLamp','State','0');
    
end

mmc.setConfig('ImageMode','BrightField');
mmc.setExposure(5);

end