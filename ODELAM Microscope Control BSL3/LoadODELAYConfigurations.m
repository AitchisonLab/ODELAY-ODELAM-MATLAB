function LoadODELAYConfigurations
    global mmc mP 

    % SidePort 0 = eyepiece; 1 = 50% Left 50% eyes; 2 = 100% Left; 3 = 100%
    % right
    
    % IL-Turret 0 = 525; 1 = 565; 2 = 605; 3 = 655; 4 = 705; 5 = GFP/ATL/BF; 
    
    %% Brightfield Image Conditions
    mmc.defineConfig('ImageMode', 'BrightField', 'IL-Shutter','State','0');
    mmc.defineConfig('ImageMode', 'BrightField', 'TL-Shutter','State','1');
    mmc.defineConfig('ImageMode', 'BrightField', 'Core','AutoShutter','0');
    mmc.defineConfig('ImageMode', 'BrightField', 'Transmitted Light','State','1');
    mmc.defineConfig('ImageMode', 'BrightField', 'Transmitted Light','Level','60');
    mmc.defineConfig('ImageMode', 'BrightField', 'IL-Turret','State', '5');
    mmc.defineConfig('ImageMode', 'BrightField', 'SidePort', 'State','2');
    mmc.defineConfig('ImageMode', 'BrightField', 'TL-FieldDiaphragm','Position', '19'); 
   
     
       %% Darkfield Image Conditions
    mmc.defineConfig('ImageMode', 'DarkField', 'IL-Shutter','State','0');
    mmc.defineConfig('ImageMode', 'DarkField', 'TL-Shutter','State','1');
    mmc.defineConfig('ImageMode', 'DarkField', 'Transmitted Light','State','1');
    mmc.defineConfig('ImageMode', 'DarkField', 'Transmitted Light','Level','255');
    mmc.defineConfig('ImageMode', 'DarkField', 'Core','AutoShutter','0');
    mmc.defineConfig('ImageMode', 'DarkField', 'TL-FieldDiaphragm','Position', '46');
    mmc.defineConfig('ImageMode', 'DarkField', 'IL-Turret','State', '5');
    mmc.defineConfig('ImageMode', 'DarkField', 'SidePort', 'State','2');
   %????
   
    
       %% FluorescentImage 605 
    mmc.defineConfig('ImageMode', 'Fluor605', 'IL-Shutter','State','0');
    mmc.defineConfig('ImageMode', 'Fluor605', 'TL-Shutter','State','0');
    mmc.defineConfig('ImageMode', 'Fluor605', 'IL-Turret','State', '2');
    mmc.defineConfig('ImageMode', 'Fluor605', 'SidePort', 'State','2');
    mmc.defineConfig('ImageMode', 'Fluor605', 'Core','Shutter','IL-Shutter');
    mmc.defineConfig('ImageMode', 'Fluor605', 'Core','AutoShutter','1');
    
         %% FluorescentImage 655 
    mmc.defineConfig('ImageMode', 'Fluor655', 'IL-Shutter','State','0');
    mmc.defineConfig('ImageMode', 'Fluor655', 'TL-Shutter','State','0');
    mmc.defineConfig('ImageMode', 'Fluor655', 'Core','Shutter','IL-Shutter');
    mmc.defineConfig('ImageMode', 'Fluor655', 'Core','AutoShutter','1');
    mmc.defineConfig('ImageMode', 'Fluor655', 'IL-Turret','State', '3');
    mmc.defineConfig('ImageMode', 'Fluor655', 'SidePort', 'State','2');
    
   
         %% FluorescentImage 705 
    mmc.defineConfig('ImageMode', 'Fluor705', 'IL-Shutter','State','0');
    mmc.defineConfig('ImageMode', 'Fluor705', 'TL-Shutter','State','0');
    mmc.defineConfig('ImageMode', 'Fluor705', 'Core','Shutter','IL-Shutter');
    mmc.defineConfig('ImageMode', 'Fluor705', 'Core','AutoShutter','1');
    mmc.defineConfig('ImageMode', 'Fluor705', 'IL-Turret','State', '4');
    mmc.defineConfig('ImageMode', 'Fluor705', 'SidePort', 'State','2');
    
end