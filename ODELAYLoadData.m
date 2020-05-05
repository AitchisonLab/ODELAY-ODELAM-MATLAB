function  varargout = ODELAYLoadData(DataFile)

%Use uigetfile to load data from datasets without using the comand line
CurrentDir = cd;
switch nargin % See if anything is entered if not get file with UI.
    case 1
        [filepath, filename, fileExt] = fileparts(DataFile);
        cd(filepath);
        
    otherwise
        [filename,filepath,~] = uigetfile('*Index_ODELAYData.mat');
         cd(filepath)
end
% if (filename==0) | (filepath==0); return;   end
       
m = matfile(filename);
varList = who(m);

for n = 1:size(varList,1)
    load(filename,varList{n});
end

if exist('ExperimentName','var')
    Experiment_Name = ExperimentName;
end

if isdir('ODELAY Well Data');
    cd('ODELAY Well Data');
else
    error('No save well data folder exists.  Please Generate one');
end
hWait = waitbar(0,'Loading ODELAY Data');
for well = 1:numwells
    ExData = load([Experiment_Name,'_Well_',num2str(well,'%0.2d'),'_ODELAYData']);
    ExFields = fieldnames(ExData);
    
    for ctr = 1:size(ExFields,1)
        switch ExFields{ctr}
            
            case 'WellDataTemp'
               
                
                WellDataFields = fieldnames(ExData.WellDataTemp);
                
                for fieldNum = 1:size(WellDataFields,1);
                    WellData(well,1).(WellDataFields{fieldNum}) = ...
                        ExData.WellDataTemp.(WellDataFields{fieldNum,:});
                    
                end
                assignin('caller','WellData',WellData);
                
            case 'WellOutTemp'
                WellOutFields = fieldnames(ExData.WellOutTemp);
                removeFields = {'prevBwImage';'prevLblImage';'prevFFT';'prevbw';'initDateNum';};
                rmVec = false(numel(WellOutFields),1);
                for n = 1:numel(removeFields)
                    rmVec = rmVec | strcmpi(WellOutFields, removeFields{n});
                end
                WellOutFields(rmVec) = [];
                
                 for fieldNum = 1:size(WellOutFields,1)
                        WellOut(well,1).(WellOutFields{fieldNum}) = ...
                            ExData.WellOutTemp.(WellOutFields{fieldNum,:});
                 end
                 assignin('caller','WellOut',WellOut);
                 
            case 'Tracks2Temp'
                Tracks2Fields = fieldnames(ExData.Tracks2Temp);
                 for fieldNum = 1:size(Tracks2Fields,1);
                        Tracks2(well,1).(Tracks2Fields{fieldNum}) = ...
                            ExData.Tracks2Temp.(Tracks2Fields{fieldNum,:});
                 end
                 assignin('caller','Tracks2',Tracks2);
                 
            case 'Tracks3Temp'
                Tracks3Fields = fieldnames(ExData.Tracks3Temp);
                 for fieldNum = 1:size(Tracks3Fields,1);
                        Tracks3(well,1).(Tracks3Fields{fieldNum}) = ...
                            ExData.Tracks3Temp.(Tracks3Fields{fieldNum,:});
                 end
                 assignin('caller','Tracks3',Tracks3);
        end
    end
    waitbar(well/numwells, hWait,'Loading ODELAY Data Set');
end
cd(CurrentDir);


clear('WellDataFields',...
      'WellOutFields',...
      'Tracks2Fields',...
      'ExFields',...
      'ExData',...
      'ctr',...
      'fieldNum',...
      'filepath',...
      'well');
% if nargout==1
%     varargout{1}.CurrentCD
% else
    if exist('CurrentCD','var');        
        assignin('caller','CurrentCD',CurrentCD); end
    if exist('Experiment_Name','var');  
        assignin('caller','Experiment_Name',Experiment_Name); end
    if exist('ImageDirectory','var');   
        assignin('caller','ImageDirectory',ImageDirectory); end
    if exist('ImageVars','var');        
        assignin('caller','ImageVars',ImageVars); end
    if exist('TimeDilute','var');       
        assignin('caller','TimeDilute',TimeDilute); end
    if exist('numwells','var');         
        assignin('caller','numwells',numwells); end
    if exist('saveWellFile','var');     
        assignin('caller','saveWellFile',saveWellFile); end
    if exist('tiffImageFiles','var');     
        assignin('caller','tiffImageFiles',tiffImageFiles); end
% end
delete(hWait);
