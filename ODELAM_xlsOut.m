function ODELAM_xlsOut(Tracks, ImageVars, BaseFileName)
%==========================================================================
%% Author: Thurston Herricks 
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% Emails: 
% Thurston.Herricks@systemsbiology.org
%==========================================================================
% to run:  Load ODELAY Experiment Data: ODELAYLoadData
% then use this command
% ODELAY_xlsOut(Tracks2, ImageVars, Experiment_Name)
%Generate Time Delay (Td) verses 
clear StrainID TimePoints

numwells = size(Tracks,1);


%Col     1   2   3   4     5     6     7     8     9       10        11        12       13

numwells = size(Tracks,1);
% Well_Data = ODELAY_StrainID(Well_Data);


for well = 1:numwells

    Well_ID(well,1) = {Tracks(well).ObjectInfo.WellID};
    StrainID(well,1) = {[ImageVars.StrainID{well,3},' ',ImageVars.StrainID{well,4},' ',ImageVars.StrainID{well,1}]};
%     TimePoints(well,:) = Well_Data(well).ObjectInfo.TimePoints';
    objectArea = Tracks(well).ObjectInfo.ObjectArea;
    numObs = sum(~isnan(objectArea), 2); % sum how many time points are not NaN per tracked colony
    
    FitData =  Tracks(well).ObjectInfo.FitDataGompDT;
    flagIndx = ~isnan(FitData(:,1)) & numObs > 20;  % filter data acording to numObs (timepoints) observed
    medianFitData(well,:) = median(FitData(flagIndx,:));
    iqrFitData(well,:)  = abs(diff(prctile(FitData(flagIndx,:),[25,75])));
    
end

%Generate Summary File
% [prms fssq Tlag Td Tex ATex Aplateau TdFlag TexFlag TVmax Tplat exitflag fssq/numtimepoints]; 
Header = {'Strain ID','a','b','tlag','dT','fval','Tlag','Td','Tex','ATex','Aplateau','TdFlag','TexFlag','TVmax','TPlat','exitflag','fssq','N/A'};

medianFitsize = size(medianFitData);
xlsSheet = cell(medianFitsize(1)+1,medianFitsize(2)+1);

xlsSheet(1,:) = Header(1,:);
xlsSheet(2:end,1) = StrainID(:,1);
cellData = num2cell(medianFitData);

for i = 2:size(xlsSheet,1);
   xlsSheet(i,2:end) = cellData(i-1,:);
end

FileName = [BaseFileName,' Summary'];
xlswrite(FileName, xlsSheet, 'Filtered Means Summary');
clear xlsSheet

% Header = {'Strain ID','a','b','tlag','texp','fval','Tlag','Td','Tex','ATex','Aplateau','TdFlag','TexFlag','TVmax'};

medianFitSize = size(medianFitData);
xlsSheet = cell(medianFitSize(1)+1,medianFitSize(2)+1);

xlsSheet(1,:) = Header(1,:);
xlsSheet(2:end,1) = StrainID(:,1);
cellData = num2cell(iqrFitData);
for i = 2:size(xlsSheet,1);
   xlsSheet(i,2:end) = cellData(i-1,:);
end

xlswrite(FileName, xlsSheet, 'Std Summary');

clear xlsSheet

%Generate Individual Strain Data

for well = 1:numwells
    FitData = Tracks(well).ObjectInfo.FitDataGompDT;
%     StrainID = (well).StrainID{:};
    Well_ID = Tracks(well).ObjectInfo.WellID;
    Header = {'Strain ID','a','b','tlag','dT','fval','Tlag','Td','Tex','ATex','Aplateau','TdFlag','TexFlag','TVmax','TPlat','exitflag','fssq','N/A'};
    flagIndx = ~isnan(FitData(:,1));
    
    Fitsize = size(FitData(flagIndx,:));
    xlsSheet = cell(Fitsize(1)+1, Fitsize(2)+1);

    xlsSheet(1,:) = Header(1,:);
    xlsSheet(2:end,1) = num2cell(1:Fitsize(1));
    cellData = num2cell(FitData(flagIndx,:));

    for i = 2:size(xlsSheet,1);
        xlsSheet(i,2:end) = cellData(i-1,:);
    end

    FileName = [BaseFileName,'_',StrainID{well,1}];
    ind = FileName == '/';
    FileName(ind) = ' ';
    ind = FileName == '.';
    FileName(ind) = '_';
    xlswrite(FileName, xlsSheet, 'Fit Data');
    clear xlsSheet

    AreaData = Tracks(well).ObjectInfo.ObjectArea(flagIndx,:);
    TimePoints = Tracks(well).ObjectInfo.TimePoints;
    Header = {'Colony ID','Time Points (minutes)', 'Area given in number of filled pixels'};
    Areasize = size(AreaData);
    xlsSheet = cell(Areasize(1)+2,Areasize(2)+1);

    xlsSheet(1,1:2) = Header(1,1:2);
    xlsSheet(2,2:end) = num2cell(TimePoints');
    xlsSheet(3:end,1) = num2cell(1:Areasize(1));
    cellData = num2cell(AreaData);

    for i = 3:size(xlsSheet,1);
        xlsSheet(i,2:end) = cellData(i-2,:);
    end

%     xlswrite(FileName, xlsSheet, 'Area Data');
    clear xlsSheet

end
