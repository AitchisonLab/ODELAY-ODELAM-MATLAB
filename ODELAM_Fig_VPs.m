function ODELAM_Fig_VPs(ImageVars,Tracks2, Experiment_Name, varargin)
%% Generate Standard Plots of colony growth histograms for 
%% comparing colony growth data
%==========================================================================
%% Author: Thurston Herricks 
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% Emails: 
% Thurston.Herricks@systemsbiology.org
%==========================================================================
% Last Modified: 2016/03/01

colInd =[1:16];
lag = 6;
dbl = 7;
a = 1;
b = 2;
% Initialize Strain ID and Experiment Parameters

StrainID =ImageVars.StrainID;

strCnt =1;
numelData = zeros(numel(Tracks2),1);
wellInd = 1:numel(Tracks2);
cntr = 1;
if isfield(Tracks2(1).ObjectInfo, 'FitDataGompDT')
    dataName = 'FitDataGompDT';
else
    dataName = 'FitDataGomp';
end


for well = wellInd
    numelData(strCnt,:) = size(Tracks2(well).ObjectInfo.(dataName),1);
    LblData(strCnt,:) = [StrainID(well,2),...
                         StrainID(well,1),...
                         StrainID(well,4),...
                         [],...
                         {well},...
                         cntr,....
                         [StrainID{well,2},' ',...
                          StrainID{well,6},' ',...
                          StrainID{well,1}]];
    strCnt = strCnt+1;
end


pRange   = 0:100;
popData  = NaN(strCnt-1,numel(colInd),max(numelData(:,1)));
numTimePoints =NaN(numel(Tracks2), 5000);
percData = zeros(strCnt-1,numel(colInd),numel(pRange));
popSize = zeros(strCnt-1,1);
strCnt =1;

for well = 1:numel(Tracks2)
    for col = 1:numel(colInd)
        [tempGomp, numTPs, popSize(well,1)] = extractFitData(Tracks2, well, colInd(col));
        popData(strCnt,col,1:size(tempGomp,1)) = tempGomp;
        percData(strCnt,col,:) = prctile(popData(strCnt,col,:),pRange);
    end
    strCnt = strCnt+1;
end

numTimePoints =NaN(numel(Tracks2), max(popSize,[],1));
for well = 1:numel(Tracks2)
    numTimePoints(well,1:popSize(well,1)) = sum(~isnan(Tracks2(well).ObjectInfo.ObjectArea),2);
end

timeDiv = 60;

for idx = 1:nargin-3
switch varargin{idx}
    
    %% FitData Column Key 
%Col     1    2    3      4     5      6      7    8     9       10        11        12       13
%Header 'a'  'b' 'vmax' 'tlag' 'fval' 'Tlag' 'Td' 'Tex' 'ATex' 'Aplateau' 'TdFlag' 'TexFlag' 'TVmax'
%Header 'a'  'b' 'tlag' 'dT'   'fval' 'Tlag' 'Td' 'Tex' 'ATex' 'Aplateau' 'TdFlag' 'TexFlag' 'TVmax'
    case 'old'
        % a    =x(1);
        % b    =x(2);
        % vmax =x(3);
        % tlag =x(4);
        % yn=a+(exp(1)*vmax*dT/Klag).*exp(-exp((Klag/dT).*(dT+tlag-tdata)));

        DT   = squeeze(popData(:,8,:));
        vMax = squeeze(popData(:,3,:));
        B    = DT.*vMax.*exp(1)./log((3+sqrt(5))/2);
        Tex  = squeeze(popData(:,8,:));
        Tlag = squeeze(popData(:,6,:));
        Td   = squeeze(popData(:,7,:));
        
    case 'new'
        % a    =x(1);
        % b    =x(2);
        % tlag =x(3);
        % dT =x(4);
        DT   = squeeze(popData(:,4,:));
        B    = squeeze(popData(:,2,:));
        Tex  = DT.*2;
        Tlag = squeeze(popData(:,3,:));
        Td   = squeeze(popData(:,7,:));
        
    case 'Mabs'
          timeDiv = 60;
          plotRange = setPlotRange('Mabs');
        
    case 'Mtb'
        timeDiv = 60;
        plotRange = setPlotRange('Mtb');
        
    case 'yeast'
        timeDiv = 1;
        plotRange = setPlotRange('yeast');    
        
    case 'sort'
end
end
flagIndx = B>0 & numTimePoints>10;
popSize = sum(flagIndx,2);
timePoint = 60*48; %Time in Minutes
% sortInd = numel(Tracks2):-1:1;
[~,sortInd]  = sort(LblData(:,6));
sortInd = flipud(sortInd);
% [~,sortInd]  = sort(squeeze(percData(:,7,50)));
%Get Indicies of various DataSets


figure1 = figure('Color',[1 1 1],...
                 'Units', 'Inches',...
                 'PaperPositionMode','manual',...
                 'PaperSize',[8.5,10],...
                 'PaperPosition',[0.25,0.25,8,9.75],...
                 'Position', [0.5 0.5 6.5 9.75]);
left   = 1.25;
bottom = 0.4;
width  = 1;
height = 8.75;
spacing = 0.05;
annotation('textbox',[0.10,9.25,5.5,0.5]./[8,9.75,8,9.75],...
           'String',Experiment_Name,...
           'FontWeight','bold',...
           'FontSize',16,... ...'FitBoxToText',
            'LineStyle','none');
        
ViolinPlot(figure1,...
           B(sortInd,:).*flagIndx(sortInd,:),...
           LblData(sortInd,6),...
           'NumDbl',...
           plotRange,...
           '',...
           [left,bottom,width,height]);

ViolinPlot(figure1,...
           Td(sortInd,:)./timeDiv.*flagIndx(sortInd,:),...
           [],...LblData(sortInd,1),...
           'Dbl',...
           plotRange,...
           '',...
           [left+(width+spacing),bottom,width,height]);

ViolinPlot(figure1,...
           Tlag(sortInd,:)./timeDiv.*flagIndx(sortInd,:),...
           [],...LblData(sortInd,1),...
           'Lag',...
           plotRange,...
           '',...
           [left+2*(width+spacing),bottom,width,height]);

ViolinPlot(figure1,...
           Tex(sortInd,:)./timeDiv.*flagIndx(sortInd,:),...
           [],...LblData(sortInd,1),...
           'Tex',...
           plotRange,...
           '',...
           [left+3*(width+spacing),bottom,width,height]);
       
PopNumber(figure1,...
           popSize(sortInd),...
          'PopNum',...
           plotRange,...
           '',...
           [left+4*(width+spacing),bottom,width,height]);


end

function [gompOut, numTimePoints, numDtInd] = extractFitData(Tracks,trackInd, fitInd)

for well = 1:size(Tracks,1)
     wellPop(well) = size(Tracks(well).ObjectInfo.FitDataGompDT,1);
end

gompOut = NaN(max(wellPop(:)),numel(trackInd));

    for cnt = 1:numel(trackInd)
        well = trackInd(cnt);
        FitData = Tracks(well).ObjectInfo.FitDataGompDT;
        DT = FitData(:,3);
        vMax = FitData(:,2);
        flagInd = true(numel(FitData(:,1)),1);
        numDtInd(cnt,1) = sum(flagInd);
        gompOut(1:numDtInd(cnt,1),cnt)= Tracks(well).ObjectInfo.FitDataGompDT(flagInd,fitInd);
       
    end
try
numTimePoints = sum(~isnan(Tracks(trackInd).ObjectInfo.ObjectArea), 2);
catch
    test = 1;
end
for well = 1:size(Tracks,1)
     wellPop(well) = size(Tracks(well).ObjectInfo.FitDataGompDT,1);
end

lineOut = NaN(max(wellPop(:)),numel(trackInd));

%     for cnt = 1:numel(trackInd)
%         well = trackInd(cnt);
%         FitData = Tracks(well).ObjectInfo.FitDataLinear;
% %         flagInd = Tracks(well).ObjectInfo.flagIndex;
%         flagInd = (FitData(:,1)>0);
%         numInd = sum(flagInd);
%         if isfield(Tracks(well).ObjectInfo, 'FitDataLinear')
% %             lineOut(1:numInd,cnt)= Tracks(well).ObjectInfo.FitDataLinear(flagInd,fitInd);
%         end
%     end    
    
    
    
end

function ViolinPlot(figure1, StrainData, StrainLabels, paramLbl, plotRange, TitleBlock, plotPos)

    xRange    = plotRange.(paramLbl).xRange;
    xStep     = plotRange.(paramLbl).xStep;
    xLabel    = plotRange.(paramLbl).xLabel;
    titleFrag = plotRange.(paramLbl).titleFrag;

sizeData = size(StrainData);
if numel(sizeData>2);
    StrainData = squeeze(StrainData);
end

titleBlock = {[titleFrag]};

Label_Font = 6;
Title_Font = 10;
Number_Font = 6;
Legend_Font = 6;

wScale = 0.9;
xLabelPos = repmat(xRange(1)-diff(xRange)*0.05,1,numel(StrainLabels));
yLabelPos = [1:numel(StrainLabels)];

violinAxes = axes('Parent',figure1,...
                  'XLim',xRange,...
                  'YLim',[0 size(StrainData,1)+1],...
                  'Units','Inches',...
                  'Position',plotPos,...
                  'box','on',...
                  'LineWidth',2,...
                  'FontSize',Legend_Font,...
                  'TickLength',[0.0025,0.0025],...
                  'YTickLabel',[]);...xLabelPos);
% xlabel(xLabel);          
hold(violinAxes, 'on') 
nbins = linspace(xRange(1),xRange(2),30);
normVirts = zeros(size(StrainData,1),numel(nbins)+2);
xVirts = zeros(numel(nbins)+2,size(StrainData,1));
xVirts(2:end-1,:) = repmat(nbins',1,size(StrainData,1));
yVirts = zeros(size(xVirts,1), size(xVirts,2));

for cnt = 1:size(StrainData,1)
    virts = hist(StrainData(cnt,:),nbins);
    normVirts(cnt,2:end-1) = [virts./max(virts(2:end-1))].*wScale; 
    yVirts(:,cnt) = normVirts(cnt,:)'+cnt-0.5;
end
xVirts(1,:)   = xVirts(2,:);
xVirts(end,:) = xVirts(end-1,:);
yVirts(1,:)   = [1:size(yVirts,2)]-0.5;
yVirts(end,:) = [1:size(yVirts,2)]-0.5;
zVirts = ones(size(yVirts));

virtC = zeros(1,size(xVirts,2),3);
for n = 1:size(xVirts,2); virtC(1,n,1:3) = [0 0 1];end
hPatches = patch(xVirts, yVirts, zVirts,virtC,...
                 'EdgeColor','none',...
                 'Parent',violinAxes);  
 hYTick        = text(xLabelPos,...
                      yLabelPos,...
                      StrainLabels,...
                      'FontWeight', 'bold',...
                      'FontSize', Label_Font,...
                      'HorizontalAlignment','right',...
                      'rotation',0);
                  
XTickLabel = get(gca,'XTickLabel');
XTick = get(gca,'XTick');
yTickPos = repmat(-1,1,numel(XTick));
xTickPos = XTick;
set(gca,'XTickLabel',[]);

if xTickPos(1) == xRange(1)
    sTick = 2;
else
    sTick=1;
end
hxTick        = text(xTickPos(sTick:end),...
                     yTickPos(sTick:end),...
                     XTickLabel(sTick:end,:),...
                     'FontWeight', 'bold',...
                     'FontSize', Number_Font,...
                     'HorizontalAlignment','left',...
                     'rotation',-90);
                 
set(violinAxes,'box','on');
title( titleBlock,...
       'FontSize',Title_Font,...
       'FontWeight','bold');
end

function PopNumber(figure1, StrainData, paramLbl, plotRange, TitleBlock, plotPos)

    xRange    = plotRange.(paramLbl).xRange;
    xStep     = plotRange.(paramLbl).xStep;
    xLabel    = [];
    titleFrag = plotRange.(paramLbl).titleFrag;


titleBlock = {'Pop Num'};

Label_Font = 6;
Title_Font = 10;
Number_Font = 6;
Legend_Font = 6;

wScale = 0.9;
% xLabelPos = repmat(xRange(1)-diff(xRange)*0.05,1,numel(StrainLabels));
% yLabelPos = [1:numel(StrainLabels)];

popAxes = axes('Parent',figure1,...
              'YLim',[0 size(StrainData,1)+1],...
              'Units','Inches',...
              'Position',plotPos);...xLabelPos);
              
hbarPlot =  bar(StrainData,'Horizontal','on',...
                           'FaceColor','b',...
                           'BarWidth',0.5,...
                           'EdgeColor','b',...
                           'LineWidth',0.1);
                       
set(popAxes,'LineWidth',2,...
            'YTickLabel',[],...
            'XLim',[1,10000],...
            'XScale','log',...
            'XTick',[10,100,1000],...
            'FontWeight','bold',...
            'XTickLabelRotation',270,...
            'TickLength',[0.0025,0.0025],...
            'FontSize',Legend_Font);
            
                  
% XTickLabel = get(gca,'XTickLabel');
% XTick = get(gca,'XTick');
% yTickPos = repmat(-1,1,numel(XTick));
% xTickPos = XTick;
% set(gca,'XTickLabel',[]);
% 
% if xTickPos(1) == xRange(1)
%     sTick = 2;
% else
%     sTick=1;
% end
% hxTick        = text(xTickPos(sTick:end),...
%                      yTickPos(sTick:end),...
%                      XTickLabel(sTick:end,:),...
%                      'FontWeight', 'bold',...
%                      'FontSize', Number_Font,...
%                      'HorizontalAlignment','left',...
%                      'rotation',-90);
%                  
title( titleBlock,...
       'FontSize',Title_Font,...
       'FontWeight','bold');
end

function plotRange = setPlotRange(organism);

    switch organism

        case 'Mtb'
            plotRange.Dbl.xRange = [0, 100];
            plotRange.Dbl.xStep = 5;
            plotRange.Dbl.xLabel = 'Hours';
            plotRange.Dbl.titleFrag = 'Dbl Time Hr';

            plotRange.Lag.xRange = [0, 100];
            plotRange.Lag.xStep = 2;
            plotRange.Lag.xLabel = 'Hours';
            plotRange.Lag.titleFrag = 'Lag Time Hr';

            plotRange.Tex.xRange = [0, 100];
            plotRange.Tex.xStep = 2;
            plotRange.Tex.xLabel = 'Hours';
            plotRange.Tex.titleFrag = 'Tex Hr';

            plotRange.Area.xRange = [0, 30];
            plotRange.Area.xStep = 0.25;
            plotRange.Area.xLabel = 'log2 Pixels';
            plotRange.Area.titleFrag = 'log2 Area';

            plotRange.NumDbl.xRange = [0, 10];
            plotRange.NumDbl.xStep = 0.25;
            plotRange.NumDbl.xLabel = 'log2 Pixels';
            plotRange.NumDbl.titleFrag = 'Num Dbl Rel';

        case 'Mabs'
             plotRange.Dbl.xRange = [0, 10];
            plotRange.Dbl.xStep = 0.5;
            plotRange.Dbl.xLabel = 'Hours';
            plotRange.Dbl.titleFrag = 'Dbl Time Hr';

            plotRange.Lag.xRange = [0, 40];
            plotRange.Lag.xStep = 1;
            plotRange.Lag.xLabel = 'Hours';
            plotRange.Lag.titleFrag = 'Lag Time Hr';

            plotRange.Tex.xRange = [0, 40];
            plotRange.Tex.xStep = 1;
            plotRange.Tex.xLabel = 'Hours';
            plotRange.Tex.titleFrag = 'Tex Hr';

            plotRange.Area.xRange = [0, 30];
            plotRange.Area.xStep = 0.25;
            plotRange.Area.xLabel = 'log2 Pixels';
            plotRange.Area.titleFrag = 'log2 Area';

            plotRange.NumDbl.xRange = [0, 10];
            plotRange.NumDbl.xStep = 0.25;
            plotRange.NumDbl.xLabel = 'log2 Pixels';
            plotRange.NumDbl.titleFrag = 'Num Dbl Rel';

        case 'yeast'

            plotRange.Dbl.xRange = [25, 400];
            plotRange.Dbl.xStep = 4;
            plotRange.Dbl.xLabel = 'minutes';
            plotRange.Dbl.titleFrag = 'Dbl Time Min';

            plotRange.Lag.xRange = [0, 500];
            plotRange.Lag.xStep = 20;
            plotRange.Lag.xLabel = 'Minutes';
            plotRange.Lag.titleFrag = 'Lag Time Min';

            plotRange.Tex.xRange = [0, 1000];
            plotRange.Tex.xStep = 20;
            plotRange.Tex.xLabel = 'Minutes';
            plotRange.Tex.titleFrag = 'Tex Min';

            plotRange.Area.xRange = [0, 40];
            plotRange.Area.xStep = 0.5;
            plotRange.Area.xLabel = 'log2 Pixels';
            plotRange.Area.titleFrag = 'log2 Area';

            plotRange.NumDbl.xRange = [0, 10];
            plotRange.NumDbl.xStep = 0.25;
            plotRange.NumDbl.xLabel = 'log2 Pixels';
            plotRange.NumDbl.titleFrag = 'Num Dbl Rel';

        otherwise

            plotRange.Dbl.xRange = [0, 50];
            plotRange.Dbl.xStep = 2;
            plotRange.Dbl.xLabel = 'Hours';
            plotRange.Dbl.titleFrag = 'Doubling Time Min';

            plotRange.Lag.xRange = [0, 100];
            plotRange.Lag.xStep = 2;
            plotRange.Lag.xLabel = 'Hours';
            plotRange.Lag.titleFrag = 'Lag Time Min';

            plotRange.Tex.xRange = [0, 100];
            plotRange.Tex.xStep = 2;
            plotRange.Tex.xLabel = 'Hours';
            plotRange.Tex.titleFrag = 'Tex Min';

            plotRange.Area.xRange = [0, 30];
            plotRange.Area.xStep = 0.25;
            plotRange.Area.xLabel = 'log2 Pixels';
            plotRange.Area.titleFrag = 'log2 Area';

            plotRange.NumDbl.xRange = [0, 10];
            plotRange.NumDbl.xStep = 0.25;
            plotRange.NumDbl.xLabel = 'log2 Pixels';
            plotRange.NumDbl.titleFrag = 'log2 Area';
    end
    
    plotRange.PopNum.xRange=[1,10000];
    plotRange.PopNum.xStep = 10;
    plotRange.PopNum.xLabel = 'log10 Pop';
    plotRange.PopNum.titleFrag = 'Pop Num';
end
