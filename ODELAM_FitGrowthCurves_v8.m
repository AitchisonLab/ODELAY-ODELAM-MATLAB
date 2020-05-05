function DataOut = ODELAM_FitGrowthCurves_v8(TracksIn, varargin)
%% Integrates the gompertsFitParam function into the image pipeline to fit 
%% growth curves of colonies linked over time to a parameterized version of 
%% the Gompertz function.
%==========================================================================
%% Author: Thurston Herricks 
% Institute for Systems Biology
% 401 Terry Ave N
% Seattle, WA 98109 USA
% Emails: 
% Thurston.Herricks@systemsbiology.org
%==========================================================================
% Last Modified: 2016/03/01
numargs = numel(varargin);

if numargs(1)~=0
    for n = 1:numel(TracksIn)
        wellLbl(n) = {TracksIn(n).ObjectInfo.WellID};
    end
    Tracks = TracksIn(strcmpi(wellLbl,varargin(1,2))).ObjectInfo;
else
    Tracks = TracksIn.ObjectInfo;
end

timePoints = Tracks.TimePoints;
imageList = 1:numel(timePoints);

objectArea = Tracks.ObjectArea;

%Check to see if the object area shrinks between the first and second
%timepoint if it does the focus may be bad and use the second time point to
%start fitting.

%Set up fitting with gompertsfit functions.       
%  FitDataGomp = NaN(size(objectArea,1),17);
 FitDataGompDT = NaN(size(objectArea,1),17);
%  FitDataGomp2 = NaN(size(objectArea,1),17);
if size(objectArea,2)>1; %Catch the problem if ObjectArea is filled
    dArea = diff(log2(objectArea(:,1:2)),1,2);
    if dArea<-1; tp1 = 2;else tp1 = 1;end
    else tp1 = 1;
end

diffAreaStd = abs(diff(log2(objectArea),1,2));
dbInds = diffAreaStd>30;
bgSteps = cumsum(dbInds,2)~=0;
objectArea(bgSteps)= NaN;

trackLength = sum(~isnan(objectArea),2);
trackVec = 1:size(objectArea,1);
trackInds = trackVec(trackLength>10);

%Set up fitting with gompertsfit functions.  
if ~isempty(trackInds);
    for n = trackInds
            inds = ~isnan(objectArea(n,:));
            tData    = timePoints(inds)';
            areaData = objectArea(n,inds);
            FitDataGompDT(n,1:end-1) = gompertzFitDT(tData, areaData, initTvTLag);
    end
end
%% FitData Column Key 
%Col     1   2   3   4  5     6     7     8     9       10        11        12       13
%Header 'a' 'b' 'T-lag' 'dT' 'fval' 'Tlag' 'Td' 'Tex' 'ATex' 'Aplateau' 'TdFlag' 'TexFlag' 'TVmax'
%     a:        Gompertz function value, a for y=a+b*exp(-exp(-ct+d))
%     b:        Gompertz function value, b for y=a+b*exp(-exp(-ct+d))
%     c:        T-lag
%     d:        dT or time between lag time and vmax 
%     fval:     Sum of squares (modeled versus experimental) from fmincom
%     Tlag:     Extracted lag time
%     Td:       Extracted doubling time
%     Texp:     Time at which modeled curve exits exponential phase
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
DataOut = Tracks;
flagIndex = ~isnan(FitDataGompDT(:,12));
DataOut(1).FitDataGompDT = FitDataGompDT;
DataOut(1).flagIndex     = flagIndex;


end

function prmsAll = gompertzFitDT(tData, iareas)
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
% October 2010; Last revision: 2020/04/30
%
    numtimepoints = sum(~isnan(iareas));

    idata = log2(iareas(~isnan(iareas)));
    tdata =  tData(~isnan(iareas));
    Nsteps=40;
    
    % Start geometric approxiamtion of initial values
    % Smooth data
    smthArea = smoothdata(idata,'movmean',5);
    
    % find line between first and last observed datapoints
    x = tdata([1,numtimepoints]);
    y = smthArea([1,numtimepoints]);
    m = diff(y)/diff(x);
    yVals = m*tdata + y(1)-m*x(1);
    
    diffVals = smthArea-yVals;
    cumVals  = cumsum(diffVals);
    % The min and max of the diff and cumulative values indicate the lag
    % time, t-exp and maximum velocity
    [~, lagInd] = min(diffVals);
    [~, texInd] = max(diffVals);
    [~ ,vmxInd] = min(cumVals);
    
   % Check values to ensure they are of the correct order.  
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
    % If the other operations fail use course grid optimization function 
    % findPrmsGompF to find a local minima 
    
        vecDT=linspace(1,100*60,Nsteps);
        bmin = 0;
        bmax = 10;
        vecTlag = linspace(1,tdata(end),Nsteps);
        vecB = linspace(bmin,bmax,Nsteps);
        [estB, estTlag,estDT,  fssq1] = findPrmsGompBDt(vecB,  vecTlag, vecDT, tdata, idata);
    
   end
    
   % assume the a parameter is the mean of the first 5 data points.
    meanArea = mean(idata(1:5));
    stdArea  = std(idata(1:5),0,2);
    Klag = log((3+sqrt(5))/2);
    prms0 = [meanArea estB  estTlag estDT];

    meanArea = mean(idata(1:5));
    stdArea  = std(idata(1:5),0,2);
    Klag = log((3+sqrt(5))/2);
    
    % Set limits for fmincon
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
    % Define teh gompertz function to be minimized
    gompFunc = @(prms)gompParBDt(prms, tdata, idata);
    options = optimoptions('fmincon','Algorithm','interior-point',...
                           'Display','off');
    Aeq = [0, 1, 0, 1]; B = [tdata(end)];
    
    % Perform minimization
    [prms, fssq, exitflag] = fmincon(gompFunc,prms0,Aeq,B,[],[],lb,ub,[],options);
    a    = prms(1);
    b    = prms(2);
    Tlag = prms(3);
    dT   = prms(4);
    % Extract perameters
    Klag = log((3+sqrt(5))/2);
    Kex  = log((3-sqrt(5))/2);
    c = Klag/dT;
    d = Tlag*c+Klag;
    
    TVmax = d/c;
    Tex = dT*2;
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
    % Place parameters in output.  
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

 

