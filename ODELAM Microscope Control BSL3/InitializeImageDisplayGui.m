function hImageDisp = InitializeImageDisplayGui;
global mP

% if isempty(mP)
% InitializeMicroscopeProperties;
% end

hImageDisp.hImagePanel = figure('Name','MATLAB Micro-Manager',...
                             'Units','Pixels',...
                             'Toolbar','none',...
                             'Colormap',mP.colorMap,...
                             'Position',mP.posDisp);
                         

hImageDisp.hDisplayAx = axes('Parent',hImageDisp.hImagePanel,...
                            'Units','Normalized',...
                            'Position',[0.15,0.15,0.825,0.825]);

hImageDisp.hImDisplay = image(mP.InitImage,...
                              'Parent',hImageDisp.hDisplayAx,...
                              'CDataMapping','scaled');
set(hImageDisp.hDisplayAx,'XTick',[],'YTick',[],'Visible','off');                        

hImageDisp.hHistAx = axes('Parent',hImageDisp.hImagePanel,...
                          'Units','Normalized',...
                          'Position',[0.05,0.15,0.08,0.825]);
                        
hImageDisp.hPlot = line(-3:0.1:3, normpdf([-3:0.1:3]),...
                        'Parent',hImageDisp.hHistAx);
                    
set(hImageDisp.hHistAx, 'View',[-90,90],...
                        'XTick',[],'YTick',[],...
                        'XLim',[1,double(4096)],...
                        'XAxisLocation','top');


end