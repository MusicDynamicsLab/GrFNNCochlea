function TCfitGUI
% function TCfitGUI - GUI for fitting model tuning curves to empirical
% results. Change parameter values by moving the sliders or by entering
% numbers in the text boxes.

figure(1);
set(gcf, 'Visible', 'off', 'Toolbar', 'figure', 'Color', [.8 .8 .8], ...
    'Position', [5, 5, 900 800], 'Name', 'Tuning Curves');
initialRun;
movegui(gcf,'center');
set(gcf, 'Visible', 'on');



%% Initial run function for when GUI first is created
function initialRun(~, ~)
clf
axis off;

ax1 = axes('Units', 'normalized', 'Position', [.06 .65 .5 .3]);
ax2 = axes('Units', 'normalized', 'Position', [.63  .36 .35 .2]);
ax3 = axes('Units', 'normalized', 'Position', [.63  .045 .35 .2]);
ax4 = axes('Units', 'normalized', 'Position', [.06  .36 .5 .2]);
ax5 = axes('Units', 'normalized', 'Position', [.06  .045 .5 .2]);
ax = [ax1 ax2 ax3 ax4 ax5];

defaultparamsBMOC = csvread('defaultparamsBMOC.csv');

%% Initial parameter values
tcnum       = 5;
alphabm     = defaultparamsBMOC(10,tcnum);
beta1bm     = defaultparamsBMOC(11,tcnum);
delta1bm    = defaultparamsBMOC(12,tcnum);
alphaoc     = defaultparamsBMOC(13,tcnum);
beta1oc     = defaultparamsBMOC(14,tcnum);
delta1oc    = defaultparamsBMOC(15,tcnum);
c21         = defaultparamsBMOC(16,tcnum);
c12         = defaultparamsBMOC(17,tcnum);
thresh      = defaultparamsBMOC(18,tcnum);

clear defaultparamsBMOC
TCnumLim   = [1 11];


%% Make text boxes and labels for parameter input

% Panel for parameter input
hp1 = uipanel('Title','Parameters',...
    'FontSize',11,'FontName','Helvetica','FontUnits','normalized',...
    'BackgroundColor',[.8 .8 .8],'Position',[.71 .66 .28 .3]);
axes('parent', hp1,'position',[0 0 1 1], 'visible', 'off');

% Label which layer the textboxes correspond to
text(.55,.96,'BM','HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
text(.85,.96,'OC','HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')

% Label and textbox for alpha_BM
text(.45,.89,'alpha:','HorizontalAlignment','right','VerticalAlignment','middle',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
handles.alphabm = uicontrol('Parent',hp1,'Style','edit', 'String', alphabm,...
    'FontSize',12,'FontUnits','normalized','Units','normalized','Position', [.45 .83 .2 .13]);
% Label and textbox for beta_BM
text(.45,.74,'beta:','HorizontalAlignment','right','VerticalAlignment','middle',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
handles.betabm = uicontrol('Parent',hp1,'Style','edit', 'String', beta1bm,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.45 .68 .2 .13],'enable','off');
% Label and textbox for delta_BM
text(.45,.59,'delta:','HorizontalAlignment','right','VerticalAlignment','middle',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
handles.deltabm = uicontrol('Parent',hp1,'Style','edit', 'String', delta1bm,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.45 .53 .2 .13],'enable','off');

% Label and textbox for alpha_OC
handles.alphaoc = uicontrol('Parent',hp1,'Style','edit', 'String', alphaoc,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.76 .83 .2 .13]);
% Label and textbox for beta_OC
handles.betaoc = uicontrol('Parent',hp1,'Style','edit', 'String', beta1oc,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.76 .68 .2 .13]);
% Label and textbox for delta_OC
handles.deltaoc = uicontrol('Parent',hp1,'Style','edit', 'String', delta1oc,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.76 .53 .2 .13]);

% Label and textbox for c12
text(.45,.45,'c12:','HorizontalAlignment','right','VerticalAlignment','middle',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
handles.c12 = uicontrol('Parent',hp1,'Style','edit', 'String', c12,...    
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.45 .38 .2 .13]);
% Label and texbox for c21
text(.76,.45,'c21:','HorizontalAlignment','right','VerticalAlignment','middle',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
handles.c21 = uicontrol('Parent',hp1,'Style','edit', 'String', c21,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.76 .38 .2 .13]);

% Label and textbox for threshold
text(.05,.83,'Threshold:','HorizontalAlignment','left',...
    'VerticalAlignment','middle','Units','normalized','FontSize',12,'FontUnits','normalized')
handles.thresh = uicontrol('Parent', hp1, 'Style','edit', 'String', thresh,...
    'FontSize',12,'Units', 'normalized', 'HorizontalAlignment', 'center',...    
    'Position', [.05 .67 .2 .13],'FontUnits','normalized');

% Pushbutton controls for saving as default parameters and choosing c12
% value
handles.buttonSaveAsDefault = uicontrol('Parent', hp1, 'Style', 'pushbutton', 'String', {'Save as Default'},...
    'TooltipString', ['Save the current parameters' char(10)...
    'as the default parameters in' char(10) 'defaultparamsBMOC.csv'],...
    'FontUnits','normalized','Units', 'normalized',  'Position', [.5 .2 .4 .18]);
handles.buttonChoosec12 = uicontrol('Parent', hp1, 'Style', 'pushbutton', 'String', 'Choose c12',...
    'FontUnits','normalized','Units', 'normalized',...    
    'TooltipString', ['Set value of c12 to make spontaneous' char(10) 'amplitude equal to half the threshold.'],...
    'Position', [0 .37 .35 .15]);

% Tuning curve slider, textbox, and label
% axes('position',[.8 .62 .05 .02],'Visible','off');
text(.33,.13,'Tuning Curve:','HorizontalAlignment','right','VerticalAlignment','baseline',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
handles.textTC = uicontrol('Style','edit', 'String', tcnum,...
    'FontSize',12,'FontUnits','normalized','Units','normalized','Position', [.805 .69 .046 .03]);
handles.sliderTC = uicontrol('Style', 'slider', 'Min', 1, 'Max', 11, 'Value', tcnum, 'SliderStep', [.1 .4], ...
    'Min', TCnumLim(1), 'Max', TCnumLim(2), 'Units', 'normalized', 'Position', [.72 .66 .15 .025]);


% Create panel for toggle controls
hp2 = uipanel('Title','Toggles',...
    'FontSize',11,'FontName','Helvetica','FontUnits','normalized',...
    'BackgroundColor',[.8 .8 .8],'Position',[.57 .66 .13 .3]);
axes('Units','normalized','Position',get(hp2,'Position'))
axis off

% Toggle controls for middle ear filter, frequency scaling, setting organ
% of corti to be the treshold, and using a single oscillator network
handles.toggleMEF = uicontrol('Parent',hp2,'Style','togglebutton', 'String','MEF',...
    'TooltipString', 'Turn ON/OFF Middle Ear Filter',...
    'Value', 1,'FontUnits','normalized','Units', 'normalized', 'Position', [.03 0.8 0.95 0.15]);
handles.toggleFS = uicontrol('Parent',hp2,'Style','togglebutton', 'String','Freq Scaling',...
    'TooltipString', 'Turn ON/OFF Frequency Scaling',...
    'Value', 0,'FontUnits','normalized','Units', 'normalized', 'Position', [.03 .55 .95 .15]);
handles.toggleThresh = uicontrol('Parent',hp2,'Style','togglebutton', 'String','Thresh: OC',...
    'TooltipString', ['ON = use OC amplitude as threshold ' char(10) 'OFF = use BM amplitude as threshold'],...
    'Value', 1,'FontUnits','normalized','Units', 'normalized', 'Position', [.03 .3 .95 .15]);
handles.toggleSingle = uicontrol('Parent',hp2,'Style','togglebutton', 'String','Two Layer',...
    'TooltipString', ['ON = two layer model' char(10) 'OFF = single layer model'],...
    'Value', 1,'FontUnits','normalized','Units', 'normalized', 'Position', [.03 .05 .95 .15]);


%% Set callback functions
handlesTC           = [handles.sliderTC handles.textTC];
handles.textGroup   = [handles.alphabm handles.betabm handles.deltabm ...
                       handles.alphaoc handles.betaoc handles.deltaoc ...
                       handles.c21 handles.c12 handles.thresh];
handles.toggleGroup = [handles.toggleMEF handles.toggleFS handles.toggleThresh handles.toggleSingle];

set(handlesTC,                  'Callback', {@Callback_TC, handles, ax});
set(handles.textGroup,          'Callback', {@text_Callback,   handles, ax});
set(handles.buttonSaveAsDefault,'Callback', {@saveAsDef_Callback, handles});
set(handles.buttonChoosec12,    'Callback', {@choosec12_Callback, handles, ax});
set(handles.toggleGroup,        'Callback', {@toggle_Callback,handles,ax});

hold off
plotTC(tcnum, alphabm, beta1bm, delta1bm, alphaoc, beta1oc, delta1oc, c21, c12, thresh, {1, 0, 1, 1}, ax);

%% Toggle pushbutton callback function
function toggle_Callback(~, ~, handles, ax)

defaultparamsBMOC = csvread('defaultparamsBMOC.csv');
tcnum = get(handles.sliderTC,'Value');
freqScaling    = get(handles.toggleFS, 'Value');
threshOC       = get(handles.toggleThresh, 'Value');
singleCritical = get(handles.toggleSingle, 'Value');

if ~singleCritical
    if freqScaling
        params = defaultparamsBMOC(37:45,tcnum);
    else
        params = defaultparamsBMOC(46:54,tcnum);
    end
else
    
    if threshOC
        if freqScaling
            params = defaultparamsBMOC(1:9, tcnum);
        else
            params = defaultparamsBMOC(10:18,tcnum);
        end
    else
        if freqScaling
            params = defaultparamsBMOC(19:27, tcnum);
        else
            params = defaultparamsBMOC(28:36,tcnum);
        end
    end
    
end

set(handles.alphabm, 'String', num2str(params(1),3));
set(handles.betabm,  'String', num2str(params(2),3));
set(handles.deltabm, 'String', num2str(params(3),3));
set(handles.alphaoc, 'String', num2str(params(4),3));
set(handles.betaoc,  'String', num2str(params(5),3));
set(handles.deltaoc, 'String', num2str(params(6),3));
set(handles.c21,     'String', num2str(params(7),3));
set(handles.c12,     'String', num2str(params(8),3));
set(handles.thresh,  'String', num2str(params(9),3));

clear defaultparamsBMOC

parameterValues = str2double(get(handles.textGroup, 'String'));
toggleValues    = get(handles.toggleGroup, 'Value');
temp = num2cell(parameterValues);
plotTC(tcnum, temp{:}, toggleValues, ax)


%% Save as Default Parameters pushbutton
function saveAsDef_Callback(~, ~, handles)

pathToParams = which('defaultparamsBMOC.csv');
defaultparamsBMOC = csvread('defaultparamsBMOC.csv');

freqScaling    = get(handles.toggleFS, 'Value');
threshOC       = get(handles.toggleThresh, 'Value');
singleCritical = get(handles.toggleSingle, 'Value');

params = str2double(get(handles.textGroup, 'String'));

tcnum  = str2double(get(handles.textTC, 'String'));

if ~singleCritical
    if freqScaling
        defaultparamsBMOC(37:45,tcnum) = params;
    else
        defaultparamsBMOC(46:54,tcnum) = params;
    end
else
    
    if threshOC
        if freqScaling
            defaultparamsBMOC(1:9,tcnum)  = params;
        else
            defaultparamsBMOC(10:18,tcnum) = params;
        end
    else
        if freqScaling
            defaultparamsBMOC(19:27,tcnum)  = params;
        else
            defaultparamsBMOC(28:36,tcnum) = params;
        end
    end
    
end

csvwrite(pathToParams,defaultparamsBMOC);
clear defaultparamsBMOC


%% Set c12=abm/c21*(aoc+boc*(thresh*factor)^2)
function choosec12_Callback(~, ~, handles, ax)

singleCritical = get(handles.toggleSingle, 'Value');

if singleCritical
    
    params = str2double(get(handles.textGroup, 'String'));
    
    abm    = params(1);
    aoc    = params(4);    
    boc    = params(5);
    c21    = params(7);
    thresh = params(9);
    factor = 0.5;
    
    c12=abm/c21*(aoc+boc*(thresh*factor)^2);
    
    textValue = num2str(c12);
    
    set(handles.c12,'String',textValue);
    
    toggleValues    = get(handles.toggleGroup, 'Value');
    
    parameterValues = str2double(get(handles.textGroup, 'String'));
    tcnum = get(handles.sliderTC,'Value');
    temp = num2cell(parameterValues);
    plotTC(tcnum, temp{:}, toggleValues, ax)
    
end


%% Tuning curve callback function
function Callback_TC(source, ~, handles, ax)

if strcmp(get(source,'Style'),'edit')
    tcnum = round(str2double(get(source,'String')));    
    set(handles.textTC, 'String', num2str(tcnum));
    set(handles.sliderTC, 'Value', tcnum)
else
    tcnum = round(get(source,'Value'));
    set(handles.sliderTC, 'Value', tcnum);
    set(handles.textTC,  'String', num2str(tcnum));
end

defaultparamsBMOC = csvread('defaultparamsBMOC.csv');
toggleValues      = get(handles.toggleGroup, 'Value');

if ~toggleValues{4}
    if toggleValues{2}
        params = defaultparamsBMOC(37:45,tcnum);
    else
        params = defaultparamsBMOC(46:54,tcnum);
    end
else
    
    if toggleValues{3}
        if toggleValues{2}
            params = defaultparamsBMOC(1:9, tcnum);
        else
            params = defaultparamsBMOC(10:18,tcnum);
        end
    else
        if toggleValues{2}
            params = defaultparamsBMOC(19:27, tcnum);
        else
            params = defaultparamsBMOC(28:36,tcnum);
        end
    end
    
end

set(handles.alphabm, 'String', num2str(params(1),3));
set(handles.betabm,  'String', num2str(params(2),3));
set(handles.deltabm, 'String', num2str(params(3),3));
set(handles.alphaoc, 'String', num2str(params(4),3));
set(handles.betaoc,  'String', num2str(params(5),3));
set(handles.deltaoc, 'String', num2str(params(6),3));
set(handles.c21,     'String', num2str(params(7),3));
set(handles.c12,     'String', num2str(params(8),3));
set(handles.thresh,  'String', num2str(params(9),3));

clear defaultparamsBMOC
inputs = num2cell(params);
plotTC(tcnum, inputs{:}, toggleValues, ax)



%% Text input callback function
function text_Callback(~, ~, handles, ax)

toggleValues    = get(handles.toggleGroup, 'Value');
parameterValues = str2double(get(handles.textGroup, 'String'));
tcnum = get(handles.sliderTC,'Value');
params = num2cell(parameterValues);
plotTC(tcnum, params{:}, toggleValues, ax)


%% Plotting tuning curves function
function plotTC(tcnum, alphabm, beta1bm, delta1bm, alphaoc, beta1oc, delta1oc, c21, c12, thresh, toggleValues, ax)

% Load the Empirical Tuning Curve data.
load('JorisTCdata.mat');

% Define Reference Value
Fref    = 0.00002;  % Reference pressure in Pa.

tcnum = round(tcnum);


cf            = JorisTCdata{tcnum}(2,1);     % Center Frequency CF.
thresholdAtCF = JorisTCdata{tcnum}(1,1);
f0            = JorisTCdata{tcnum}(2,2:end); % Input/Stimulus Frequency
thresholds    = JorisTCdata{tcnum}(1,2:end);

idxcf         = f0 == cf;

f        = cf;
fsquared = f^2;

omega        =  2*pi*(f0 - f);
omegasquared = (2*pi*(f0 - f)).^2;
omega2       = omegasquared/fsquared;

MEFdBgains   = MEFdBgain(f0);  % In dB

pAt0dB = dB2Pa( MEFdBgain( [f-1 f f+1] ));
pAt0dB = pAt0dB(2);


if ~toggleValues{4}  % If single critical model
    
    if toggleValues{2}  % If frequency scaling
        F = thresh/c21*sqrt((alphaoc+beta1oc*thresh^2)^2 + (omega/f + delta1oc*thresh^2).^2);
    else  % Else not frequency scaling
        F = thresh/c21*sqrt((alphaoc+beta1oc*thresh^2)^2 + (omega + delta1oc*thresh^2).^2);
    end
    
else  % Else coupled oscillator model
    
    if toggleValues{3}  % If taking roc as threshold
        
        if toggleValues{2}  % If frequency scaling
            
            rbm = thresh/c21*sqrt((alphaoc + beta1oc*thresh^2)^2 + (omega/f + delta1oc*thresh^2).^2);
            sinPsioc = thresh*(omega/f + delta1oc*thresh^2)/(c21*rbm);
            cosPsioc = sqrt(1-sinPsioc.^2);
            F = sqrt((alphabm*rbm + c12*thresh*cosPsioc).^2 + (omega/f.*rbm + c12*thresh*sinPsioc).^2);    % Ji Chul's calculation for bidirectional coupling
            
        else  % Else not frequency scaling
            
            rbm = thresh/c21*sqrt((alphaoc + beta1oc*thresh^2)^2 + (omega + delta1oc*thresh^2).^2);
            sinPsioc = thresh*(omega + delta1oc*thresh^2)/(c21*rbm);
            cosPsioc = sqrt(1-sinPsioc.^2);
            F = sqrt((alphabm*rbm + c12*thresh*cosPsioc).^2 + (omega.*rbm + c12*thresh*sinPsioc).^2);    % Ji Chul's calculation for bidirectional coupling
            
        end
        
    else  % else taking rbm as threshold
        
        if toggleValues{2}  % If frequency scaling
            
            F = NaN(3,length(omega));
            for j = 1:length(omega)
                
                roc = sqrt(roots([beta1oc^2 + omega2(j), 2*(alphaoc*beta1oc + delta1oc*omega(j)/f), alphaoc^2+omega2(j), -(c21^2*thresh^2)]));
                cosPsioc = -(alphaoc*roc+beta1oc*roc.^3)/(c21*thresh);
                sinPsioc = (omega(j)*roc/f+delta1oc*roc.^3)/(c21*thresh);
                A = alphabm*thresh + c12*roc.*cosPsioc;
                B = omega(j)*thresh/f + c12*roc.*sinPsioc;
                temp = sqrt(A.^2 + B.^2);
                ind = abs(imag(temp)) < 10^(-5);
                temp = real(temp(ind));
                for i = 1:length(temp)
                    F(i,j) = temp(i);
                end
                
            end
            
        else  % Else not frequency scaling
            
            F = NaN(3,length(omega));
            for j = 1:length(omega)
                
                roc = sqrt(roots([beta1oc^2 + omegasquared(j), 2*(alphaoc*beta1oc + delta1oc*omega(j)), alphaoc^2+omegasquared(j), -(c21^2*thresh^2)]));
                cosPsioc = -(alphaoc*roc+beta1oc*roc.^3)/(c21*thresh);
                sinPsioc = (omega(j)*roc+delta1oc*roc.^3)/(c21*thresh);
                A = alphabm*thresh + c12*roc.*cosPsioc;
                B = omega(j)*thresh + c12*roc.*sinPsioc;
                temp = sqrt(A.^2 + B.^2);
                ind = abs(imag(temp)) < 10^(-5);
                temp = real(temp(ind));
                for i = 1:length(temp)
                    F(i,j) = temp(i);
                end
                
            end
            
        end
        
    end
    
end

FStable   = NaN(1,length(omega));
FUnstable = NaN(3,length(omega));

if ~toggleValues{4}  % Single critical model
    
    FStable = F;
    
else  % Coupled oscillator model
    
    if toggleValues{3}  % OC is threshold
        for j=1:length(omega)
            for i=1:size(F,1)
                if ~isnan(F(i,j))
                    if toggleValues{2}  % If frequency scaling
                        [~,r2,~,~,stable] = rStarCochlea(alphabm,alphaoc,beta1oc,delta1oc,F(i,j),c21,c12,omega(j)/f,1);
                    else               % Else not frequency scaling
                        [~,r2,~,~,stable] = rStarCochlea(alphabm,alphaoc,beta1oc,delta1oc,F(i,j),c21,c12,omega(j),1);
                    end
                    [~,index] = min(abs(thresh-r2));
                    if stable(index)
                        FStable(j)     = F(i,j);
                    else
                        FUnstable(i,j) = F(i,j);
                    end
                end
            end
        end
        
    else  % else BM is threshold
        
        for j=1:length(omega)
            for i=1:size(F,1)
                if ~isnan(F(i,j))
                    if toggleValues{2}  % If frequency scaling
                        [r1,~,~,~,stable] = rStarCochlea(alphabm,alphaoc,beta1oc,delta1oc,F(i,j),c21,c12,omega(j)/f,1);
                    else               % Else not frequency scaling
                        [r1,~,~,~,stable] = rStarCochlea(alphabm,alphaoc,beta1oc,delta1oc,F(i,j),c21,c12,omega(j),1);
                    end
                    [~,index] = min(abs(thresh-r1));
                    if stable(index)
                        FStable(j)     = F(i,j);
                    else
                        FUnstable(i,j) = F(i,j);
                    end
                end
            end
        end
        
    end
    
end

% Note: 1/Fref = 50000
if toggleValues{1}  % If middle ear filtering
    LstimStable   = 20*log10(FStable/Fref) - MEFdBgains;
    LstimUnstable = 20*log10(FUnstable/Fref) - repmat(MEFdBgains,3,1);
else               % Else not middle ear flitering
    LstimStable   = 20*log10(FStable/Fref);
    LstimUnstable = 20*log10(FUnstable/Fref);
end

semilogx(ax(1), f0, thresholds, 'b.-', 'Linewidth', 2);  hold(ax(1), 'on');
semilogx(ax(1), f0, LstimStable,   'Color', [.2 .8 .2], 'LineWidth', 2, 'Marker', '.');
semilogx(ax(1), f0, LstimUnstable, 'Color', [.8 .2 .8], 'LineWidth', 2, 'Marker', '.');

meanSquareError = mean((thresholds-LstimStable).^2);

xlim(ax(1), [50 30000]);
ylim(ax(1), [0 90]);
grid(ax(1), 'on');
ylabel(ax(1),'Threshold (dB SPL)','FontUnits','normalized');
xlabel(ax(1),'Forcing frequency (Hz)','FontUnits','normalized');
title(ax(1), ['TC ' num2str(tcnum) ' Fit  with  Tip = (' num2str(f) ', ' num2str(LstimStable(idxcf)) '), RMSE = ' num2str(sqrt(meanSquareError)) ' dB'],...
        'FontUnits','normalized');
set(ax(1), 'FontUnits','normalized');
hold(ax(1), 'off');


%% Amplitude Curves Method 2 
% Compute rocs by specifying Fs then use then use the rocs to compute
% rbms

F = logspace(log10(Fref),log10(20),201);

rBM = [];
rOC = [];
F_rBM_at_rOCthreshold = [];

clear f0;
n = 0;

threef0 = [f f/2 2*f];

for f0 = threef0
    n = n+1;
    
    rstarbm = NaN(size(F));
    rstaroc = NaN(size(F));
    for nF = 1:length(F)
        if ~toggleValues{4}  % If single critical model
            if toggleValues{2}  % If frequency scaling
                temp = sqrt(roots([beta1oc^2+delta1oc^2,...
                    2*(alphaoc*beta1oc+delta1oc*2*pi*(f-f0)/f),...
                    alphaoc^2+(2*pi*(f-f0)/f)^2,...
                    -F(nF)^2]));
                rstaroc(nF) = max(temp(~imag(temp)));
            else  % Else not frequency scaling
                temp = sqrt(roots([beta1oc^2+delta1oc^2,...
                    2*(alphaoc*beta1oc+delta1oc*2*pi*(f-f0)),...
                    alphaoc^2+(2*pi*(f-f0))^2,...
                    -F(nF)^2]));
                rstaroc(nF) = max(temp(~imag(temp)));
            end
        else  % Else coupled oscillator model
            if toggleValues{2}  % If frequency scaling
                [rbmn, rocn] = rStarCochlea(alphabm, alphaoc, beta1oc, delta1oc, ...
                    F(nF), c21, c12, 2*pi*(f-f0)/f);
            else  % Else not frequency scaling
                [rbmn, rocn] = rStarCochlea(alphabm, alphaoc, beta1oc, delta1oc, ...
                    F(nF), c21, c12, 2*pi*(f-f0));
            end
            if isempty(rbmn)
                rstarbm(nF) = NaN;
                rstaroc(nF) = NaN;
            else
                rstarbm(nF) = rbmn(1);
                rstaroc(nF) = rocn(1);
            end
        end
    end
    rBM = [rBM; rstarbm];
    rOC = [rOC; rstaroc];
end

clear rstaroc rstarbm

if c21&&c12   % both non-zero
    
    [spontBM, spontOC] = spontAmpCochlea(alphabm, alphaoc, beta1oc, delta1oc, c21, c12);
    
end

%% Plot Amplitude Curves 

F = Pa2dB(F); % For a log-log plot.
%F_rBM_at_rOCthreshold(:,1) = Pa2dB(F_rBM_at_rOCthreshold(:,1));

if ~toggleValues{4}  % If single critical oscillator
    
    semilogy(ax(2), F, rOC(1,:), 'g', F, rOC(2,:), 'b', F, rOC(3,:), 'r', [Fref 120], [thresh thresh], 'k--', 'Linewidth', 2); %, [0 max(F)] , [rbm rbm], 'k--');
    plot(ax(3),NaN);
    
else  % Else coupled oscillator model
    
    if toggleValues{3} % if oc==thresh
        
        if c21&&c12 % if afferent and efferent connections are both nonzero
            semilogy(ax(2), F, rOC(1,:), 'g', F, rOC(2,:), 'b', F, rOC(3,:), 'r', [Fref 120], [thresh thresh], 'k--', [Fref 120] , [spontOC spontOC], 'm--', 'Linewidth', 2); %, [0 max(F)] , [rbm rbm], 'k--');
            semilogy(ax(3), F, rBM(1,:), 'g', F, rBM(2,:), 'b', F, rBM(3,:), 'r', [Fref 120], [spontBM spontBM], 'm--', 'Linewidth', 2);
        else
            semilogy(ax(2), F, rOC(1,:), 'g', F, rOC(2,:), 'b', F, rOC(3,:), 'r', [Fref 120], [thresh thresh], 'k--', 'Linewidth', 2); %, [0 max(F)] , [rbm rbm], 'k--');
            semilogy(ax(3), F, rBM(1,:), 'g', F, rBM(2,:), 'b', F, rBM(3,:), 'r', 'Linewidth', 2);
        end
        
    else
        
        if c21&&c12
            semilogy(ax(2), F, rOC(1,:), 'g', F, rOC(2,:), 'b', F, rOC(3,:), 'r', [Fref 120], [spontOC spontOC], 'm--', 'Linewidth', 2); %, [0 max(F)] , [rbm rbm], 'k--');
            semilogy(ax(3), F, rBM(1,:), 'g', F, rBM(2,:), 'b', F, rBM(3,:), 'r', [Fref 120], [thresh thresh], 'k--', [Fref 120], [spontBM spontBM], 'm--', 'Linewidth', 2);
        else
            semilogy(ax(2), F, rOC(1,:), 'g', F, rOC(2,:), 'b', F, rOC(3,:), 'r', 'Linewidth', 2); %, [0 max(F)] , [rbm rbm], 'k--');
            semilogy(ax(3), F, rBM(1,:), 'g', F, rBM(2,:), 'b', F, rBM(3,:), 'r', [Fref 120], [thresh thresh], 'k--', 'Linewidth', 2);
        end
        
    end
    
end

xlabel(ax(2), 'Forcing F (dB SPL)','FontUnits','normalized');
ylabel(ax(2), 'r_{OC}^*','FontUnits','normalized');
title(ax(2), 'Steady State Amplitude Curve for OC Oscillator.',...
    'FontUnits','normalized');
set(ax(2), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
set(ax(2), 'FontUnits','normalized');
ax2ylim = get(ax(2), 'YLim');
if exist('rOC','var') && ~isnan(rOC(1,1))
    ylim(ax(2), [rOC(1,1)/10 ax2ylim(2)]);
end

grid(ax(3), 'on');

xlabel(ax(3), 'Forcing F (dB SPL)','FontUnits','normalized');
ylabel(ax(3), 'r_{BM}^*','FontUnits','normalized');
title(ax(3), 'Steady State Amplitude Curve for BM Oscillator.',...
    'FontUnits','normalized');
set(ax(3), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
set(ax(3), 'FontUnits','normalized');
ax3ylim = get(ax(3), 'YLim');
if exist('rBM', 'var') && ~isnan(rBM(1,1))
    ylim(ax(3), [rBM(1,1)/10 ax3ylim(2)]);
end

%% OC & BM Compression Curves 

f0sOrig = logspace(log10(50),log10(30000),201)';
temp    = f0sOrig - f;
[~,ind] = min(abs(temp));
f0s     = f0sOrig - temp(ind);
F       = dB2Pa([120 100 80 60 40 20 0]);
if toggleValues{2}  % If frequency scaling
    W = 2*pi*(f-f0s)/f;
else               % Else not frequency scaling
    W = 2*pi*(f-f0s);
end
BMcomp = zeros(length(W),length(F));
OCcomp = zeros(length(W),length(F));

for nF = 1:length(F)
    for nW = 1:length(W)
        if ~toggleValues{4}  % If single critical oscillator
            temp = sqrt(roots([beta1oc^2+delta1oc^2,...
                2*(alphaoc*beta1oc+delta1oc*W(nW)),...
                alphaoc^2+W(nW)^2,...
                -F(nF)^2]));
            OCcomp(nW,nF) = max(temp(~imag(temp)));
        else  % Else coupled oscillator model
            
            [rbmn, rocn] = rStarCochlea(alphabm, alphaoc, beta1oc, delta1oc, ...
                F(nF), c21, c12, W(nW));
            if isempty(rbmn)
                BMcomp(nW,nF) = NaN;
                OCcomp(nW,nF) = NaN;
            else
                BMcomp(nW,nF) = rbmn(1);
                OCcomp(nW,nF) = rocn(1);
            end
        end
    end
end

if ~toggleValues{4}  % If single critical oscillator
    
    loglog(ax(4), f0s(f0s>0), OCcomp(f0s>0,:), 'Linewidth', 2); hold(ax(4), 'on');
    loglog(ax(4), [min(f0sOrig) max(f0sOrig)], [thresh thresh], 'k--', 'Linewidth', 2); hold(ax(4), 'off');
    loglog(ax(5), NaN);
    
else
    
    if toggleValues{3}  % If using OC as threshold
        
        if c21&&c12
            loglog(ax(5), f0s(f0s>0), BMcomp(f0s>0,:), 'Linewidth', 2);hold(ax(5), 'on');
            loglog(ax(5), [min(f0sOrig) max(f0sOrig)], [spontBM spontBM], 'm--', 'Linewidth', 2);hold(ax(5), 'off');
            loglog(ax(4), f0s(f0s>0), OCcomp(f0s>0,:), 'Linewidth', 2); hold(ax(4), 'on');
            loglog(ax(4), [min(f0sOrig) max(f0sOrig)], [thresh thresh], 'k--', 'Linewidth', 2);
            loglog(ax(4), [min(f0sOrig) max(f0sOrig)], [spontOC spontOC], 'm--', 'Linewidth', 2); hold(ax(4), 'off');
        else
            loglog(ax(5), f0s(f0s>0), BMcomp(f0s>0,:), 'Linewidth', 2);
            loglog(ax(4), f0s(f0s>0), OCcomp(f0s>0,:), 'Linewidth', 2); hold(ax(4), 'on');
            loglog(ax(4), [min(f0sOrig) max(f0sOrig)], [thresh thresh], 'k--', 'Linewidth', 2); hold(ax(4), 'off');
        end
        
    else  % Else using BM as threshold
        
        if c21&&c12
            loglog(ax(5), f0s(f0s>0), BMcomp(f0s>0,:), 'Linewidth', 2);hold(ax(5), 'on');
            loglog(ax(5), [min(f0sOrig) max(f0sOrig)], [thresh thresh], 'k--', 'Linewidth', 2);
            loglog(ax(5), [min(f0sOrig) max(f0sOrig)], [spontBM spontBM], 'm--', 'Linewidth', 2);hold(ax(5), 'off');
            loglog(ax(4), f0s(f0s>0), OCcomp(f0s>0,:), 'Linewidth', 2); hold(ax(4), 'on');
            loglog(ax(4), [min(f0sOrig) max(f0sOrig)], [spontOC spontOC], 'm--', 'Linewidth', 2); hold(ax(4), 'off');
        else
            loglog(ax(5), f0s(f0s>0), BMcomp(f0s>0,:), 'Linewidth', 2);hold(ax(5), 'on');
            loglog(ax(5), [min(f0sOrig) max(f0sOrig)], [thresh thresh], 'k--', 'Linewidth', 2); hold(ax(5), 'off');
            loglog(ax(4), f0s(f0s>0), OCcomp(f0s>0,:), 'Linewidth', 2);
        end
        
    end
    
end


xlim(ax(5), [min(f0sOrig) max(f0sOrig)]);
set(ax(5), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
xlabel(ax(5), 'Forcing frequency (Hz)','FontUnits','normalized');
ylabel(ax(5), 'r_{BM}^*','FontUnits','normalized');
title(ax(5), 'BM Compression Curves','FontUnits','normalized');
set(ax(5), 'FontUnits','normalized');

xlim(ax(4), [min(f0sOrig) max(f0sOrig)]);

ylim(ax(4), [thresh/10 2*max(rOC(1,:))] );
set(ax(4), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
xlabel(ax(4), 'Forcing frequency (Hz)','FontUnits','normalized');
ylabel(ax(4), 'r_{OC}^*','FontUnits','normalized');
title(ax(4), 'OC Compression Curves','FontUnits','normalized');
set(ax(4), 'FontUnits','normalized');

%% Decide what needs to be plotted based on toggle button values

stringL1 = {'120','100','80','60','40','20','0dB','thresh'};
stringL2 = {'f0 = CF = f','f0 = f/2','f0 = 2f'};
    
if c21&&c12 % if afferent and efferent connections are both nonzero 
    stringL1 = [stringL1 {'Spont Amp'}];
end

hl2 = legend(ax(2),stringL2,'Orientation', 'horizontal','fontsize', 12,...
    'units', 'normalized');
hl4 = legend(ax(4),stringL1,'Orientation', 'horizontal','fontsize', 12,...
    'units', 'normalized');

% Set position, orientation of legends 
temp2 = get(hl2, 'Position'); hl2size = temp2(3:4);
temp4 = get(hl4, 'Position'); hl4size = temp4(3:4);

set(hl4, 'Position', [.02 .28 hl4size(:)'])
set(hl2, 'Position', [.7 .28 hl2size(:)'])

