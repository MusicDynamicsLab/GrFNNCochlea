function TCfitGUI
% function TCfitGUI - GUI for fitting model tuning curves to empirical
% results. Change parameter values by moving the sliders or by entering
% numbers in the text boxes.

figure(999);
set(gcf, 'Visible', 'off', 'Toolbar', 'figure', 'Color', [.8 .8 .8], ...
    'Position', [5, 5, 900 800]);
initialRun;
movegui(999, 'center');
set(999, 'Visible', 'on');



%% Initial run function for when GUI first is created
function initialRun(~, ~)
clf
axis off;

ax1 = axes('Units', 'normalized', 'Position', [.06 .65 .5 .3]);
ax2 = axes('Units', 'normalized', 'Position', [.63  288/800 .35 .2]);
ax3 = axes('Units', 'normalized', 'Position', [.63   36/800 .35 .2]);
ax4 = axes('Units', 'normalized', 'Position', [.06  288/800 .5 .2]);
ax5 = axes('Units', 'normalized', 'Position', [.06   36/800 .5 .2]);
ax = [ax1 ax2 ax3 ax4 ax5];

xlim(ax1, [50 20000]);
ylim(ax1, [0 100]);
grid(ax1, 'on');
ylabel(ax1, 'Threshold (dB SPL)');
xlabel(ax1, 'Frequency (Hz)');
title(ax1, 'TC Fit');

xlabel(ax2, 'F (Pa)');
ylabel(ax2, 'r* of OC oscillator');
title(ax2, 'Steady State Amplitude Curve for OC Oscillator.');


xlabel(ax3, 'Forcing F (Pa) = r* OC');
ylabel(ax3, 'r* of BM oscillator');
title(ax3, 'Steady State Amplitude Curve for BM Oscillator.');

xlabel(ax4, 'f_{0} (Hz)');
ylabel(ax4, 'rOC^*');
title(ax4, 'OC Compression Curves: 0 - 120 dB SPL in 20 dB steps');

xlabel(ax5, 'f_{0} (Hz)');
ylabel(ax5, 'rBM^*');
title(ax5, 'BM Compression Curves: 0 - 120 dB SPL in 20 dB steps');

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

% Tuning curve slider, textbox, and label
text(.71,.62,'Tuning Curve:','HorizontalAlignment','right','VerticalAlignment','baseline',...
    'Units','normalized','FontSize',15)
handles.sliderTC = uicontrol('Style', 'slider', 'Min', 1, 'Max', 11, 'Value', tcnum, 'SliderStep', [.1 .1], ...
    'Min', TCnumLim(1), 'Max', TCnumLim(2), 'Units', 'normalized', 'Position', [.71 .59 .222 .025]);
handles.textTC = uicontrol('Style','edit', 'String', tcnum,...
    'FontSize',12,'FontUnits','normalized','Units', 'normalized','Position', [.81 .62 .046 .03]);

% Panel for parameter input
hp1 = uipanel('Title','Parameters','FontSize',11,'FontName','Helvetica',...
    'BackgroundColor',[.8 .8 .8],'Position',[.67 .66 .28 .3]);
axes('Units','normalized','Position',get(hp1,'Position'))
axis off

% Label and textbox for alpha_BM
text(.18,.83,'$$\alpha_{bm}:$$','Interpreter','Latex',...
    'HorizontalAlignment','right','VerticalAlignment','baseline',...
    'Units','normalized','FontSize',15)
handles.alphabm = uicontrol('Parent',hp1,'Style','edit', 'String', alphabm,...
    'FontSize',12,'Units', 'normalized','Position', [.19 .83 .2 .13]);
% Label and textbox for beta_BM
text(.18,.68,'$$\beta_{bm}:$$','Interpreter','Latex',...
    'HorizontalAlignment','right','VerticalAlignment','baseline',...
    'Units','normalized','FontSize',15)
handles.betabm = uicontrol('Parent',hp1,'Style','edit', 'String', beta1bm,...
    'FontSize',12,'Units', 'normalized','Position', [.19 .68 .2 .13],'enable','off');
% Label and textbox for delta_BM
text(.18,.55,'$$\delta_{bm}:$$','Interpreter','Latex',...
    'HorizontalAlignment','right','VerticalAlignment','baseline',...
    'Units','normalized','FontSize',15)
handles.deltabm = uicontrol('Parent',hp1,'Style','edit', 'String', delta1bm,...
    'FontSize',12,'Units', 'normalized','Position', [.19 .53 .2 .13],'enable','off');
% Label and textbox for c12
text(.18,.4,'$$c_{12}:$$','Interpreter','Latex',...
    'HorizontalAlignment','right','VerticalAlignment','baseline',...
    'Units','normalized','FontSize',15)
handles.c12 = uicontrol('Parent',hp1,'Style','edit', 'String', c12,...    
    'FontSize',12,'Units', 'normalized','Position', [.19 .38 .2 .13]);

% Label and textbox for alpha_OC
text(.54,.82,'$$\alpha_{oc}:$$','Interpreter','Latex',...
    'HorizontalAlignment','right','VerticalAlignment','baseline',...
    'Units','normalized','FontSize',15)
handles.alphaoc = uicontrol('Parent',hp1,'Style','edit', 'String', alphaoc,...
    'FontSize',12,'Units', 'normalized','Position', [.55 .83 .22 .13]);
% Label and textbox for beta_OC
text(.54,.67,'$$\beta_{oc}:$$','Interpreter','Latex',...
    'HorizontalAlignment','right','VerticalAlignment','baseline',...
    'Units','normalized','FontSize',15)
handles.betaoc = uicontrol('Parent',hp1,'Style','edit', 'String', beta1oc,...
    'FontSize',12,'Units', 'normalized','Position', [.55 .68 .22 .13]);
% Label and textbox for delta_OC
text(.54,.52,'$$\delta_{oc}:$$','Interpreter','Latex',...
    'HorizontalAlignment','right','VerticalAlignment','baseline',...
    'Units','normalized','FontSize',15)
handles.deltaoc = uicontrol('Parent',hp1,'Style','edit', 'String', delta1oc,...
    'FontSize',12,'Units', 'normalized','Position', [.55 .53 .22 .13]);
% Label and texbox for c21
text(.54,.37,'$$c_{21}:$$','Interpreter','Latex',...
    'HorizontalAlignment','right','VerticalAlignment','baseline',...
    'Units','normalized','FontSize',15)
handles.c21 = uicontrol('Parent',hp1,'Style','edit', 'String', c21,...
    'FontSize',12,'Units', 'normalized','Position', [.55 .38 .22 .13]);

% Label and textbox for threshold
text(.3,.08,'Threshold:','HorizontalAlignment','right',...
    'VerticalAlignment','baseline','Units','normalized','FontSize',15)
handles.thresh = uicontrol('Parent', hp1, 'Style','edit', 'String', thresh,...
    'FontSize',12,'Units', 'normalized', 'HorizontalAlignment', 'center',...
    'Position', [.31 .05 .22 .13]);

% Pushbutton controls for saving as default parameters and choosing c12
% value
handles.buttonSaveAsDefault = uicontrol('Parent', hp1, 'Style', 'pushbutton', 'String', {'Save as Default'},...
    'TooltipString', ['Save the current parameters' char(10)...
    'as the default parameters in' char(10) 'defaultparamsBMOC.csv'],...
    'FontUnits','normalized','Units', 'normalized',  'Position', [.58 .02 .4 .18]);
handles.buttonChoosec12 = uicontrol('Parent', hp1, 'Style', 'pushbutton', 'String', 'Choose c12',...
    'FontUnits','normalized','Units', 'normalized',...    
    'TooltipString', ['Set value of c12 to make spontaneous' char(10) 'amplitude equal to half the threshold.'],...
    'Position', [.10 .2 .28 .15]);

% Create panel for toggle controls
hp2 = uipanel('Title','Toggles','FontSize',11,'FontName','Helvetica',...
    'BackgroundColor',[.8 .8 .8],'Position',[.57 .66 .09 .3]);
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
handles.toggleThresh = uicontrol('Parent',hp2,'Style','togglebutton', 'String','OC==thresh',...
    'TooltipString', ['ON = use OC amplitude as threshold ' char(10) 'OFF = use BM amplitude as threshold'],...
    'Value', 1,'FontUnits','normalized','Units', 'normalized', 'Position', [.03 .3 .95 .15]);
handles.toggleSingle = uicontrol('Parent',hp2,'Style','togglebutton', 'String','Single Osc',...
    'TooltipString', ['ON = single layer model ' char(10) 'OFF = two layer model'],...
    'Value', 0,'FontUnits','normalized','Units', 'normalized', 'Position', [.03 .05 .95 .15]);


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
plotTC(tcnum, alphabm, beta1bm, delta1bm, alphaoc, beta1oc, delta1oc, c21, c12, thresh, {1, 0, 1, 0}, ax);

%% Toggle pushbutton callback function
function toggle_Callback(~, ~, handles, ax)

defaultparamsBMOC = csvread('defaultparamsBMOC.csv');
tcnum = get(handles.sliderTC,'Value');
freqScaling    = get(handles.toggleFS, 'Value');
threshOC       = get(handles.toggleThresh, 'Value');
singleCritical = get(handles.toggleSingle, 'Value');

if singleCritical
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

if singleCritical
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

if ~singleCritical
    
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


%%  Tuning curve callback function
function Callback_TC(source, ~, handles, ax)

if strcmp(get(source,'Style'),'edit')
    tcnum = str2double(get(source,'String'));
    set(handles.sliderTC, 'Value', tcnum)
else
    tcnum = get(source,'Value');
    set(handles.textTC,  'String', num2str(tcnum));
end

defaultparamsBMOC = csvread('defaultparamsBMOC.csv');
toggleValues      = get(handles.toggleGroup, 'Value');

if toggleValues{4}
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



% =========================================================================
function text_Callback(~, ~, handles, ax)

toggleValues    = get(handles.toggleGroup, 'Value');
parameterValues = str2double(get(handles.textGroup, 'String'));
tcnum = get(handles.sliderTC,'Value');
params = num2cell(parameterValues);
plotTC(tcnum, params{:}, toggleValues, ax)


% =========================================================================
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


if toggleValues{4}  % If single critical model
    
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

if toggleValues{4}  % Single critical model
    
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

xlim(ax(1), [50 30000]);
ylim(ax(1), [0 90]);
grid(ax(1), 'on');
ylabel(ax(1),'Threshold (dB SPL)');
xlabel(ax(1),'Forcing frequency (Hz)');
title(ax(1), ['TC ' num2str(tcnum) ' Fit  with  Tip = (' num2str(f) ', ' num2str(LstimStable(idxcf)) ')']);
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
        if toggleValues{4}  % If single critical model
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

if toggleValues{4}  % If single critical oscillator
    
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

xlabel(ax(2), 'Forcing F (dB SPL)');
ylabel(ax(2), 'r* of OC oscillator');
title(ax(2), 'Steady State Amplitude Curve for OC Oscillator.');
set(ax(2), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
set(ax(2), 'FontUnits','normalized');

grid(ax(3), 'on');

xlabel(ax(3), 'Forcing F (dB SPL)');
ylabel(ax(3), 'r* of BM oscillator');
title(ax(3), 'Steady State Amplitude Curve for BM Oscillator.');
set(ax(3), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
set(ax(3), 'FontUnits','normalized');




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
        if toggleValues{4}  % If single critical oscillator
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

if toggleValues{4}  % If single critical oscillator
    
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
xlabel(ax(5), 'Forcing frequency (Hz)');
ylabel(ax(5), 'rBM^*');
title(ax(5), 'BM Compression Curves');
set(ax(5), 'FontUnits','normalized');

xlim(ax(4), [min(f0sOrig) max(f0sOrig)]);

ylim(ax(4), [thresh/10 2*max(rOC(1,:))] );
set(ax(4), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
xlabel(ax(4), 'Forcing frequency (Hz)');
ylabel(ax(4), 'rOC^*');
title(ax(4), 'OC Compression Curves');
set(ax(4), 'FontUnits','normalized');

%% Decide what needs to be plotted based on toggle button values

stringL2 = {'f0 = CF = f','f0 = f/2','f0 = 2f'};
stringL3 = {'f0 = CF = f','f0 = f/2','f0 = 2f'};
stringL4 = {'120 dB SPL','100 dB SPL','80 dB SPL','60 dB SPL',...
    '40 dB SPL','20 dB SPL','0 dB SPL'};
stringL5 = {'120 dB SPL','100 dB SPL','80 dB SPL','60 dB SPL',...
    '40 dB SPL','20 dB SPL','0 dB SPL'};

if toggleValues{3} % if oc==thresh
    stringL2 = [stringL2 {'Threshold'}];
    stringL4 = [stringL4 {'Threshold'}];
else % else bm==thresh
    stringL3 = [stringL3 {'Threshold'}];
    stringL5 = [stringL5 {'Threshold'}];
end
if c21&&c12 % if afferent and efferent connections are both nonzero 
    stringL2 = [stringL2 {'Spontaneous amp'}];
    stringL3 = [stringL3 {'Spontaneous amp'}];
    stringL4 = [stringL4 {'Spontaneous amp'}];
    stringL5 = [stringL5 {'Spontaneous amp'}];
end

if toggleValues{4} % if using single oscillator model
    hl2 = legend(ax(2),stringL2,'location','best');
    hl4 = legend(ax(4),stringL4);
    
    if tcnum < 5
        set(hl4,'Location','northeast');
    else
        set(hl4,'Location','northwest');
    end
else
    hl2 = legend(ax(2),stringL2,'location','best');
    hl3 = legend(ax(3),stringL3,'location','best');
    hl4 = legend(ax(4),stringL4);
    hl5 = legend(ax(5),stringL5);
    
    if tcnum < 5
        set([hl4,hl5],'Location','northeast');
    else
        set([hl4,hl5],'Location','northwest');
    end
end

% Set font size of legend text
set([hl2 hl3 hl4 hl5],'FontSize',12);



