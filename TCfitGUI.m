function TCfitGUI
% function TCfitGUI
% GUI for fiiting model tuning curves to empirical results.
% Change parameter values by moving the sliders or by entering numbers
% in the text boxes. Min and max for the sliders can also be changed.

% clear all;
% close all;
% clc;

%global MEFdBgains thresholds omegasquared f f0;

figure(999);
set(gcf, 'Visible', 'off', 'Toolbar', 'figure', 'Color', [.8 .8 .8], ...
    'Position', [5, 5, 1200 800]);
initialRun;
movegui(999, 'center');
set(999, 'Visible', 'on');



% =========================================================================
function initialRun(source, eventdata)
clf
axis off;

ax1 = axes('Units', 'normalized', 'Position', [ 50/1200  520/800 600/1200 240/800]);
ax2 = axes('Units', 'normalized', 'Position', [ 50/1200  288/800 400/1200 160/800]);
ax3 = axes('Units', 'normalized', 'Position', [ 50/1200   36/800 400/1200 160/800]);
ax4 = axes('Units', 'normalized', 'Position', [520/1200  288/800 480/1200 160/800]);
ax5 = axes('Units', 'normalized', 'Position', [520/1200   36/800 480/1200 160/800]);


ax = [ax1 ax2 ax3 0 ax4 ax5];


xlim(ax1, [20 20000]);
ylim(ax1, [0 100]);
grid(ax1, 'on');
ylabel(ax1, 'Threshold (dB SPL)');
xlabel(ax1, 'Frequency (Hz)');
title(ax1, 'TC Fit');

xlabel(ax2, 'F (Pa)');
ylabel(ax2, 'r* of OC oscillator');
title(ax2, {['Steady State Amplitude Curve for OC Oscillator.'], ['Blue @ f0 = CF = f, Green @ f0 = f/2,  Red @ f0 = 2f,  Black = rBM Threshold']});

xlabel(ax3, 'Forcing F (Pa) = r* OC');
ylabel(ax3, 'r* of BM oscillator');
title(ax3, {['Steady State Amplitude Curve for BM Oscillator.'], ['Blue @ f0 = CF = f, Green @ f0 = f/2,  Red @ f0 = 2f,  Black = rBM Threshold']});

xlabel(ax4, 'f_{0} (Hz)');
ylabel(ax4, 'rOC^*');
title(ax4, 'OC Compression Curves: 0 - 120 dB SPL in 20 dB steps');

xlabel(ax5, 'f_{0} (Hz)');
ylabel(ax5, 'rBM^*');
title(ax5, 'BM Compression Curves: 0 - 120 dB SPL in 20 dB steps');

defaultparamsBMOC = csvread('defaultparamsBMOC.csv');

% Initial parameter values
tcnum       = 5;
alphaoc     = defaultparamsBMOC(9,tcnum);
alphabm     = defaultparamsBMOC(10,tcnum);
beta1oc     = defaultparamsBMOC(11,tcnum);
delta       = defaultparamsBMOC(12,tcnum);
c21         = defaultparamsBMOC(13,tcnum);
c12         = defaultparamsBMOC(14,tcnum);
thresh      = defaultparamsBMOC(15,tcnum);
TCfreqshift = defaultparamsBMOC(16,tcnum);
Lstim       = 0; % Stimulus/input level in dB SPL


clear defaultparamsBMOC


TCnumLim   = [1 11];


slider1 = uicontrol('Style', 'slider', 'Min', 1, 'Max', 11, 'Value', tcnum, 'SliderStep', [.1 .1], ...
    'Min', TCnumLim(1), 'Max', TCnumLim(2), 'Units', 'normalized', 'Position', [980/1200 704/800 200/1200 20/800]);


text1 = uicontrol('Style','edit', 'String', alphaoc,...
    'FontSize',12,'FontUnits','normalized','Units', 'normalized','Position', [1056/1200 600/800 56/1200 20/800]);
text2 = uicontrol('Style','edit', 'String', alphabm,...
    'FontSize',12,'FontUnits','normalized','Units', 'normalized','Position', [1056/1200 520/800 56/1200 20/800]);
text3 = uicontrol('Style','edit', 'String', beta1oc,...
    'FontSize',12,'FontUnits','normalized','Units', 'normalized','Position', [1056/1200 440/800 56/1200 20/800]);
text4 = uicontrol('Style','edit', 'String', delta,...
    'FontSize',12,'FontUnits','normalized','Units', 'normalized','Position', [1056/1200 360/800 56/1200 20/800]);
text5 = uicontrol('Style','edit', 'String', c21,...
    'FontSize',12,'FontUnits','normalized','Units', 'normalized','Position', [1056/1200 280/800 56/1200 20/800]);
text6 = uicontrol('Style','edit', 'String', c12,...
    'FontSize',12,'FontUnits','normalized','Units', 'normalized','Position', [1056/1200 200/800 56/1200 20/800]);
text7 = uicontrol('Style','edit', 'String', thresh,...
    'FontSize',12,'FontUnits','normalized','Units', 'normalized','Position', [1056/1200 120/800 56/1200 20/800]);
text8 = uicontrol('Style','edit', 'String', TCfreqshift,...
    'FontSize',12,'FontUnits','normalized','Units', 'normalized','Position', [1056/1200  32/800 56/1200 20/800]);
text9 = uicontrol('Style','edit', 'String', tcnum,...
    'FontSize',12,'FontUnits','normalized','Units', 'normalized','Position', [1056/1200 680/800 56/1200 20/800]);


toggle1 = uicontrol('Style','togglebutton', 'String','Middle ear filter On/Off',...
    'Value', 1,'FontUnits','normalized','Units', 'normalized', 'Position', [660/1200 740/800 240/1200 20/800]);

toggle2 = uicontrol('Style','togglebutton', 'String','Freq Scaling On/Off',...
    'Value', 0,'FontUnits','normalized','Units', 'normalized', 'Position', [660/1200 680/800 240/1200 20/800]);

toggle3 = uicontrol('Style','togglebutton', 'String','On for OC==thresh, off for BM==thresh',...
    'Value', 1,'FontUnits','normalized','Units', 'normalized', 'Position', [660/1200 620/800 240/1200 20/800]);

toggle4 = uicontrol('Style','togglebutton', 'String','On for single oscillator model',...
    'Value', 0,'FontUnits','normalized','Units', 'normalized', 'Position', [660/1200 560/800 240/1200 20/800]);


label1 = uicontrol('Style','text','String', {'Alpha OC'},...
    'BackgroundColor',[.8 .8 .9],'FontSize', 12,'FontUnits','normalized','Units', 'normalized',  'Position', [1040/1200 620/800 80/1200 20/800]);
label2 = uicontrol('Style','text','String', {'Alpha BM'},...
    'BackgroundColor',[.8 .8 .9],'FontSize', 12,'FontUnits','normalized','Units', 'normalized',  'Position', [1040/1200 540/800 80/1200 20/800]);
label3 = uicontrol('Style','text','String', {'Beta'},...
    'BackgroundColor',[.8 .8 .9],'FontSize', 12,'FontUnits','normalized','Units', 'normalized',  'Position', [1040/1200 460/800 80/1200 20/800]);
label4 = uicontrol('Style','text','String', {'Delta'},...
    'BackgroundColor',[.8 .8 .9],'FontSize', 12,'FontUnits','normalized','Units', 'normalized',  'Position', [1040/1200 380/800 80/1200 20/800]);
label5 = uicontrol('Style','text','String', {'c21'},...
    'BackgroundColor',[.8 .8 .9],'FontSize', 12,'FontUnits','normalized','Units', 'normalized',  'Position', [1040/1200 300/800 80/1200 20/800]);
label6 = uicontrol('Style','text','String', {'c12'},...
    'BackgroundColor',[.8 .8 .9],'FontSize', 12,'FontUnits','normalized','Units', 'normalized',  'Position', [1040/1200 220/800 80/1200 20/800]);
label7 = uicontrol('Style','text','String', {'Threshold'},...
    'BackgroundColor',[.8 .8 .9],'FontSize', 12,'FontUnits','normalized','Units', 'normalized',  'Position', [1040/1200 140/800 80/1200 20/800]);
label8 = uicontrol('Style','text','String', {'TC Freq. Shift (HZ)'},...
    'BackgroundColor',[.8 .8 .9],'FontSize', 12,'FontUnits','normalized','Units', 'normalized',  'Position', [1040/1200 60/800 80/1200 28/800]);
label9 = uicontrol('Style','text','String', {'TC Number'},...
    'BackgroundColor',[.8 .8 .9],'FontSize', 12,'FontUnits','normalized','Units', 'normalized',  'Position', [1040/1200 724/800 80/1200 20/800]);


pushbutton1 = uicontrol('Style', 'pushbutton', 'String', 'Save as Default Params',...
    'FontUnits','normalized','Units', 'normalized',  'Position', [660/1200 500/800 240/1200 20/800]);
pushbutton2 = uicontrol('Style', 'pushbutton', 'String', 'Choose',...
    'FontUnits','normalized','Units', 'normalized',  'Position', [1116/1200 200/800 62/1200 20/800]);



% Create an array of handles

handles = [ 0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0;
    slider1 0    0;
    text1   0    0;
    text2   0    0;
    text3   0    0;
    text4   0    0;
    text5   0    0;
    text6   0    0;
    text7   0    0;
    text8   0    0;
    text9   0    0;
    toggle1 0    1;
    toggle2 0    1;
    toggle3 0    1;
    toggle4 0    1];


% Set callback functions
set(handles(9,1),    'Callback', {@slider_Callback_TCNUM, handles, ax});
set(handles(10:17,1), 'Callback', {@text_Callback,   handles, ax});
set(handles(18,1),   'Callback', {@text_CallbackTCNUM,  handles, ax});
set(pushbutton1,     'Callback', {@pushbutton1_Callback, handles});
set(pushbutton2,     'Callback', {@pushbutton2_Callback, handles, ax});
set(handles(19,1), 'Callback',{@toggle_Callback,handles,ax});
set(handles(20,1), 'Callback',{@toggle_Callback,handles,ax});
set(handles(21,1), 'Callback',{@toggle_Callback,handles,ax});
set(handles(22,1), 'Callback',{@toggle_Callback,handles,ax});


hold off
plotTC(tcnum, alphaoc, alphabm, beta1oc, delta, c21, c12, TCfreqshift, thresh, {1, 0, 1, 0}, ax);


% =========================================================================
function toggle_Callback(source, eventdata, handles, ax)

defaultparamsBMOC = csvread('defaultparamsBMOC.csv');
tcnum = get(handles(9,1),'Value');
freqScaling    = get(handles(20,1), 'Value');
threshOC       = get(handles(21,1), 'Value');
singleCritical = get(handles(22,1), 'Value');

if singleCritical
    if freqScaling
        params = defaultparamsBMOC(33:40,tcnum);
    else
        params = defaultparamsBMOC(41:48,tcnum);
    end
else
    
    if threshOC
        if freqScaling
            params = defaultparamsBMOC(1:8, tcnum);
        else
            params = defaultparamsBMOC(9:16,tcnum);
        end
    else
        if freqScaling
            params = defaultparamsBMOC(17:24, tcnum);
        else
            params = defaultparamsBMOC(25:32,tcnum);
        end
    end
    
end

set(handles(10,1),  'String', num2str(params(1),3));
set(handles(11,1),  'String', num2str(params(2),3));
set(handles(12,1),  'String', num2str(params(3),3));
set(handles(13,1),  'String', num2str(params(4),3));
set(handles(14,1),  'String', num2str(params(5),3));
set(handles(15,1),  'String', num2str(params(6),3));
set(handles(16,1),  'String', num2str(params(7),3));
set(handles(17,1),  'String', num2str(params(8),3));
set(handles(18,1),  'String', num2str(tcnum));

toggleValues    = get(handles(19:22,1), 'Value');
parameterValues = str2double(get(handles(10:18,1), 'String'));
clear defaultparamsBMOC

plotTC(parameterValues(9), parameterValues(1), parameterValues(2), parameterValues(3), ...
    parameterValues(4), parameterValues(5), parameterValues(6), parameterValues(8), parameterValues(7), toggleValues, ax)




% =========================================================================
function pushbutton1_Callback(source, eventdata, handles)

pathToParams = which('defaultparamsBMOC.csv');
defaultparamsBMOC = csvread('defaultparamsBMOC.csv');

freqScaling    = get(handles(20,1), 'Value');
threshOC       = get(handles(21,1), 'Value');
singleCritical = get(handles(22,1), 'Value');

params = str2double(get(handles(10:17,1), 'String'));

tcnum  = str2double(get(handles(18,1), 'String'));

if singleCritical
    if freqScaling
        defaultparamsBMOC(33:40,tcnum) = params;
    else
        defaultparamsBMOC(41:48,tcnum) = params;
    end
else
    
    if threshOC
        if freqScaling
            defaultparamsBMOC(1:8,tcnum)  = params;
        else
            defaultparamsBMOC(9:16,tcnum) = params;
        end
    else
        if freqScaling
            defaultparamsBMOC(17:24,tcnum)  = params;
        else
            defaultparamsBMOC(25:32,tcnum) = params;
        end
    end
    
end

csvwrite(pathToParams,defaultparamsBMOC);
clear defaultparamsBMOC


% =========================================================================
function pushbutton2_Callback(source, eventdata, handles, ax)

singleCritical = get(handles(22,1), 'Value');

params = str2double(get(handles(10:17,1), 'String'));

if ~singleCritical
    
    aoc   =params(1);
    abm   =params(2);
    boc   =params(3);
    c21   =params(5);
    thresh=params(7);
    factor=.5;
    
    c12=abm/c21*(aoc+boc*(thresh*factor)^2);
    
    textValue = num2str(c12);
    
    set(handles(15,1),'String',textValue);
    
    toggleValues    = get(handles(19:22,1), 'Value');
    parameterValues = str2double(get(handles(10:18,1), 'String'));
    
    plotTC(parameterValues(9), parameterValues(1), parameterValues(2), parameterValues(3), ...
        parameterValues(4), parameterValues(5), parameterValues(6), parameterValues(8), parameterValues(7), toggleValues, ax)
end


% =========================================================================
function slider_Callback_TCNUM(source, eventdata, handles, ax)

defaultparamsBMOC = csvread('defaultparamsBMOC.csv');
tcnum = get(source,'Value');
freqScaling    = get(handles(20,1), 'Value');
threshOC       = get(handles(21,1), 'Value');
singleCritical = get(handles(22,1), 'Value');

if singleCritical
    if freqScaling
        params = defaultparamsBMOC(33:40,tcnum);
    else
        params = defaultparamsBMOC(41:48,tcnum);
    end
else
    
    if threshOC
        if freqScaling
            params = defaultparamsBMOC(1:8, tcnum);
        else
            params = defaultparamsBMOC(9:16,tcnum);
        end
    else
        if freqScaling
            params = defaultparamsBMOC(17:24, tcnum);
        else
            params = defaultparamsBMOC(25:32,tcnum);
        end
    end
    
end

set(handles(9,1),   'Value',  tcnum);
set(handles(10,1),  'String', num2str(params(1),3));
set(handles(11,1),  'String', num2str(params(2),3));
set(handles(12,1),  'String', num2str(params(3),3));
set(handles(13,1),  'String', num2str(params(4),3));
set(handles(14,1),  'String', num2str(params(5),3));
set(handles(15,1),  'String', num2str(params(6),3));
set(handles(16,1),  'String', num2str(params(7),3));
set(handles(17,1),  'String', num2str(params(8),3));
set(handles(18,1),  'String', num2str(tcnum));

toggleValues    = get(handles(19:22,1), 'Value');
parameterValues = str2double(get(handles(10:18,1), 'String'));
clear defaultparamsBMOC

plotTC(parameterValues(9), parameterValues(1), parameterValues(2), parameterValues(3), ...
    parameterValues(4), parameterValues(5), parameterValues(6), parameterValues(8), parameterValues(7), toggleValues, ax)




% =========================================================================
function text_Callback(source, eventdata, handles, ax)
[row,col] = find(handles == source);

textValue = get(source,'String');

set(handles(row,1),'String',textValue);

toggleValues    = get(handles(19:22,1), 'Value');
parameterValues = str2double(get(handles(10:18,1), 'String'));

plotTC(parameterValues(9), parameterValues(1), parameterValues(2), parameterValues(3), ...
    parameterValues(4), parameterValues(5), parameterValues(6), parameterValues(8), parameterValues(7), toggleValues, ax)




% =========================================================================
function text_CallbackTCNUM(source, eventdata, handles, ax)

tcnum = str2double(get(source,'String'));

set(handles(9,1), 'Value', tcnum)

defaultparamsBMOC = csvread('defaultparamsBMOC.csv');
params = defaultparamsBMOC(:,tcnum);
freqScaling = get(handles(20,1), 'Value');
threshOC    = get(handles(21,1), 'Value');
singleCritical = get(handles(22,1), 'Value');

if singleCritical
    if freqScaling
        params = defaultparamsBMOC(33:40,tcnum);
    else
        params = defaultparamsBMOC(41:48,tcnum);
    end
else
    
    if threshOC
        if freqScaling
            params = defaultparamsBMOC(1:8, tcnum);
        else
            params = defaultparamsBMOC(9:16, tcnum);
        end
    else
        if freqScaling
            params = defaultparamsBMOC(17:24, tcnum);
        else
            params = defaultparamsBMOC(25:32, tcnum);
        end
    end
    
end

set(handles(10,1),  'String', num2str(params(1),3));
set(handles(11,1),  'String', num2str(params(2),3));
set(handles(12,1),  'String', num2str(params(3),3));
set(handles(13,1),  'String', num2str(params(4),3));
set(handles(14,1),  'String', num2str(params(5),3));
set(handles(15,1),  'String', num2str(params(6),3));
set(handles(16,1),  'String', num2str(params(7),3));
set(handles(17,1),  'String', num2str(params(8),3));
set(handles(18,1),  'String', num2str(tcnum));

toggleValues    = get(handles(19:22,1), 'Value');
parameterValues = [params;tcnum];

clear defaultparamsBMOC

plotTC(parameterValues(9), parameterValues(1), parameterValues(2), parameterValues(3), ...
    parameterValues(4), parameterValues(5), parameterValues(6), parameterValues(8), parameterValues(7), toggleValues, ax)



% =========================================================================
% =========================================================================
function plotTC(tcnum, alphaoc, alphabm, beta1oc, delta, c21, c12, TCfreqshift, thresh, toggleValue, ax)

% Load the Empirical Tuning Curve data.
load('JorisTCdata.mat');

% Define Reference Values
pmin    = dB2Pa(  0);
pmax    = dB2Pa(120);
epsilon = 1/pmax^2;
Fref    = 0.00002;  % Reference pressure in Pa.

% Store TC
TCtips = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tcnum = 3;  % Choose a Tuning Curve.  1 - 11
tcnum = round(tcnum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


if toggleValue{4}  % If single critical model
    
    if toggleValue{2}  % If frequency scaling
        F = thresh/c21*sqrt((alphaoc+beta1oc*thresh^2)^2 + (omega/f + delta*thresh^2).^2);
    else  % Else not frequency scaling
        F = thresh/c21*sqrt((alphaoc+beta1oc*thresh^2)^2 + (omega + delta*thresh^2).^2);        
    end
    
else  % Else coupled oscillator model
    
    if toggleValue{3}  % If taking roc as threshold
        
        if toggleValue{2}  % If frequency scaling
            
            rbm = thresh/c21*sqrt((alphaoc + beta1oc*thresh^2)^2 + (omega/f + delta*thresh^2).^2);
            sinPsioc = thresh*(omega/f + delta*thresh^2)/(c21*rbm);
            cosPsioc = sqrt(1-sinPsioc.^2);
            F = sqrt((alphabm*rbm + c12*thresh*cosPsioc).^2 + (omega/f.*rbm + c12*thresh*sinPsioc).^2);    % Ji Chul's calculation for bidirectional coupling
            
        else  % Else not frequency scaling
            
            rbm = thresh/c21*sqrt((alphaoc + beta1oc*thresh^2)^2 + (omega + delta*thresh^2).^2);
            sinPsioc = thresh*(omega + delta*thresh^2)/(c21*rbm);
            cosPsioc = sqrt(1-sinPsioc.^2);
            F = sqrt((alphabm*rbm + c12*thresh*cosPsioc).^2 + (omega.*rbm + c12*thresh*sinPsioc).^2);    % Ji Chul's calculation for bidirectional coupling
            
        end
        
    else  % else taking rbm as threshold
        
        if toggleValue{2}  % If frequency scaling
            
            F = NaN(3,length(omega));
            for j = 1:length(omega)
                
                roc = sqrt(roots([beta1oc^2 + omega2(j), 2*(alphaoc*beta1oc + delta*omega(j)/f), alphaoc^2+omega2(j), -(c21^2*thresh^2)]));
                cosPsioc = -(alphaoc*roc+beta1oc*roc.^3)/(c21*thresh);
                sinPsioc = (omega(j)*roc/f+delta*roc.^3)/(c21*thresh);
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
                
                roc = sqrt(roots([beta1oc^2 + omegasquared(j), 2*(alphaoc*beta1oc + delta*omega(j)), alphaoc^2+omegasquared(j), -(c21^2*thresh^2)]));
                cosPsioc = -(alphaoc*roc+beta1oc*roc.^3)/(c21*thresh);
                sinPsioc = (omega(j)*roc+delta*roc.^3)/(c21*thresh);
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

if toggleValue{4}  % Single critical model
    
    FStable = F;
    
else  % Coupled oscillator model
    
    if toggleValue{3}  % OC is threshold
        for j=1:length(omega)
            for i=1:size(F,1)
                if ~isnan(F(i,j))
                    if toggleValue{2}  % If frequency scaling
                        [~,r2,~,~,stable] = rStarCochlea(alphabm,alphaoc,beta1oc,delta,F(i,j),c21,c12,omega(j)/f,1);
                    else               % Else not frequency scaling
                        [~,r2,~,~,stable] = rStarCochlea(alphabm,alphaoc,beta1oc,delta,F(i,j),c21,c12,omega(j),1);
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
                    if toggleValue{2}  % If frequency scaling
                        [r1,~,~,~,stable] = rStarCochlea(alphabm,alphaoc,beta1oc,delta,F(i,j),c21,c12,omega(j)/f,1);
                    else               % Else not frequency scaling
                        [r1,~,~,~,stable] = rStarCochlea(alphabm,alphaoc,beta1oc,delta,F(i,j),c21,c12,omega(j),1);
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
if toggleValue{1}  % If middle ear filtering
    LstimStable   = 20*log10(FStable/Fref) - MEFdBgains;
    LstimUnstable = 20*log10(FUnstable/Fref) - repmat(MEFdBgains,3,1);
else               % Else not middle ear flitering
    LstimStable   = 20*log10(FStable/Fref);
    LstimUnstable = 20*log10(FUnstable/Fref);
end

semilogx(ax(1), f0, thresholds, 'b.-', 'Linewidth', 2);  hold(ax(1), 'on');
semilogx(ax(1), f0 + TCfreqshift, LstimStable,   'Color', [.2 .8 .2], 'LineWidth', 2, 'Marker', '.');
semilogx(ax(1), f0 + TCfreqshift, LstimUnstable, 'Color', [.8 .2 .8], 'LineWidth', 2, 'Marker', '.');

% semilogx(ax(1), f0, thresholds, 'b');  hold(ax(1), 'on');
% semilogx(ax(1), f0 + TCfreqshift, LstimStable,   'Color', [.2 .8 .2], 'LineWidth', .9);
% semilogx(ax(1), f0 + TCfreqshift, LstimUnstable, 'Color', [.8 .2 .8], 'LineWidth', .9);

xlim(ax(1), [44.4 30000]);
ylim(ax(1), [0 90]);
grid(ax(1), 'on');
ylabel(ax(1),'Threshold (dB SPL)');
xlabel(ax(1),'Forcing frequency (Hz)');
title(ax(1), ['TC ' num2str(tcnum) ' Fit  with  Tip = (' num2str(f) ', ' num2str(LstimStable(idxcf)) ')']);
set(ax(1), 'FontUnits','normalized');

hold(ax(1), 'off');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude Curves Method 2:  Compute roc's by specifying F's             %
%                             then use the roc's to compute rbm's         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = logspace(log10(Fref),log10(20),201);

rBM = [];
rOC = [];
F_rBM_at_rOCthreshold = [];

clear f0;
n = 0;

threef0 = [0.5*f f 2*f];

for f0 = threef0
    n = n+1;
    
    rstarbm = NaN(size(F));
    rstaroc = NaN(size(F));
    for nF = 1:length(F)
        if toggleValue{4}  % If single critical model
            if toggleValue{2}  % If frequency scaling
                temp = sqrt(roots([beta1oc^2+delta^2,...
                                   2*(alphaoc*beta1oc+delta*2*pi*(f-f0)/f),...
                                   alphaoc^2+(2*pi*(f-f0)/f)^2,...
                                   -F(nF)^2]));
                rstaroc(nF) = max(temp(~imag(temp)));
            else  % Else not frequency scaling
                temp = sqrt(roots([beta1oc^2+delta^2,...
                                   2*(alphaoc*beta1oc+delta*2*pi*(f-f0)),...
                                   alphaoc^2+(2*pi*(f-f0))^2,...
                                   -F(nF)^2]));
                rstaroc(nF) = max(temp(~imag(temp)));         
            end
        else  % Else coupled oscillator model
            if toggleValue{2}  % If frequency scaling
                [rbmn, rocn] = rStarCochlea(alphabm, alphaoc, beta1oc, delta, ...
                    F(nF), c21, c12, 2*pi*(f-f0)/f);
            else  % Else not frequency scaling
                [rbmn, rocn] = rStarCochlea(alphabm, alphaoc, beta1oc, delta, ...
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

if c21*c12   % both non-zero
    
    [spontBM, spontOC] = spontAmpCochlea(alphabm, alphaoc, beta1oc, delta, c21, c12);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Amplitude Curves %
%%%%%%%%%%%%%%%%%%%%%%%%%

F = Pa2dB(F); % For a log-log plot.
%F_rBM_at_rOCthreshold(:,1) = Pa2dB(F_rBM_at_rOCthreshold(:,1));

if toggleValue{4}  % If single critical oscillator

    semilogy(ax(2), F, rOC(1,:), 'g', F, rOC(2,:), 'b', F, rOC(3,:), 'r', [Fref 120], [thresh thresh], 'k', 'Linewidth', 2); %, [0 max(F)] , [rbm rbm], 'k');
    plot(ax(3),NaN);
    
else  % Else coupled oscillator model
    
    if toggleValue{3}
        
        if c21*c12
            semilogy(ax(2), F, rOC(1,:), 'g', F, rOC(2,:), 'b', F, rOC(3,:), 'r', [Fref 120], [thresh thresh], 'k', [Fref 120] , [spontOC spontOC], 'm', 'Linewidth', 2); %, [0 max(F)] , [rbm rbm], 'k');
            semilogy(ax(3), F, rBM(1,:), 'g', F, rBM(2,:), 'b', F, rBM(3,:), 'r', [Fref 120], [spontBM spontBM], 'm', 'Linewidth', 2);
        else
            semilogy(ax(2), F, rOC(1,:), 'g', F, rOC(2,:), 'b', F, rOC(3,:), 'r', [Fref 120], [thresh thresh], 'k', 'Linewidth', 2); %, [0 max(F)] , [rbm rbm], 'k');
            semilogy(ax(3), F, rBM(1,:), 'g', F, rBM(2,:), 'b', F, rBM(3,:), 'r', 'Linewidth', 2);
        end
        
    else
        
        if c21*c12
            semilogy(ax(2), F, rOC(1,:), 'g', F, rOC(2,:), 'b', F, rOC(3,:), 'r', [Fref 120], [spontOC spontOC], 'm', 'Linewidth', 2); %, [0 max(F)] , [rbm rbm], 'k');
            semilogy(ax(3), F, rBM(1,:), 'g', F, rBM(2,:), 'b', F, rBM(3,:), 'r', [Fref 120], [thresh thresh], 'k', [Fref 120], [spontBM spontBM], 'm', 'Linewidth', 2);
        else
            semilogy(ax(2), F, rOC(1,:), 'g', F, rOC(2,:), 'b', F, rOC(3,:), 'r', 'Linewidth', 2); %, [0 max(F)] , [rbm rbm], 'k');
            semilogy(ax(3), F, rBM(1,:), 'g', F, rBM(2,:), 'b', F, rBM(3,:), 'r', [Fref 120], [thresh thresh], 'k', 'Linewidth', 2);
        end
        
    end
    
end

xlabel(ax(2), 'Forcing F (dB SPL)');
ylabel(ax(2), 'r* of OC oscillator');
title(ax(2), {['Steady State Amplitude Curve for OC Oscillator.'], ['Blue @ f0 = CF = f, Green @ f0 = f/2,  Red @ f0 = 2f,  Black = Threshold']});
set(ax(2), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
set(ax(2), 'FontUnits','normalized');

grid(ax(3), 'on');

xlabel(ax(3), 'Forcing F (dB SPL)');
ylabel(ax(3), 'r* of BM oscillator');
title(ax(3), {['Steady State Amplitude Curve for BM Oscillator.'], ['Blue @ f0 = CF = f, Green @ f0 = f/2,  Red @ f0 = 2f']});
set(ax(3), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
set(ax(3), 'FontUnits','normalized');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OC & BM Compression Curves %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f0sOrig = logspace(log10(20),log10(30000),201)';
temp    = f0sOrig - f;
[~,ind] = min(abs(temp));
f0s     = f0sOrig - temp(ind);
F       = dB2Pa([0 20 40 60 80 100 120]);
if toggleValue{2}  % If frequency scaling
    W = 2*pi*(f-f0s)/f;
else               % Else not frequency scaling
    W = 2*pi*(f-f0s);
end
BMcomp = zeros(length(W),length(F));
OCcomp = zeros(length(W),length(F));

for nF = 1:length(F)
    for nW = 1:length(W)
        if toggleValue{4}  % If single critical oscillator
            temp = sqrt(roots([beta1oc^2+delta^2,...
                              2*(alphaoc*beta1oc+delta*W(nW)),...
                              alphaoc^2+W(nW)^2,...
                              -F(nF)^2]));
            OCcomp(nW,nF) = max(temp(~imag(temp)));
        else  % Else coupled oscillator model
            
            [rbmn, rocn] = rStarCochlea(alphabm, alphaoc, beta1oc, delta, ...
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

if toggleValue{4}  % If single critical oscillator
    
    loglog(ax(5), f0s(f0s>0), OCcomp(f0s>0,:), 'Linewidth', 2); hold(ax(5), 'on');
    loglog(ax(5), [min(f0sOrig) max(f0sOrig)], [thresh thresh], 'k', 'Linewidth', 2); hold(ax(5), 'off');
    loglog(ax(6), NaN);
    
else
    
    if toggleValue{3}  % If using OC as threshold
        
        if c21*c12
            loglog(ax(6), f0s(f0s>0), BMcomp(f0s>0,:), 'Linewidth', 2);hold(ax(6), 'on');
            loglog(ax(6), [min(f0sOrig) max(f0sOrig)], [spontBM spontBM], 'm', 'Linewidth', 2);hold(ax(6), 'off');
            loglog(ax(5), f0s(f0s>0), OCcomp(f0s>0,:), 'Linewidth', 2); hold(ax(5), 'on');
            loglog(ax(5), [min(f0sOrig) max(f0sOrig)], [thresh thresh], 'k', 'Linewidth', 2);
            loglog(ax(5), [min(f0sOrig) max(f0sOrig)], [spontOC spontOC], 'm', 'Linewidth', 2); hold(ax(5), 'off');
        else
            loglog(ax(6), f0s(f0s>0), BMcomp(f0s>0,:), 'Linewidth', 2);
            loglog(ax(5), f0s(f0s>0), OCcomp(f0s>0,:), 'Linewidth', 2); hold(ax(5), 'on');
            loglog(ax(5), [min(f0sOrig) max(f0sOrig)], [thresh thresh], 'k', 'Linewidth', 2); hold(ax(5), 'off');
        end
        
    else  % Else using BM as threshold
        
        if c21*c12
            loglog(ax(6), f0s(f0s>0), BMcomp(f0s>0,:), 'Linewidth', 2);hold(ax(6), 'on');
            loglog(ax(6), [min(f0sOrig) max(f0sOrig)], [thresh thresh], 'k', 'Linewidth', 2);
            loglog(ax(6), [min(f0sOrig) max(f0sOrig)], [spontBM spontBM], 'm', 'Linewidth', 2);hold(ax(6), 'off');
            loglog(ax(5), f0s(f0s>0), OCcomp(f0s>0,:), 'Linewidth', 2); hold(ax(5), 'on');
            loglog(ax(5), [min(f0sOrig) max(f0sOrig)], [spontOC spontOC], 'm', 'Linewidth', 2); hold(ax(5), 'off');
        else
            loglog(ax(6), f0s(f0s>0), BMcomp(f0s>0,:), 'Linewidth', 2);hold(ax(6), 'on');
            loglog(ax(6), [min(f0sOrig) max(f0sOrig)], [thresh thresh], 'k', 'Linewidth', 2); hold(ax(6), 'off');
            loglog(ax(5), f0s(f0s>0), OCcomp(f0s>0,:), 'Linewidth', 2);
        end
        
    end
    
end


xlim(ax(6), [min(f0sOrig) max(f0sOrig)]);
set(ax(6), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
xlabel(ax(6), 'Forcing frequency (Hz)');
ylabel(ax(6), 'rBM^*');
title(ax(6), 'BM Compression Curves.  0 - 120 dB SPL in 20 dB steps.');
set(ax(6), 'FontUnits','normalized');

xlim(ax(5), [min(f0sOrig) max(f0sOrig)]);
ylim(ax(5), [thresh/10 2*max(rOC(2,:))] );
set(ax(5), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
xlabel(ax(5), 'Forcing frequency (Hz)');
ylabel(ax(5), 'rOC^*');
title(ax(5), 'OC Compression Curves:  Forcing is provided by the BM oscillator dependent on f0.');
set(ax(5), 'FontUnits','normalized');

