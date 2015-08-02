function G = MEFdBgain(f, varargin)

%
% Middle Ear Filter Gain Function in dB SPL.
%
% f = Frequency (Hz) Row Vector
% animal is either Cat or Human where
% Cat = 1, Human = 2
%
% gain = The gain of the digital filter, defaults to 1.
%


if length(varargin) < 1
    animal = 2;  % Cat = 1, Human = 2
    gain   = 1;
elseif  length(varargin) == 1
    animal = varargin(1);
    gain   = 1;
else
    animal = varargin(1);
    gain   = varargin(2);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Middle-Ear Filter
% Using the Zilany et al. 2006 filter.

% Fs = 500e3;  % Sampling Frequency
% 
% % Specify  second-order section coefficients of the
% % digital filter.
% sos = [1  1 0 1 -0.9986 0; ... 
%        1 -1.9998 0.9998 1 -1.9777 0.9781; ... 
%        1 -1.9943 0.9973 1 -1.9856 0.9892];
% 
% % Note: g = 0.0127 is the gain
% [B, A]    = sos2tf(sos, 0.0127);
% MEFdBgain = 20*log10(abs(freqz(B, A, 1:40000, Fs)));
% clear A B;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Middle-Ear Filter

cat_or_human = animal; % Choose 1 for cat and 2 for human

% Sampling Frequency (this code should work/stable for Fs > 70e3)
Fs = 100e3;

tdres = 1/Fs;
fp    = 1e3;  % prewarping frequency 1 kHz 
TWOPI = 2*pi;
C     = TWOPI*fp/tan((TWOPI/2)*fp*tdres);

if cat_or_human == 1  % cat
    
    %%% cat middle ear function %%%%%%%%%%%%
    
    m11 = C/(C + 693.48); 
    m12 = (693.48 - C)/C; 
    m21 = 1/(power(C,2) + 11053*C + 1.163e8);
    m22 = -2*power(C,2) + 2.326e8;
    m23 = power(C,2) - 11053*C + 1.163e8;
    m24 = power(C,2) + 1356.3*C + 7.4417e8;
    m25 = -2*power(C,2) + 14.8834e8;
    m26 = power(C,2) - 1356.3*C + 7.4417e8; 
    m31 = 1/(power(C,2) + 4620*C + 909059944);
    m32 = -2*power(C,2) + 2*909059944; 
    m33 = power(C,2) - 4620*C + 909059944; 
    m34 = 5.7585e5*C + 7.1665e7; 
    m35 = 14.333e7; 
    m36 = 7.1665e7 - 5.7585e5*C;

    sos = [1 -1 0 1 m11*m12 0;...
           m21*m24 m21*m25 m21*m26 1 m21*m22 m21*m23;...
           m31*m34 m31*m35 m31*m36 1 m31*m32 m31*m33];
    %%% end of cat middle ear function

elseif cat_or_human == 2 % human

    %%% Human middle ear function %%%%%%%%%%%%%%%%%%%%

    m11 = 1/(power(C,2) + 5.9761e+003*C + 2.5255e+007);
    m12 = (-2*power(C,2) + 2*2.5255e+007);
    m13 = (power(C,2) - 5.9761e+003*C + 2.5255e+007);
    m14 = (power(C,2) + 5.6665e+003*C);
    m15 = -2*power(C,2);
    m16 = (power(C,2) - 5.6665e+003*C);
    
    m21 = 1/(power(C,2) + 6.4255e+003*C + 1.3975e+008);
    m22 = (-2*power(C,2) + 2*1.3975e+008);
    m23 = (power(C,2) - 6.4255e+003*C + 1.3975e+008);
    m24 = (power(C,2) + 5.8934e+003*C + 1.7926e+008);
    m25 = (-2*power(C,2) + 2*1.7926e+008);
    m26 = (power(C,2)-5.8934e+003*C+1.7926e+008);
    
    m31 = 1/(power(C,2) + 2.4891e+004*C+1.2700e+009);
    m32 = (-2*power(C,2) + 2*1.2700e+009);
    m33 = (power(C,2) - 2.4891e+004*C + 1.2700e+009);
    m34 = (3.1137e+003*C + 6.9768e+008);
    m35 = 2*6.9768e+008;
    m36 = (-3.1137e+003*C + 6.9768e+008);    
    
    sos = [m11*m14 m11*m15 m11*m16 1 m11*m12 m11*m13;...
           m21*m24 m21*m25 m21*m26 1 m21*m22 m21*m23;...
           m31*m34 m31*m35 m31*m36 1 m31*m32 m31*m33];

    %%% end of Human middle ear function %%%%%%%%%%%%%%
end

% Gain of 1 is typically used
[B, A]    = sos2tf(sos, gain);

%MEFdBgain = 20*log10(abs(freqz(B, A, 1:40000, Fs)));
G = 20*log10(abs(freqz(B, A, f, Fs)));




% f = 20:20000;
% cat_or_human = 2;  gain = 1;
% 
% 
% figure(1);
% subplot(2,1,1);
% semilogx(f, G, 'b');  grid on;
% %set(gca, 'xlim', [0 20000]);
% title('Gain, i.e., Frequency Response/Spectrum of Middle Ear Filter (MEF)');
% xlabel('Frequency (Hz)');
% ylabel('Gain (dB)');
% xlim([20 20000]);
% 
% subplot(2,1,2);
% semilogx(f, angle(freqz(B, A, f, Fs)), 'r');  grid on;
% title('Phase Response of Middle Ear Filter (MEF)');
% xlabel('Frequency (Hz)')
% ylabel('Phase (radians)');
% xlim([20 20000]);
% ylim([-pi pi]);

