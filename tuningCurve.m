function [LstimStable, LstimUnstable] = tuningCurve(cf, f0, alphabm, alphaoc, beta1oc, delta1oc, c21, c12, thresh, toggleValues)
% toggleValues must be a cell array of zeros and/or ones as follows:
% toggleValues{1} == 0 for no middle ear filter, == 1 for middle ear filter
% toggleValues{2} == 0 for unscaled, == 1 for scaled
% toggleValues{3} == 0 for taking rbm as threshold, == 1 for taking roc
% toggleValues{4} == 0 for single critical, == 1 for coupled model


% Define Reference Value
Fref    = 0.00002;  % Reference pressure in Pa.

f        = cf;
fsquared = f^2;

omega        =  2*pi*(f0 - f);
omegasquared = (2*pi*(f0 - f)).^2;
omega2       = omegasquared/fsquared;

MEFdBgains   = MEFdBgain(f0);  % In dB

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
