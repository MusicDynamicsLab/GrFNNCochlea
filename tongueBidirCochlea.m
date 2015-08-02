%%% Arnold tongue plot for the cochlea paper
applyMEF = 1;
resolution = 201; % resolution of plot (for fine resolution, use 2001)

%fz = 152; % osc frequency
%fz = 1850;
fz = 7262;

abm = -412; % BM parameters
aoc = 0; b1oc = -40816; d1oc = 0; % OC parameters
c21 = 197960*fz; % BM-to-OC (afferent) coupling strength

rocThr = .1; % threshold roc
c12 = abm*(aoc+b1oc*rocThr^2)/c21*.99; % OC-to-BM (efferent) coupling strength

f0 = fz*2.^linspace(-3,3,resolution); % input frequency
W = (fz-f0)*2*pi; % radian frequency difference

FdBSPL = linspace(0,90,resolution); % forcing in dB SPL
MEFdBgains = MEFdBgain(f0); % middle ear filter gain

STABTYPE = zeros(length(FdBSPL),length(W));
RSTAROC = NaN(length(FdBSPL),length(W));

tic
for nf = 1:length(FdBSPL)
  for nw = 1:length(W)
    if applyMEF
      F = dB2Pa(FdBSPL(nf)+MEFdBgains(nw)); % forcing in Pa
    else
      F = dB2Pa(FdBSPL(nf)); % forcing in Pa
    end
    [rStarBM,rStarOC,psiStarBM,psiStarOC,stab,type] ...
      = rStarCochlea(abm,aoc,b1oc,d1oc,F,c21,c12,W(nw),1);
    if sum(stab) < 2
      STABTYPE(nf,nw) = max(type); % show stable type if any
    elseif sum(stab) == 2
      STABTYPE(nf,nw) = 5;
    else
      error('Three or more stable fixed points?')
    end
    RSTAROC(nf,nw) = max(rStarOC);
  end
end
toc

%% contour plot
colorPlot = 1; % 1 for color plot, 0 for gray plot

uniqueTypes = unique(STABTYPE);

figure
set(gcf,'Position',[450 320 560 420])

contourf(f0,FdBSPL,STABTYPE,uniqueTypes-.5)
set(gca,'XScale','log')
xlabel('Forcing frequency (Hz)','FontSize',12)
ylabel('Forcing amplitude (dB SPL)','FontSize',12)
grid on

if colorPlot
  cmap = [jet(5);.5 0 .5]; % stability color map
  colormap(cmap(uniqueTypes+1,:))
  caxis([min(uniqueTypes) max(uniqueTypes)])
else
  colormap(flipud(gray))
  caxis([0 7])
end

if rocThr % if threshold roc is set
  hold on
  [~,h] = contour(f0,FdBSPL,RSTAROC,rocThr,'r'); % draw contour
  set(h,'ShowText','on') % annotate with threshold value
  hold off
end
