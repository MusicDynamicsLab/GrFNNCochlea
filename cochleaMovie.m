% function cochleaMovie(n, s, 'fn', frate)
n2 = M.n{2};

figsize = [2 52 1440 810];
moviename = 'cochleaMEF1000-080.avi'

%% Set up global display parameters
mx = max(max(abs(n2.Z')))
% mx = .0147;
% mx = .14;

fig = figure(2001);
set(fig, 'Position',  figsize)

%% Signal Axis
as = subplot('Position', [0.1 0.7 0.75 0.2]);
set(gca, 'FontSize', 14);
ylabel('Pressure (Pa)');
xlabel('Time');

plot(s.t, real(s.x), 'k'); axis tight;
grid on;
yl = get(gca, 'YLim');
hold on; 
timeline = plot(s.t(1)*[1 1], yl, 'r');
hold off;
set(timeline, 'EraseMode', 'xor')

%% Text Display
at = subplot('Position', [0.875 0.65 0.25 0.25]);
set(gca, 'FontSize', 14, 'Visible', 'off')
xlim([ 0 1]); ylim([-1 1])

% if s.analytic == 0
    timetext = text(0, 1, sprintf('t = %6.3f', 0), 'FontSize', 14)
    % fctext   = text(0, .6, sprintf('f_c = %6.3f', 0), 'FontSize', 14)
    % dBctext  = text(0, .2, sprintf('dB_c = %6.3f', 0), 'FontSize', 14)
%else
%end
set(timetext, 'EraseMode', 'xor')

%% Cochlea Display
ac = subplot('Position', [0.1 0.135 0.8 0.5]);
set(gca, 'FontSize', 14)
xlabel('Frequency (Hz)');
ylabel('Displacement');

% envpos =  semilogx(n2.f,  abs(n2.Z(:,1)), 'r-', 'LineWidth', 2);
% hold on
% envneg =  semilogx(n2.f, -abs(n2.Z(:,1)), 'r-', 'LineWidth', 2);
cochleaplot =  semilogx(n2.f, real(n2.Z(:,1)), '-', 'LineWidth', 2);
% hold off
set(gca, 'XLim', [30 10000])                    % temporary change
set(gca, 'XDir', 'reverse')
set(gca, 'YLim', [-1.5 1.5]); % 1.20*[-mx, mx])     % temporary change
grid on;

%% Set up movie
writerObj = VideoWriter(moviename, 'Uncompressed AVI');
writerObj.FrameRate = 50;

set(fig, 'Position',  figsize)
% return

open(writerObj);

%% Loop through time to create movie
for k = 1:10:length(n2.t)
    
    % Signal
    set(timeline, 'XData', s.t(k)*[1 1])
    drawnow

    % Text
    set(timetext, 'String', sprintf('t = %6.3f', s.t(k)))

    % Cochlea
    set(cochleaplot, 'YData', real(n2.Z(:,k)))
%     set(envpos, 'YData',  abs(n2.Z(:,k)))
%     set(envneg, 'YData', -abs(n2.Z(:,k)))

    % Record frames
    drawnow
    frame = getframe(fig);
    writeVideo(writerObj,frame);

end

close(writerObj);
%
