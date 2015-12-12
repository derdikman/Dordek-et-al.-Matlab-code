function [f_dominant, fSize_dominant, periodicity]=dominantFreq(x,numOfElements)

% Fs = 1/mean(abs(diff(x)));%1000; % sampling frequency 1 kHz
% t = 0 : 1/Fs : 0.296; % time scale
% f = 200; % Hz, embedded dominant frequency
% x = cos(2*pi*f*t) + randn(size(t)); % time series
%plot(t,x), axis('tight'), grid('on'), title('Time series'), figure
%nfft = 512; % next larger power of 2
%nfft = 2^ceil(log2(length(x)));% next larger power of 2
y = fft(x-mean(x)); % Fast Fourier Transform + extract the mean to get rid of the DC
y = abs(y.^2); % raw power spectrum density
y = y(1:1+length(x)/2); % half-spectrum %DO NOT START FROM ZERO ! maybe extract mean?
[sortedValues,sortIndex] = sort(y,'descend');
maxIndex = sortIndex(1:numOfElements);
maxSize = sortedValues(1);
periodicity = maxSize/mean(sortedValues);
%[~,k] = max(y); % find maximum
%f_scale = (0:length(x)/2)* Fs/length(x); % frequency scale
%plot(f_scale, y),axis('tight'),grid('on'),title('Dominant Frequency')
f_dominant= (sortIndex(1:numOfElements)); % dominant frequency estimate
fSize_dominant  = sortedValues(1:numOfElements); %size of largest freq elements
%fprintf('Dominant freq.: true %f Hz, estimated %f Hz\n', f, f_est)
%fprintf('Frequency step (resolution) = %f Hz\n', f_scale(2))
end

% n = 0:199;
% x = cos(0.257*pi*n) + sin(0.2*pi*n) + 0.01*randn(size(n));
% [S,f] = pmusic(x,4)   ;   % Set p to 4 because there are two real inputs
% figure; plot(f*pi/max(f*pi),20*log10(S)), grid on

