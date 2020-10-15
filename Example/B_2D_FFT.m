% 2-D Transform
% The 2-D Fourier transform is useful for processing 2-D signals and other 2-D data such as images.
% Create and plot 2-D data with repeated blocks.
clear all;
close all;

P = peaks(20);
X = repmat(P,[5 10]);
figure(1)
imagesc(X)

% TODO : Compute the 2-D Fourier transform of the data.  
% Shift the zero-frequency component to the center of the output, and 
% plot the resulting 100-by-200 matrix, which is the same size as X.


% In the case of Radar signal processing. Convert the signal in MxN matrix, where M is the size of Range FFT samples and N is the size of Doppler FFT samples:
figure(2)
signal  = reshape(X, [100, 200]);
imagesc(signal);

% Run the 2D FFT across both the dimensions.
figure(3)
signal_fft = fft2(signal, 100, 200);
imagesc(abs(signal_fft));

% Shifat zero-frequency terms to the center of the array
signal_fft = fftshift(signal_fft);

% Take the absolute value
signal_fft = abs(signal_fft);

% Here since it is a 2D output, it can be plotted as an image. Hence, we use the imagesc function
figure(5)
imagesc(signal_fft);