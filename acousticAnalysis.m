[y,Fs] = audioread('loudness5location3dur1sec.wav');

%y = y(1000:end);
avgnum = 30;
N = Fs;
N1 = length(y);
fs = Fs;
fsamp = Fs;
dt = 1/fsamp;
dt1 = 1/fs;
df1 = 1/N/dt1;
f = (0:N1/2-1)*(1/N1/dt);
f1 = (-fs/2:(fs/2)-1)*df1;
Tp = 1;
N = fs*Tp;
fmax = 20000;
fmin = 0.1;
t = linspace(0,1,Fs);
philog = 2*pi*((fmin*((fmax/fmin).^(t/Tp)))*(Tp/log(fmax/fmin)));
x = sin(philog);
fftrealx = fft(x).*(1/fs);
fftrealx = fftshift(fftrealx);
[fftrealy,~,~,~,~,~,~] = lin_avg(y',fs,avgnum,fsamp);
impresp = (ifft(fftrealy./fftrealx));

%shift for each location, since every location has a different silence
%duration before the impulse response begins
silence_dur = 0.028;
shift = silence_dur*length(impresp);

impresps = impresp(shift:end);
early50 = sum(impresps(1:length(impresps)*0.05).^2);
late50 = sum(impresps(length(impresps)*0.05:end).^2);
C50_total = 10*log10(early50/late50)
early80 = sum(impresps(1:length(impresps)*0.08).^2);
late80 = sum(impresps(length(impresps)*0.08:end).^2);
C80_total = 10*log10(early80/late80)

figure(1)
plot(t,impresp);
title('General Impulse Response (Location 3)')
xlabel('Time (s)')
ylabel('Amplitude (WU)')


%big boi butter (filter) bank
%tested around with the order of the filter - higher makes the cutoff
%frequency better defined (makes the filter steeper), but higher orders did
%not work well with this dataset. n=3 showed the impulse response best.
[b1,a1] = butter(3,[88 177]/(Fs/2),'bandpass'); %125 Hz
[b2,a2] = butter(3,[177 355]/(Fs/2),'bandpass');%250 Hz
[b3,a3] = butter(3,[355 710]/(Fs/2),'bandpass');%500 Hz
[b4,a4] = butter(3,[710 1420]/(Fs/2),'bandpass');%1 kHz
[b5,a5] = butter(3,[1420 2840]/(Fs/2),'bandpass');%2 kHz
[b6,a6] = butter(3,[2840 5680]/(Fs/2),'bandpass');%4 kHz
%but do the filters look good in the frequency domain?!?! Need big boi
%filter freq response plots

figure(2)
freqz(b1,a1,1024,fs)
hold on
freqz(b2,a2,1024,fs)
hold on
freqz(b3,a3,1024,fs)
hold on
freqz(b4,a4,1024,fs)
hold on
freqz(b5,a5,1024,fs)
hold on
freqz(b6,a6,1024,fs)

ax = findall(gcf, 'Type', 'axes');
set(ax, 'XScale', 'log');
title('Frequency Response of each Filter')
%on the logarithmic scale, each filter has around an equal size passband
%region, indicating the octave bands. The transition band for each is
%rather spread, which is a consequence of the low order of the filters, but
%having a high order filter yields undesirable results, perhaps due to
%stability?

%overall, acceptable filters.

%big boi filtering
filteredy125 = filter(b1,a1,y);
filteredy250 = filter(b2,a2,y);
filteredy500 = filter(b3,a3,y);
filteredy1k = filter(b4,a4,y);
filteredy2k = filter(b5,a5,y);
filteredy4k = filter(b6,a6,y);
%big boi ffts for the filtered bois
[fftrealy125,~,~,~,~,~,~] = lin_avg(filteredy125',fs,avgnum,fsamp);
[fftrealy250,~,~,~,~,~,~] = lin_avg(filteredy250',fs,avgnum,fsamp);
[fftrealy500,~,~,~,~,~,~] = lin_avg(filteredy500',fs,avgnum,fsamp);
[fftrealy1k,~,~,~,~,~,~] = lin_avg(filteredy1k',fs,avgnum,fsamp);
[fftrealy2k,~,~,~,~,~,~] = lin_avg(filteredy2k',fs,avgnum,fsamp);
[fftrealy4k,~,~,~,~,~,~] = lin_avg(filteredy4k',fs,avgnum,fsamp);
%big boi impulse response calcs. I fftshifted everything before, so I did not need
%to ifftshift back.
impresp125 = (ifft(fftrealy125./fftrealx));
impresp250 = (ifft(fftrealy250./fftrealx));
impresp500 = (ifft(fftrealy500./fftrealx));
impresp1k = (ifft(fftrealy1k./fftrealx));
impresp2k = (ifft(fftrealy2k./fftrealx));
impresp4k = (ifft(fftrealy4k./fftrealx));


impresp125s = impresp125(shift:end);
impresp250s = impresp250(shift:end);
impresp500s = impresp500(shift:end);
impresp1ks = impresp1k(shift:end);
impresp2ks = impresp2k(shift:end);
impresp4ks = impresp4k(shift:end);

%for clarity index calculations, abuse the fact that the impulse response
%is 1 second and the the first 50 milliseconds is just a fraction of it
early50_125 = sum(impresp125s(1:length(impresp125s)*0.05).^2);
late50_125 = sum(impresp125s(length(impresp125s)*0.05:end).^2);
C50_125 = 10*log10(early50_125/late50_125);
early80_125 = sum(impresp125s(1:length(impresp125s)*0.08).^2);
late80_125 = sum(impresp125s(length(impresp125s)*0.08:end).^2);
C80_125 = 10*log10(early80_125/late80_125);

early50_250 = sum(impresp250s(1:length(impresp250s)*0.05).^2);
late50_250 = sum(impresp250s(length(impresp250s)*0.05:end).^2);
C50_250 = 10*log10(early50_250/late50_250);
early80_250 = sum(impresp250s(1:length(impresp250s)*0.08).^2);
late80_250 = sum(impresp250s(length(impresp250s)*0.08:end).^2);
C80_250 = 10*log10(early80_250/late80_250);

early50_500 = sum(impresp500s(1:length(impresp500s)*0.05).^2);
late50_500 = sum(impresp500s(length(impresp500s)*0.05:end).^2);
C50_500 = 10*log10(early50_500/late50_500);
early80_500 = sum(impresp500s(1:length(impresp500s)*0.08).^2);
late80_500 = sum(impresp500s(length(impresp500s)*0.08:end).^2);
C80_500 = 10*log10(early80_500/late80_500);

early50_1k = sum(impresp1ks(1:length(impresp1ks)*0.05).^2);
late50_1k = sum(impresp1ks(length(impresp1ks)*0.05:end).^2);
C50_1k = 10*log10(early50_1k/late50_1k);
early80_1k = sum(impresp1ks(1:length(impresp1ks)*0.08).^2);
late80_1k = sum(impresp1ks(length(impresp1ks)*0.08:end).^2);
C80_1k = 10*log10(early80_1k/late80_1k);

early50_2k = sum(impresp2ks(1:length(impresp2ks)*0.05).^2);
late50_2k = sum(impresp2ks(length(impresp2ks)*0.05:end).^2);
C50_2k = 10*log10(early50_2k/late50_2k);
early80_2k = sum(impresp2ks(1:length(impresp2ks)*0.08).^2);
late80_2k = sum(impresp2ks(length(impresp2ks)*0.08:end).^2);
C80_2k = 10*log10(early80_2k/late80_2k);

early50_4k = sum(impresp4ks(1:length(impresp4ks)*0.05).^2);
late50_4k = sum(impresp4ks(length(impresp4ks)*0.05:end).^2);
C50_4k = 10*log10(early50_4k/late50_4k);
early80_4k = sum(impresp4ks(1:length(impresp4ks)*0.08).^2);
late80_4k = sum(impresp4ks(length(impresp4ks)*0.08:end).^2);
C80_4k = 10*log10(early80_4k/late80_4k);

C50 = [C50_125 C50_250 C50_500 C50_1k C50_2k C50_4k]
C80 = [C80_125 C80_250 C80_500 C80_1k C80_2k C80_4k]
bands = [125 250 500 1000 2000 4000];


%big boi plot
figure(3)
plot(t,(impresp125))
hold on
plot(t,(impresp250))
hold on
plot(t,(impresp500))
hold on
plot(t,(impresp1k))
hold on
plot(t,(impresp2k))
hold on
plot(t,(impresp4k))
legend('125Hz','250Hz','500Hz','1kHz','2kHz','4kHz')
title('Impulse Responses of Octave Bands(Location 3)')
xlabel('Time (s)')
ylabel('Amplitude (WU)')

%attempted exponential filter stuff :(
%frac = 44100;
%deltaT = floor(Fs/frac);
%Tc = 0.2;
%bexpfilt = [0,1];
%aexpfilt = [1,Tc];
%impresp250expfiltered = filter(b,a,impresp250);
%figure(2)
%plot(t,impresp250expfiltered)

%getting the schroeder curves for the reverb time calculations. 
%sadly does not reach the noise floor because we only did 1s
impresp125 = fliplr(impresp125);
impresp125 = cumsum(abs(impresp125));
impresp125 = fliplr(impresp125);

impresp250 = fliplr(impresp250);
impresp250 = cumsum(abs(impresp250));
impresp250 = fliplr(impresp250);

impresp500 = fliplr(impresp500);
impresp500 = cumsum(abs(impresp500));
impresp500 = fliplr(impresp500);

impresp1k = fliplr(impresp1k);
impresp1k = cumsum(abs(impresp1k));
impresp1k = fliplr(impresp1k);

impresp2k = fliplr(impresp2k);
impresp2k = cumsum(abs(impresp2k));
impresp2k = fliplr(impresp2k);

impresp4k = fliplr(impresp4k);
impresp4k = cumsum(abs(impresp4k));
impresp4k = fliplr(impresp4k);

figure(4)
plot(t,10*log((impresp125)));
hold on
plot(t,10*log(impresp250))
hold on
plot(t,10*log(impresp500))
hold on
plot(t,10*log(impresp1k))
hold on
plot(t,10*log(impresp2k))
hold on
plot(t,10*log(impresp4k))
legend('125Hz','250Hz','500Hz','1kHz','2kHz','4kHz')
title('Schroeder Curves for each octave band (Location 3)')
xlabel('Time (s)')
ylabel('Amplitude (dB)')




%y - input signal
%fs - how many samples should be taken for each record
    %note: T = fs/fsamp (where T is the duration of each record)
%numAvg - number of averages/records to take
%fsamp - sampling frequency for the entire input signal
%NOTE: since fs is the number of samples and is mathematically found from
%duration and sampling frequency, using my PSD function which takes fsamp
%and T, it is equivalent to put fs (number of samples) for the fsamp input
%and 1 as T, since it creates the same length vector.
%the PSD function can be simplified further because of this
function [fftreal,time_AVG,rmsmean,fullGxx,vect_avg,time_avg,N] = lin_avg(y,fs,numAvg,fsamp)
%creating the length of the vector given sampling frequency and time
T = fs/fsamp;
N = fs*T;
N1 = fsamp*(length(y)/fsamp);
[fullGxx,~,~] = psd(y,fsamp,length(y)/fsamp);
GXX = [];
avg_fft = [];
Y = [];
%storing the averages and input records in 2D matrices
for i = 1:numAvg
    [Gxx,Sxx,fftreal] = psd(y((i-1)*fs+1:i*fs),fs,1);
    avg_fft = [avg_fft;fftreal]; %FFT of records
    GXX = [GXX;Gxx]; %for RMS average
    Y = [Y;y((i-1)*fs+1:i*fs)]; %Input records
end
%averaging for vector average and RMS
rmsmean = mean(GXX);
avg_fft = mean(avg_fft);
avg_Sxx = (1)*(abs(avg_fft)).^2; %T = 1 (one unit of time), not necessarily a second
vect_avg = 2*avg_Sxx(fs/2+1:fs);
%averaging for time average
time_AVG = mean(Y);
[time_avg,Sxx,fftreal] = psd(time_AVG,fs,1);

%generating the frequency(Hz) vectors to plot against
dt = 1/fsamp;
dt1 = 1/fs;
df1 = 1/N/dt1;
f = (0:N1/2-1)*(1/N1/dt);
f1 = (0:(fs/2)-1)*df1;

%{
%plotting a zoomed in (y axis) view of both RMS and Synchronous Averaging to see the
%overtone clarity
figure(1)
subplot(2,1,1)
plot(f1,vect_avg,'Linewidth',1)
hold on
plot(f1,time_avg,'--','Linewidth',2)
legend('Vector Averaging','Time Averaging')
title('Comparison of Overtones')
xlabel('Frequencies (Hz)')
ylabel('Amplitude(WU)')
axis([0 fsamp/2 0 2e-4]); 
    
subplot(2,1,2)
plot(f1,rmsmean,'Linewidth',1)
title('RMSmean Overtones')
xlabel('Frequencies (Hz)')
ylabel('Amplitude(WU)')
axis([0 fsamp/2 0 3e-4]); 

%plotting the general power spectrum. Not sure if this is right.
figure(2)
plot(f,fullGxx)
title('Gxx without averaging the signal')
xlabel('Frequencies (Hz)')
ylabel('Amplitude(WU)')

%plotting the full height of both synchronous method and RMS to see
%amplitude differences
figure(3)
semilogx(f1,10*log(vect_avg),'Linewidth',1)
hold on
semilogx(f1,10*log(time_avg),'--','Linewidth',2)
%hold on
%plot(f1,rmsmean)
legend('Vector Averaging','Time Averaging','RMS averaging')
title('Comparison of Synchronous Methods')
xlabel('Frequencies (Hz)')
ylabel('Amplitude(WU)')

%comparison with RMS average, used if there is no synchronization.
figure(5)
plot(f1,rmsmean)
title('RMS Averaging (no Synchronization)')
xlabel('Frequencies (Hz)')
ylabel('Amplitude(WU)')
%}
end

function [Gxx,Sxx,fftreal] = psd(x,fs,T)
N = T*fs;
fftreal = fft(x).*(1/fs); %scaling each fft value with the dt (total time/num samples)
fftreal = fftshift(fftreal);
Sxx = (1/T)*(abs(fftreal)).^2; %power spectral density from fft
Gxx = 2*Sxx(N/2+1:N); %one sided PSD (taking the right handed side of the Sxx)
end
