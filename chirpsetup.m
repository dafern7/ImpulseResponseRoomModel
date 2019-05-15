clear all

fs = 44100;
Tp = 0.5;
N = fs*Tp;
fmax = 20000;
fmin = 0.1;
t = linspace(0,Tp,N);
avgnum = 30;
%3ticks

%linear and logarithmic chirps
phi = 2*pi*((1/2)*((fmax-fmin)/Tp)*t.^2 + fmin*t);
x = sin(phi);
%figure(1)
%256 averages, overlap of 250 samples for each for higher frequency
%resolution
%spectrogram(x,256,250,256,fs,'yaxis')
%soundsc(x)
%figure(1)
%plot(tlin,x)

philog = 2*pi*((fmin*((fmax/fmin).^(t/Tp)))*(Tp/log(fmax/fmin)));
y = sin(philog);
y = [y zeros(1,length(y))];
y = repmat(y,1,avgnum);
%figure(2)
%too lazy to use the one i made
%spectrogram(y,256,250,256,fs,'yaxis')
%soundsc(y,fs)
%figure(2)
%plot(tlog,y)

%backgroundnoise measurement
z = zeros(1,400000);

%testing stuff
s = daq.createSession('ni')
d = daq.getDevices


addAnalogOutputChannel(s,'Dev1','ao0','Voltage')
queueOutputData(s,y')
addAnalogInputChannel(s,'Dev1','ai0','IEPE');
s.Rate = 44100;
[data,timestamps] = startForeground(s);

audiowrite('loudness5location1dur0.5sec.wav',data,fs)
t1 = linspace(0,Tp*avgnum*2,N*avgnum*2);
plot(t1,data)


%next steps
%so set up repeating chirps, average the measured response PSD (Gxx) every Tp. Can
%reduce Tp if necessary. (10s is pretty long).

%H(f) = Y(f)/X(f) -> fft(y(t))/fft(x(t)), y is received response, x is
%original
%imp response = h(t) -> ifft(H(f))

%filter entire signal at each octaveband, cut in time domain.
%bandpass(x, [flow fhigh], fs)  
%flow,fhigh dependent on the octave band center frequency.
