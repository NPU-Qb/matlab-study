clc;
dis=[]
for i=1:6
fs=4096;
x=dis(:,i);   
N=8192;
df=fs/N;
f=0:df:N*df-df;
FFT=fft(x);
margin(:,i)=abs(FFT)*2/N;
end
figure(1);
plot(f,margin(:,5));
xlabel('频率/Hz')
ylabel('x /m')
xlim([100 4000])
clear all