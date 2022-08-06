
clear all
%�������� �ɵĿ������㷨 �Ӿ���������
a=-0.01;b=0.1;c=0.5;%�������6.1������֤
n=1;acce(1)=2;ts2(1)=0;vol(1)=0.9;disp(1)=0;ts1(1)=0;dis(1)=0;ts3(1)=0;
vol(2)=1+(-0.1+2*0.002)*cos((-0.1+0.002)*0.002);
while ts2<20  %���ٶȼ��㣬����Ƶ��Ϊ1000
ts2(n+1)=ts2(n)+0.001;
n=n+1;
acce(n)=2*c*cos((a+c*ts2(n))*ts2(n))-(a+2*c*ts2(n))^2*sin((a+c*ts2(n))*ts2(n));
end
n=1;
while ts1<20  %λ�Ƽ��㣬����Ƶ��Ϊ100
ts1(n+1)=ts1(n)+0.01;
n=n+1;
disp(n)=sin((a+c*ts1(n))*ts1(n))+b*ts1(n);
end
ha=1/1000;
hd=1/100;

% Nacce=load('acc.txt');
% Ndisp=load('dis1.txt');
% dis=load('dis.txt');
% adis=dis(100:149,1);
% acce=Nacce(496:746,1);
% disp=Ndisp(100:149,1);
% 
% rmsA=randn(1,length(acce)); 
% rmsA=rmsA/std(rmsA); 
% rmsA=rmsA-mean(rmsA); 
% rmsa=sqrt(0.1); 
% rmsA=rmsA*rmsa;
% nacce=acce+rmsA';
% 
% rmsD=randn(1,length(disp)); 
% rmsD=rmsD/std(rmsD); 
% rmsD=rmsD-mean(rmsD); 
% rmsd=sqrt(0.1); 
% rmsD=rmsD*rmsd;
% ndisp=disp+rmsD';
% %rmsd=0.0000001;

%%����������ʱ����
rmsa=sum(acce.*acce);
rmsa=(rmsa/length(acce))^0.5;
nacce=acce+0.1*rmsa*randn(1,length(acce));%   nacce - noised acceleration measuremnts

rmsd=sum(disp.*disp);
rmsd=(rmsd/length(disp))^0.5;
ndisp=disp+0.1*rmsd*randn(1,length(disp));%   ndisp - noised displacment measuremnts

%%��ʼ�����趨
dt=ha; % ���ٶȲ���ʱ�䣨����Ƶ�ʵ�����
dd=hd;% λ�Ʋ���ʱ�䣨����Ƶ�ʵ�����
sampletimes=hd/ha;% �������ʱ��
mk=0;% ��ʼʱ��
x0=0; % λ�Ƴ�ʼԤ��ֵ
v0=0;% �ٶȳ�ʼԤ��ֵ

%%���·��̡���������ģ���й� ��ʱ����
A = [1 dt; 0 1] ; 
B = [dt^2/2; dt]; 
C = [1 0];       
%%��Ҫ��������***
Q_loc_estimate = [];% smooth of displacment
vel_smooth = [];% smooth of velocity
Q_estimate=[x0;v0];% x_estimate of initial location estimation
QuailAccel_noise_mag =var(nacce-acce);% process noise
NinjaVision_noise_mag =var(ndisp-disp); % measurement noise
Ex = QuailAccel_noise_mag * [dt^3/3 dt^2/2; dt^2/2 dt];
Ez = NinjaVision_noise_mag ; 
Pcorr = [1 0;1 0];

loc_meas = ndisp;% displacment measurement

loc_estimate = []; % displacment estimate
vel_estimate = [];  % velocity estimate  
loc_estimate_pred = [];
vel_estimate_pred = [];
Pcorr_estimate = [];
Ppred_estimate = [];
Q = [x0; v0]; 

for t = 1:1:length(nacce)-1
   
   u = nacce(t);
   Q_estimate = A * Q_estimate + B * u;
   Ppred = A * Pcorr * A' + Ex;
   loc_estimate_pred = [loc_estimate_pred; Q_estimate(1)];
   vel_estimate_pred = [vel_estimate_pred; Q_estimate(2)];
   Ppred_estimate = [Ppred_estimate;Ppred];
   Pcorr_estimate = [Pcorr_estimate;Ppred]; 
   loc_estimate = [loc_estimate; Q_estimate(1)];
%    vel_estimate = [vel_estimate; Q_estimate(2)]; 
   while t == round(mk*sampletimes+1)
    K = Ppred*C'* inv(C*Ppred*C'+ Ez/dd);
    Q_estimate = Q_estimate + K * (loc_meas(mk+1) - C * Q_estimate);
    Pcorr =  (eye(2)-K*C)*Ppred;
    mk = mk+1;
    Pcorr_estimate(2*t-1:2*t,:) = Pcorr;
    loc_estimate(t,1) = Q_estimate(1);
    vel_estimate(t,1) =  Q_estimate(2); 
   end  
    Ppred_estimate = [Ppred_estimate;Ppred];
    Pcorr_estimate = [Pcorr_estimate;Pcorr];
    Q_loc_estimate = [Q_loc_estimate; Q_estimate(1)];
    vel_estimate = [vel_estimate; Q_estimate(2)];
end

%%�˺�����ͼ����Ӧ��ûʲô����
t_1(1)=0;
for n=1:1:length(Q_loc_estimate)-1
t_1(n+1)=t_1(n)+ha;
end
figure(1) %ƽ�����λ��ͼ��
plot(t_1',Q_loc_estimate,'--');
hold on
t_2(1)=0;
for i=1:1:length(disp)-1
t_2(i+1)=t_2(i)+hd;
end
plot(t_2,disp);
% ts2(1)=0;
% for i=1:1:length(adis)-1
% ts2(i+1)=ts2(i)+ha;
% end
% plot(ts2,adis);
xlabel('time(sec)');ylabel('Disp');title('���ٶȡ�λ��������Ϊ10%');
legend('��ʵλ��','�˲����')

for j=1:length(disp)
    rmss(j)=(dis(j)-Q_loc_estimate(j));
end
mean(abs(rmss))
for j=1:length(dis)
    rms_estimate(j)=(dis(j)-loc_estimate(j));
end
rms_loc_extimate=sum(rms_estimate.*rms_estimate);
rms_loc_extimate=(rms_loc_extimate/length(rms_loc_extimate))^0.5
figure(2)
t_3(1)=0;
for i=1:1:length(ndisp)-1
t_3(i+1)=t_3(i)+hd;
end
plot(t_3,ndisp-disp,'--');
hold on
plot(t_1,rms_estimate','Color','black','LineWidth',2);
xlabel('time(sec)');ylabel('Disp error');title('���ٶȡ�λ��������Ϊ10%');
legend('����','�˲����');
