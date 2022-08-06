function [freq,modal_shape,damping] = cov_ssi_simplify()
% cov_ssi方法实现模态识别的程序（环境激励下的模态识别）
% 数据输入在程序内部进行，需要通过读入数据表来实现
% 输出：freq,modal_shape,damping 依次为频率、振型和阻尼比矩阵（向量）
% 读入数据
% 下面这三行代码在读入excel数据时，需要数据格式为：第一列为时间序列，单位为秒，之后每一列表示一个测点的加速度、速度或者位移数据
[FileName,PathName]=uigetfile('*.xlsx','请选择要整合数据的挠度时程表格');  %显示选择时程文件的对话框
File=strcat(PathName,FileName);  %把数横向连接成单个字符串
time_history_0 = xlsread(File);   % 该参数需要数据格式为：第一列为时间序列，单位为秒，之后每一列表示一个测点的加速度、速度或者位移数据
[m,n] = size(time_history_0);
if m < 10 || n < 2
    key_of_source = questdlg('检测到输入数据长度过小或没有测点数据，是否依旧进行参数识别？','是否依旧进行参数识别？','是','否','是');
    if key_of_source == '是'
        time_history = time_history_0;
    else
        freq = [];modal_shape = [];damping = [];
        return
    end
else
    time_history = time_history_0;
end

T = time_history(3,1) - time_history(2,1);    % 计算采样时间间隔
fs = 1/T;         % 计算采样频率

% cov_ssi 模态识别的主程序，来源于网络，时间太长，忘记出处了，还请见谅
[TimePointNum,MesNodeTotalNum]=size(time_history_0);
MesNodeTotalNum=MesNodeTotalNum-1;

prompt={'TOEPLITZ矩阵块数M(>=采样频率/(2*结构基频)):','系统阶数N(<i*实际测点数):'};
dlg_title='请输入COV_SSI方法计算的必要参数(最好为2的整倍数)';
num_lines=[1 100;1 100];
def={'50','30'};                                 %系统阶数N=2*n
answer=inputdlg(prompt,dlg_title,num_lines,def,'on');

M=str2double(answer{1});
N=str2double(answer{2});
disp('----------模态参数识别开始----------');
disp('          构建Hankel矩阵中...');
pause(0.000001);  %让程序暂停一点点时间，以刷新界面
j=TimePointNum-2*M;
Hankel=zeros((2*M+1)*MesNodeTotalNum,j);
for ti=1:j  %这个位置用for循环的效率比较低，用数组来实现矩阵的构建可能效率会高一点
    m=0;
    for tii=1:(2*M+1)
        km=ti+tii-1;
        m=m+1;
        %计算每次赋值在所在行的起始位置
        beginy=(m-1)*MesNodeTotalNum+1;
        %计算每次赋值在所在行的结束位置
        endy=m*MesNodeTotalNum;
        Hankel(beginy:endy,ti)=(time_history_0(km,2:(MesNodeTotalNum+1)))';
    end
end

Yp=Hankel(1:M*MesNodeTotalNum,:);
Yf=Hankel(M*MesNodeTotalNum+1:2*M*MesNodeTotalNum,:);
Yf2=Hankel((M+1)*MesNodeTotalNum+1:(2*M+1)*MesNodeTotalNum,:);
%组成第一个Toeplitz矩阵
Teop1=Yf*Yp'/j;
%组成第二个Toeplitz矩阵
Teop2=Yf2*Yp'/j;
disp('          正在进行奇异值分解...');
pause(0.000001);  %让程序暂停一点点时间，以刷新界面
[U,S,V]=svd(Teop1);

%提取前N阶奇异值分解结果
U_de=U(:,1:N);
S_de=S(1:N,1:N);
V_de=V(:,1:N);

%计算C,A矩阵的估计
O_GJ=U_de*S_de^0.5;
C_GJ=O_GJ(1:MesNodeTotalNum,:);
A_GJ=S_de^-0.5*(U_de')*Teop2*V_de*S_de^-0.5;

disp('         奇异值分解完成，进行模态参数识别...');
pause(0.000001);  %让程序暂停一点点时间，以刷新界面
%计算状态矩阵A的特征值和特征向量
[V_of_A,D_of_A]=eig(A_GJ);


%计算系统的特征值和特征向量
% fs=user.data.sim_freq;
% T=1/user.data.sim_freq;
Total_mode_jieShu=N;
%计算振型
MODE=C_GJ*V_of_A;
%计算频率
LN_lamuda = zeros(N,1);
Re_Abs = zeros(N,1);
Im_Abs = zeros(N,1);
Freq = zeros(N,1);
Dpro = zeros(N,1);
for i = 1:N
    LN_lamuda(i,1)=log(D_of_A(i,i))/T;
    Re_Abs(i,1)=abs(real(LN_lamuda(i,1)));
    Im_Abs(i,1)=abs(imag(LN_lamuda(i,1)));
    Freq(i,1)=abs(LN_lamuda(i,1))/2/pi;      %这里计算得到的频率单位为HZ
    %计算阻尼比
    Dpro(i,1)=Re_Abs(i,1)/abs(LN_lamuda(i,1));
end
disp('          识别结果判别中...');
pause(0.000001);  %让程序暂停一点点时间，以刷新界面
%将实部和虚步分别归一化
for opp_index=1:Total_mode_jieShu
    min_of_real=min(real(MODE(:,opp_index)));
    max_of_real=max(real(MODE(:,opp_index)));
    min_of_imag=min(imag(MODE(:,opp_index)));
    max_of_imag=max(imag(MODE(:,opp_index)));
    if max_of_real~=min_of_real
        Mode_value_Total_of_real_GYH=2*(real(MODE(:,opp_index))-min_of_real)/(max_of_real-min_of_real)-1;
    else
        Mode_value_Total_of_real_GYH=zeros(MesNodeTotalNum,1);
    end
    if max_of_imag~=min_of_imag
        Mode_value_Total_of_imag_GYH=2*(imag(MODE(:,opp_index))-min_of_imag)/(max_of_imag-min_of_imag)-1;
    else
        Mode_value_Total_of_imag_GYH=zeros(MesNodeTotalNum,1);
    end
    MODE(:,opp_index)=Mode_value_Total_of_real_GYH+Mode_value_Total_of_imag_GYH*sqrt(-1);
    if imag(MODE(:,opp_index))~=0
        MODE(:,opp_index)=imag(MODE(:,opp_index));
    else
        MODE(:,opp_index)=real(MODE(:,opp_index));
    end
end

% 真实模态判断
% 用阻尼比剔除虚假模态，若是不知道阻尼比大范围，可以输入一个比较大的值
prompt={'请输入真实模态阻尼比最大值（用于剔除“虚假模态”）:'};  %j>20M
dlg_title='Input Parameter of Calculate';
num_lines=1;
def={'0.05'};
answer=inputdlg(prompt,dlg_title,num_lines,def);
user.max_of_damping=str2double(answer{1,1});
M2=length(Freq);
mm2=1;
for m2=1:M2
    if (Freq(m2)<=(0.5*fs))&&(Freq(m2)>0)&&(Dpro(m2)>0)&&(Dpro(m2)<user.max_of_damping)      % 这里0.05需要由用户输入
        Freq2(mm2)=Freq(m2);
        Dpro2(mm2)=Dpro(m2);
        MODE2(:,mm2)=MODE(:,m2);
        mm2=mm2+1;
    end
end

[Freq3,IX]=sort(Freq2);
for m3=1:length(Freq3)
    Dpro3(m3)=Dpro2(IX(m3));
    MODE3(:,m3)=MODE2(:,IX(m3));
end

% 重频判断系数
crtifreq = 1;              % 频率判断标准
crtidamping = 5;           % 阻尼比判断标准
crtimodeshape = 99.9;      % 模态判断标准（MAC）
M3=length(Freq3);
dualfreq=ones(1,M3);
for m3=1:M3
    for n3=m3+1:M3
        X=MODE3(:,m3);
        Y=MODE3(:,n3);
        MAC_shape=(X'*Y)^2/((X'*X)*(Y'*Y));             % 模态保证准则，表征两阶模态的相关性
        if (((Freq3(n3)-Freq3(m3))*100/Freq3(m3))<crtifreq)...
                && (abs(((Dpro3(m3)-Dpro3(n3))*100/Dpro3(m3))) < crtidamping)...
                && (abs(MAC_shape*100) > crtimodeshape)
            dualfreq(n3)=0;
        end
    end
end
indicefreq=find(dualfreq);
M4=length(indicefreq);
for m4=1:M4
    Freq4(m4)=Freq3(indicefreq(m4));
    Dpro4(m4)=Dpro3(indicefreq(m4));
    MODE4(:,m4)=MODE3(:,indicefreq(m4));
end
disp('----------识别完成----------');
pause(0.000001);  %让程序暂停一点点时间，以刷新界面


% 模态参数识别结果列表显示
user.results.freq=Freq4; user.results.damping=Dpro4; user.results.mode_shape=MODE4;
% 输出全部自振频率
[~,n]=size(user.results.freq);
c_n={('阶数'),('频率(Hz)'),('阻尼比(%)')};
b=zeros(n,3);
for i=1:n
    b(i,1)=i;
    b(i,2)=Freq4(i);
    b(i,3)=100*Dpro4(i);
end
uitable('Tag', 'out_table',...
    'Units', 'Normalized',...
    'Data',b,...
    'ColumnName',c_n,...
    'FontSize',12,...
    'Position', [0.02,0.15,0.8,0.8]);
pause(0.000001);  %让程序暂停一点点时间，以刷新界面
freq = user.results.freq;
damping = user.results.damping;
modal_shape = user.results.mode_shape;

% 以下为稳定图计算与绘制程序
key_of_source = questdlg('检测到识别已经完成，是否进行稳定图绘制？','是否进行稳定图绘制？','是','否','是');
if key_of_source == '否'
    return
else
    p = 1;        % 稳定图数据存储变量的行计数变量
    N_min = 2;     % 绘制稳定图的最小阶数
    if rem(N,2) == 0
        N_max = N;    % 绘制稳定图的最大阶数
    else
        N_max = N+1;    % 绘制稳定图的最大阶数
    end
    
    S_freq=-ones((N_max-N_min)/2+1,N_max);
    S_damp=-ones((N_max-N_min)/2+1,N_max);
    
    for N_system=N_min:2:N_max
        U_N=U(:,1:N_system);
        V_N=V(:,1:N_system);
        S_N=S(1:N_system,1:N_system);
        A_N=S_N^-0.5*(U_N')*Teop2*V_N*S_N^-0.5;
        [~,D_N]=eig(A_N);
        [m,~]=size(D_N);
        Eignvalue = zeros(1,m);
        for i=1:m
            Eignvalue(1,i)=D_N(i,i);
        end
        Eignvalue_Ac=log(Eignvalue)/T;
        Freq_N=abs(Eignvalue_Ac)/2/pi;
        Damp_N=-real(Eignvalue_Ac)./abs(Eignvalue_Ac);
        for q=1:N_system
            S_freq(p,q)=Freq_N(1,q);
            S_damp(p,q)=Damp_N(1,q);
        end
        p=p+1;
    end
    
    
    % 下面这段之前是注释掉的，2021年4月28日修改后取消注释,用阻尼比控制来剔除虚假模态
    [b_max,a_max] = size(S_damp);
    try
        for a=1:a_max
            for b=1:b_max
                if S_damp(b,a)>=user.max_of_damping || S_damp(b,a)<=0
                    S_freq(b,a)=-100;
                end
            end
        end
    catch
        prompt={'请输入真实模态阻尼比最大值（用于剔除“虚假模态”）:'};  %j>20M
        dlg_title='Input Parameter of Calculate';
        num_lines=1;
        def={'0.05'};
        answer=inputdlg(prompt,dlg_title,num_lines,def);
        user.max_of_damping=str2double(answer{1,1});
        for a=1:a_max
            for b=1:b_max
                if S_damp(b,a)>=user.max_of_damping || S_damp(b,a)<=0
                    S_freq(b,a)=-100;
                end
            end
        end
    end
    
    
    % *******************************绘制稳定图********************************
    [m,n] = size(time_history_0);
    y_axis=N_min:2:N_max;        % 稳定图y轴的坐标
    while true
        prompt={'请输入PSD谱计算时的窗长度:'};
        dlg_title='Input Parameter of Calculate';
        num_lines=1;
        def={'1024'};
        answer=inputdlg(prompt,dlg_title,num_lines,def);
        correction=str2double(answer{1,1});
        if correction > 0 && 2*correction < m
            nfft = correction;
            break
        else
            disp('          PSD谱计算时的窗长度输入值错误！');
            disp('          请重新输入小于识别数据长度一半的正整数！');
        end
    end
    
    [Pxx,f]=pwelch(time_history_0(:,2:n),hanning(nfft),nfft/2,nfft,fs);
    pxx = mean(Pxx,2);
    figure0 = figure;
    hp = axes('parent',figure0,...
        'FontSize',10,...
        'Units', 'Normalized',...
        'Nextplot', 'Add',...
        'Box', 'on',...
        'Position',[0.15 0.15 0.8 0.8]);
    yyaxis(hp,'left')
    plot(hp,S_freq,y_axis,'*');
    set(hp,'xlim',[0 fs/2]);   %设置坐标轴刻度
    set(hp,'ylim',[0 N_max]);
    set(hp.XLabel,'string','Frequency(Hz)','FontSize',10);   %设置坐标轴标签
    set(hp.YLabel,'string','Number of Singular Values(N)','FontSize',10);
    yyaxis(hp,'right')
    plot(hp,f,pxx,'-');
    set(hp,'xlim',[0 fs/2]);   %设置坐标轴刻度
    set(hp.XLabel,'string','Frequency(Hz)','FontSize',10);   %设置坐标轴标签
    set(hp.YLabel,'string','Power Spectral Density(dB)','FontSize',10);
    grid on
    title('稳定图','FontSize',10);
    disp('----------稳定图绘制完成！----------')
end

end