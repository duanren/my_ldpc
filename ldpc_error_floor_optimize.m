% 清除历史输出和工作区变量
clc,clear,close all;
% 打乱随机种子
rng('shuffle');

%从当前目录下H加载矩阵
% load("H_1_3.mat");
% load("H_1_6.mat");
% load("H_1_10.mat");
load("H_1_20.mat");

%若H不是QC-LDPC或典型形式，需化作典型形式
H=double(H);
[rearranged_H,rearranged_cols]=typical_H(H);

%配置LDPC编译码器
cfgLDPCEnc = ldpcEncoderConfig(rearranged_H);
cfgLDPCDec1 = ldpcDecoderConfig(rearranged_H,'norm-min-sum');
cfgLDPCDec2 = ldpcDecoderConfig(rearranged_H,'offset-min-sum');

%仿真次数times
times=1e5;

%BIAWGN信道
%信道参数Eb/n0
Ebn0=1;%仿真起点
alpha=[0.1 0.2 0.3 0.4 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.9];
beta=alpha;
n=length(beta);
%记录误比特数be、误包数ble、误比特率ber,误包率bler
be1=zeros(1,n);
ble1=zeros(1,n);
ber1=zeros(1,n);
bler1=zeros(1,n);
be2=zeros(1,n);
ble2=zeros(1,n);
ber2=zeros(1,n);
bler2=zeros(1,n);
minble=100; % 错100个包时停止

%BPSK调制
bpskmod = comm.BPSKModulator;
%BP最大迭代次数
maxnumiter=100;

%Eb/N0与SNR的转换
if cfgLDPCEnc.NumInformationBits==1040
    snr=Ebn0+10*log10(1/3);
else
    snr=Ebn0+10*log10(cfgLDPCEnc.CodeRate);
end
%BPSK解调
bpskdemod = comm.BPSKDemodulator('DecisionMethod','Approximate log-likelihood ratio', ...
    'Variance',1/10^(snr/10));

for kk=1:n
    aa=alpha(kk);
    bb=beta(kk);
    %仿真点数
    for ii=1:times
        %产生随机信息矩阵
        data = randi([0 1],cfgLDPCEnc.NumInformationBits,1);
        %缩短比特初始化(5G)
        if cfgLDPCEnc.NumInformationBits==1040
            for hh=1025:1040
                data(hh)=0;
            end
        end
        %LDPC编码
        encodedData=ldpcEncode(data,cfgLDPCEnc);
        %BPSK调制
        modSignal = bpskmod(encodedData);
        %AWGN信道
        Bob_Signal = awgn(modSignal,snr);
        %BPSK解调
        Bob_demodSignal = bpskdemod(Bob_Signal);

        %预处理(5G)
        if cfgLDPCEnc.NumInformationBits==1040
            %打孔
            for hh=1:208
                Bob_demodSignal(rearranged_cols==hh)=0;
            end
            for hh=3297:3328
                Bob_demodSignal(rearranged_cols==hh)=0;
            end
            %缩短
            for hh=1025:1040
                Bob_demodSignal(rearranged_cols==hh)=realmax;
            end
        end

        %LDPC译码，采用LLR-MS译码算法
        Bob_data1 = ldpcDecode(Bob_demodSignal,cfgLDPCDec1,maxnumiter,'MinSumScalingFactor',aa);
        %累计单次误比特率
        if cfgLDPCEnc.NumInformationBits==1040
            [~,err11]= biterr(data(1:1024),Bob_data1(1:1024));
        else
            [~,err11]= biterr(data,Bob_data1);
        end
        be1(kk)=be1(kk)+err11;

        %累计误包数
        if err11==0
            err12=0;
        else
            err12=1;
        end
        ble1(kk)=ble1(kk)+err12;

        %计算误码率,误包率和平均迭代次数
        ber1(kk)=be1(kk)/ii;
        bler1(kk)=ble1(kk)/ii;


        %每1000个点打印当前结果
        if mod(ii,1000)==0
            disp(strcat("times:",num2str(ii)))
            disp("ber1:")
            disp(ber1)
            disp("bler1:")
            disp(bler1)
            disp("ber2:")
            disp(ber2)
            disp("bler2:")
            disp(bler2)
        end

        if ble1(kk)>minble
            break;
        end
    end

    for ii=1:times
        %产生随机信息矩阵
        data = randi([0 1],cfgLDPCEnc.NumInformationBits,1);
        %缩短比特初始化(5G)
        if cfgLDPCEnc.NumInformationBits==1040
            for hh=1025:1040
                data(hh)=0;
            end
        end
        %LDPC编码
        encodedData=ldpcEncode(data,cfgLDPCEnc);
        %BPSK调制
        modSignal = bpskmod(encodedData);
        %AWGN信道
        Bob_Signal = awgn(modSignal,snr);
        %BPSK解调
        Bob_demodSignal = bpskdemod(Bob_Signal);

        %预处理(5G)
        if cfgLDPCEnc.NumInformationBits==1040
            %打孔
            for hh=1:208
                Bob_demodSignal(rearranged_cols==hh)=0;
            end
            for hh=3297:3328
                Bob_demodSignal(rearranged_cols==hh)=0;
            end
            %缩短
            for hh=1025:1040
                Bob_demodSignal(rearranged_cols==hh)=realmax;
            end
        end

        %LDPC译码，采用LLR-MS译码算法
        Bob_data2 = ldpcDecode(Bob_demodSignal,cfgLDPCDec2,maxnumiter,'MinSumOffset',bb);
        %累计单次误比特率

        if cfgLDPCEnc.NumInformationBits==1040
            [~,err21]= biterr(data(1:1024),Bob_data2(1:1024));
        else
            [~,err21]= biterr(data,Bob_data2);
        end
        be2(kk)=be2(kk)+err21;
        %累计误包数

        if err21==0
            err22=0;
        else
            err22=1;
        end
        ble2(kk)=ble2(kk)+err22;
        %计算误码率,误包率和平均迭代次数

        ber2(kk)=be2(kk)/ii;
        bler2(kk)=ble2(kk)/ii;


        %每1000个点打印当前结果
        if mod(ii,1000)==0
            disp(strcat("times:",num2str(ii)))
            disp("ber1:")
            disp(ber1)
            disp("bler1:")
            disp(bler1)
            disp("ber2:")
            disp(ber2)
            disp("bler2:")
            disp(bler2)
        end

        if ble2(kk)>minble
            break;
        end
    end
    % 保存数据
    disp(strcat("times:",num2str(ii)))
    disp("ber1:")
    disp(ber1)
    disp("bler1:")
    disp(bler1)
    disp("ber2:")
    disp(ber2)
    disp("bler2:")
    disp(bler2)
    save("error_floor_optimize.mat","minble","Ebn0","alpha","beta","be1","ble1","ber1","bler1","be2","ble2","ber2","bler2");
end

%绘制曲线
load("error_floor_optimize.mat");
%%
figure;
semilogy(alpha,ber1,beta,ber2);
xlabel('alpha/beta');
ylabel('BER');
title(strcat('ble:',num2str(minble),' Ebn0:',num2str(Ebn0),'dB'));
legend('NMS','OMS');
grid on
hold on
savefig("error_floor_optimize_ber.fig")

figure;
semilogy(alpha,bler1,beta,bler2);
xlabel('alpha/beta');
ylabel('BLER');
title(strcat('ble:',num2str(minble),' Ebn0:',num2str(Ebn0),'dB'));
legend('NMS','OMS');
grid on
hold on
savefig("error_floor_optimize_bler.fig")