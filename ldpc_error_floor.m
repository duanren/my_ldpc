% 清除历史输出和工作区变量
clc,clear,close all;
% 打乱随机种子
rng('shuffle');

%从当前目录下H加载矩阵
load("H_1_3.mat");
% load("H_1_6.mat");
% load("H_1_10.mat");
% load("H_1_20.mat");

%若H不是QC-LDPC或典型形式，需化作典型形式
H=double(H);
[rearranged_H,rearranged_cols]=typical_H(H);

%配置LDPC编译码器
cfgLDPCEnc = ldpcEncoderConfig(rearranged_H);
cfgLDPCDec = ldpcDecoderConfig(rearranged_H,'offset-min-sum');

%仿真次数times
times=1e9;

%BIAWGN信道
%信道参数Eb/n0
Ebn0=0:0.1:1.5;%仿真起点
n=length(Ebn0);
%记录误比特数be、误包数ble、误比特率ber,误包率bler,总迭代次数ite和平均迭代次数iter
be=zeros(1,n);
ble=zeros(1,n);
ber=zeros(1,n);
bler=zeros(1,n);
minble=100; % 错100个包时停止
ite=zeros(1,n);
iter=zeros(1,n);

%BPSK调制
bpskmod = comm.BPSKModulator;
%BP最大迭代次数
maxnumiter=200;

%遍历Ebn0
for kk=1:n
    %Eb/N0与SNR的转换
    if cfgLDPCEnc.NumInformationBits==1040
        snr=Ebn0(kk)+10*log10(1/3);
    else
        snr=Ebn0(kk)+10*log10(cfgLDPCEnc.CodeRate);
    end
    %BPSK解调
    bpskdemod = comm.BPSKDemodulator('DecisionMethod','Approximate log-likelihood ratio', ...
        'Variance',1/10^(snr/10));
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

        %LDPC译码，采用LLR-BP译码算法
        [Bob_data,Bob_iter,~] = ldpcDecode(Bob_demodSignal,cfgLDPCDec,maxnumiter);
        %累计单次误比特率
        if cfgLDPCEnc.NumInformationBits==1040
            [~,err1]= biterr(data(1:1024),Bob_data(1:1024));
        else
            [~,err1]= biterr(data,Bob_data);
        end
        be(kk)=be(kk)+err1;
        %累计误包数
        if err1==0
            err2=0;
        else
            err2=1;
        end
        ble(kk)=ble(kk)+err2;
        %累计迭代次数
        ite(kk)=ite(kk)+Bob_iter;
        %计算误码率,误包率和平均迭代次数
        ber(kk)=be(kk)/ii;
        bler(kk)=ble(kk)/ii;
        iter(kk)=ite(kk)/ii;

        %每1000个点打印当前结果
        if mod(ii,1000)==0
            disp(strcat("times:",num2str(ii)))
            disp("ber:")
            disp(ber)
            disp("bler:")
            disp(bler)
            disp("iter:")
            disp(iter)
        end
        % 保存数据
        if ble(kk)>minble
            disp(strcat("times:",num2str(ii)))
            disp("ber:")
            disp(ber)
            disp("bler:")
            disp(bler)
            disp("iter:")
            disp(iter)
            save("error_floor.mat","minble","Ebn0","be","ble","ber","bler","ite","iter");
            break;
        end
    end
end

%绘制曲线
load("error_floor.mat");
%%
figure;
semilogy(Ebn0,ber);
xlabel('Eb/n0(dB)');
ylabel('BER');
title('BIAWGN');
legend(strcat('ble:',num2str(minble)));
grid on
hold on
savefig("error_floor_ber.fig")

figure;
semilogy(Ebn0,bler);
xlabel('Eb/n0(dB)');
ylabel('BLER');
title('BIAWGN');
legend(strcat('ble:',num2str(minble)));
grid on
hold on
savefig("error_floor_bler.fig")

figure;
plot(Ebn0,iter);
xlabel('Eb/n0(dB)');
ylabel('ITER');
title('BIAWGN');
legend(strcat('ble:',num2str(minble)));
grid on
hold on
savefig("error_floor_iter.fig")