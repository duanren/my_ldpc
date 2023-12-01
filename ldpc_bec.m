function [ber2,bler2]=ldpc_bec(H,epsilon,times,rearranged_cols)
%本函数为LDPC码在BEC信道下进行仿真
%默认使用BPSK调制方法
%输入参数为
%校验矩阵H 信道参数epsilon 仿真次数times 重排列位置rearranged_cols（只针对打孔）
%输出结果为
%误码率ber 误包率bler

%配置LDPC编译码器
cfgLDPCEnc = ldpcEncoderConfig(H);
cfgLDPCDec = ldpcDecoderConfig(H);

%准备记录误码率ber和误包率bler
ber2=zeros(1,length(epsilon));
bler2=zeros(1,length(epsilon));

%BPSK调制
bpskmod = comm.BPSKModulator;
%BP最大迭代次数
maxnumiter=100;

%进行times次仿真
for kk=1:length(epsilon)
    %此处snr赋值较高
    snr=20;
    bpskdemod = comm.BPSKDemodulator('DecisionMethod','Approximate log-likelihood ratio', ...
        'Variance',1/10^(snr/10));
    ep=epsilon(kk);
    %打表
    be=0;
    ble=0;

    for ii=1:times
        %产生随机信息矩阵
        data = randi([0 1],cfgLDPCEnc.NumInformationBits,1);
        %LDPC编码
        encodedData=ldpcEncode(data,cfgLDPCEnc);
        %BPSK调制
        modSignal = bpskmod(encodedData);
        %BEC信道
        Bob_Signal = mybec(modSignal,ep);
        %BPSK解调
        Bob_demodSignal = bpskdemod(Bob_Signal);

% %         %打孔
%         for hh=1:208
%             Bob_demodSignal(rearranged_cols==hh)=0;
%         end

        %LDPC译码，采用LLR-BP译码算法
        Bob_data = ldpcDecode(Bob_demodSignal,cfgLDPCDec,maxnumiter);
        %计算BER
        [~,err1]= biterr(data,Bob_data);
        be=be+err1;
        %若BER不等于0，则误包
        if err1==0
            err2=0;
        else
            err2=1;
        end
        ble=ble+err2;
    end
    %获得误码率和误包率
    ber2(kk)=be/times+1e-8;
    bler2(kk)=ble/times+1e-8;
end


