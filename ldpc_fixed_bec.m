function [ber5,bler5]=ldpc_fixed_bec(H,epsilon,times,w,rearranged_cols)
%本函数为LDPC码在BIAWGN信道下进行定点仿真
%默认使用BPSK调制方法
%输入参数为
%校验矩阵H 信道参数Eb/n0 仿真次数times 定点位宽w 重排列位置rearranged_cols
%输出结果为
%误码率ber 误包率bler

%配置LDPC编译码器
cfgLDPCEnc = ldpcEncoderConfig(H);
cfgLDPCDec = ldpcDecoderConfig(H,'norm-min-sum');

%准备记录误码率ber和误包率bler
ber5=zeros(1,length(epsilon));
bler5=zeros(1,length(epsilon));

%BPSK调制
bpskmod = comm.BPSKModulator;
%BP最大迭代次数
maxnumiter=200;

%进行times次仿真
for kk=1:length(epsilon)
    %Eb/N0与SNR的转换
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

        %定点化
        Bob_FixedSignal=fi(Bob_Signal,1,w);


        %BPSK解调+定点化
        Bob_demodSignal = fi(bpskdemod(Bob_FixedSignal.double),1,w);
        %打孔
        for hh=1:208
            Bob_demodSignal(rearranged_cols==hh)=0;
        end
        for hh=3281:3312
            Bob_demodSignal(rearranged_cols==hh)=0;
        end


        %LDPC译码，采用LLR-MS译码算法
        Bob_data = ldpcDecode(Bob_demodSignal.double,cfgLDPCDec,maxnumiter);
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
    ber5(kk)=be/times+1e-8;
    bler5(kk)=ble/times+1e-8;
end

