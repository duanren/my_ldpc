function [ber3,bler3]=ldpc_biawgn_bec(H,Ebn0,epsilon,times,rearranged_cols)
%本函数为LDPC码在BIAWGN-BEC信道下进行仿真
%默认使用BPSK调制方法
%输入参数为
%校验矩阵H 信道参数Eb/n0 epsilon 仿真次数times 重排列位置rearranged_cols
%输出结果为
%误码率ber 误包率bler

%配置LDPC编译码器
cfgLDPCEnc = ldpcEncoderConfig(H);
cfgLDPCDec = ldpcDecoderConfig(H);

%准备记录误码率ber和误包率bler
ber3=zeros(length(epsilon),length(Ebn0));
bler3=zeros(length(epsilon),length(Ebn0));


%BPSK调制
bpskmod = comm.BPSKModulator;
%BP最大迭代次数
maxnumiter=100;

%进行times次仿真
for jj=1:length(epsilon)
    ep=epsilon(jj);
    for kk=1:length(Ebn0)
        %Eb/N0与SNR的转换
        snr=Ebn0(kk)+10*log10(cfgLDPCEnc.CodeRate);
       bpskdemod = comm.BPSKDemodulator('DecisionMethod','Approximate log-likelihood ratio', ...
        'Variance',1/10^(snr/10));
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
            %AWGN信道
            Bob_Signal1 = awgn(modSignal,snr);
            %BEC信道
            Bob_Signal2=mybec(Bob_Signal1,ep);
            %BPSK解调
            Bob_demodSignal = bpskdemod(Bob_Signal2);

%             %打孔
%             for hh=1025:2048
%                 Bob_demodSignal(rearranged_cols==hh)=0;
%             end
%             
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
        ber3(jj,kk)=be/times+1e-8;
        bler3(jj,kk)=ble/times+1e-8;
    end
end


