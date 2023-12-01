%%

%本文件为LDPC码测试平台，可在BIAWGN信道、BEC信道和带删除BIAWGN-BEC信道下进行仿真
%默认使用BPSK调制方法
%输入参数为
%校验矩阵H 信道参数Eb/n0或epsilon 仿真次数times
%输出结果为
%误码率ber 误包率bler
%并绘制曲线
%若只需仿真某一种信道，可使用"运行节"选项
%可定点仿真BIAWGN信道，需输入定点位宽w

% clc,clear,close all;
% rng('shuffle');

%校验矩阵H 目标为信息位k=1024,码率为1/3，1/6，1/10，1/20

% %QC-LDPC码可利用基矩阵P和扩展参数blockSize得到，其余LDPC码需直接给出H
% %此处以802.11LDPC(1944,1/2)为例
% P=[57 -1 -1 -1 50 -1 11 -1 50 -1 79 -1 1 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
%     3 -1 28 -1 0 -1 -1 -1 55 7 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1 -1
%     30 -1 -1 -1 24 37 -1 -1 56 14 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1 -1
%     62 53 -1 -1 53 -1 -1 3 35 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1 -1
%     40 -1 -1 20 66 -1 -1 22 28 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1 -1
%     0 -1 -1 -1 8 -1 42 -1 50 -1 -1 8 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1 -1
%     69 79 79 -1 -1 -1 56 -1 52 -1 -1 -1 0 -1 -1 -1 -1 -1 0 0 -1 -1 -1 -1
%     65 -1 -1 -1 38 57 -1 -1 72 -1 27 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1 -1
%     64 -1 -1 -1 14 52 -1 -1 30 -1 -1 32 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1 -1
%     -1 45 -1 70 0 -1 -1 -1 77 9 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0 -1
%     2 56 -1 57 35 -1 -1 -1 -1 -1 12 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 0
%     24 -1 61 -1 60 -1 -1 27 51 -1 -1 16 1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0];
% blockSize=81;
% H=ldpcQuasiCyclicMatrix(blockSize,P);

%% 

%此处可由文件中双击导入H
%若H不是QC-LDPC或典型形式，需化作典型形式
H=double(H);
[rearranged_H,rearranged_cols]=typical_H(H);

%仿真次数times 
times=20000;
%

%%
%BIAWGN信道 
%信道参数Eb/n0
Ebn0=0.8:0.1:1;

[ber1,bler1]=ldpc_biawgn(rearranged_H,Ebn0,times,rearranged_cols);
%绘制图像
figure;
semilogy(Ebn0,ber1);
xlabel('Eb/n0(dB)');
ylabel('BER');
title('BIAWGN');
%legend('802.11(1944,1/2)');
grid on
hold on

figure;
semilogy(Ebn0,bler1);
xlabel('Eb/n0(dB)');
ylabel('BLER');
title('BIAWGN');
%legend('802.11(1944,1/2)');
grid on
hold on

%%


%BEC信道
%信道参数epsilon 0~1之间
epsilon=0.75:0.01:0.85;
%

[ber2,bler2]=ldpc_bec(rearranged_H,epsilon,times,rearranged_cols);
%绘制图像
figure;
semilogy(epsilon,ber2);
xlabel('epsilon');
ylabel('BER');
title('BEC');
%legend('802.11(1944,1/2)');
grid on
hold on

figure;
semilogy(epsilon,bler2);
xlabel('epsilon');
ylabel('BLER');
title('BEC');
%legend('802.11(1944,1/2)');
grid on
hold on



%%


%带删除的BIAWGN-BEC信道
%信道参数Eb/n0
%信道参数epsilon 0~1之间
Ebn0=8;
epsilon=0.78:0.01:0.84;
%

[ber3,bler3]=ldpc_biawgn_bec(rearranged_H,Ebn0,epsilon,times,rearranged_cols);
%绘制图像
figure;
for ii=1:length(Ebn0)
    semilogy(epsilon,ber3(:,ii));
    grid on
    hold on
end
xlabel('epsilon');
ylabel('BER');
title('BIAWGN-BEC');

figure;
for ii=1:length(Ebn0)
    semilogy(epsilon,bler3(:,ii));
    grid on
    hold on
end
xlabel('epsilon');
ylabel('BLER');
title('BIAWGN-BEC');

%% 

%定点BIAWGN信道
%定点位宽w
%信道参数Eb/n0
w=12;
Ebn0=[0.9 1 1.1];
%

[ber4,bler4]=ldpc_fixed_biawgn(rearranged_H,Ebn0,times,w,rearranged_cols);
%绘制图像
figure;
semilogy(Ebn0,ber4);
xlabel('Eb/n0(dB)');
ylabel('BER');
title(strcat('FIXED BIAWGN w=',num2str(w)));
grid on
hold on

figure;
semilogy(Ebn0,bler4);
xlabel('Eb/n0(dB)');
ylabel('BLER');
title(strcat('FIXED BIAWGN w=',num2str(w)));
grid on
hold on




%% 

%定点BEC信道
%定点位宽w
%信道参数epsilon
w=6;
epsilon=0.5:0.01:0.6;
%

[ber5,bler5]=ldpc_fixed_bec(rearranged_H,epsilon,times,w,rearranged_cols);
%绘制图像
figure;
semilogy(epsilon,ber5);
xlabel('epsilon');
ylabel('BER');
title('FIXED BEC w=6');
%legend('802.11(1944,1/2)');
grid on
hold on

figure;
semilogy(epsilon,bler5);
xlabel('epsilon');
ylabel('BLER');
title('FIXED BEC w=6');
%legend('802.11(1944,1/2)');
grid on
hold on

%% 

%定点BIAWGN-BEC信道
%定点位宽w
%信道参数Eb/n0
%信道参数epsilon 0~1之间
w=6;
Ebn0=8;
epsilon=0.75:0.01:0.8;
%

[ber6,bler6]=ldpc_fixed_biawgn_bec(rearranged_H,Ebn0,epsilon,times,w,rearranged_cols);
%绘制图像
figure;
for ii=1:length(Ebn0)
    semilogy(epsilon,ber6(:,ii));
    grid on
    hold on
end
xlabel('epsilon');
ylabel('BER');
title('FIXED BIAWGN-BEC w=6');

figure;
for ii=1:length(Ebn0)
    semilogy(epsilon,bler6(:,ii));
    grid on
    hold on
end
xlabel('epsilon');
ylabel('BLER');
title('FIXED BIAWGN-BEC w=6');

