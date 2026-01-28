clc;
close all;
clear all;

data=readtable('data.csv');
matrix_miu=load('Temp\matrix_miu.mat').matrix_miu;
matrix_SLS=load('S.mat').matrix_SLS;
matrix_SNL=load('S.mat').matrix_SNL;
matrix_SC=load('Temp\matrix_SC.mat').matrix_SC;
matrix_C=load('Temp\matrix_C.mat').matrix_C;
matrix_SD=load('Temp\matrix_SD.mat').matrix_Sp;
train_n=125 %样本量需根据不同数据变化

% 存储结果
SPM_matrix=zeros(500,1);
SC_matrix=zeros(500,1);
Lshrink_matrix=zeros(500,1);
NLshrink_matrix=zeros(500,1);
SD_matrix=zeros(500,1);


for i_1 = 1:500
    i_1
 test_sample= data((train_n + i_1):(train_n-1+i_1), :);
 test_sample = table2array(test_sample);
 % 读取矩阵 
 SLS = matrix_SLS.(sprintf('SLS_%d', i_1));  
 SNL = matrix_SNL.(sprintf('SNL_%d', i_1));

 miu = matrix_miu{i_1};  % 500x1 cell需要指定列
 C=matrix_C{i_1};  % 500x1 cell需要指定列
 SD=matrix_SD{i_1};
 SC=matrix_SC{i_1};

I_SLS=inv(SLS);
I_SNL=inv(SNL);
I_SD=inv(SD);
I_Sc=C;
I_SPM=inv(SC);
L=5;

p=100;
one_p=ones(p,1);
[LWEUP] = calculate_W(I_SLS,L,one_p,miu);
[NLWEUP] = calculate_W(I_SNL,L,one_p,miu);
[SCWEUP] = calculate_W(I_Sc,L,one_p,miu);
[SDWEUP] = calculate_W(I_SD,L,one_p,miu);
[SPMWEUP]=calculate_W(I_SPM,L,one_p,miu);
%%

[Lshrink1, NLshrink1, SPM1, SD1,SC1] = calculate_EUP(test_sample, LWEUP, NLWEUP, SPMWEUP, SDWEUP,SCWEUP, L);
    SPM_matrix(i_1)=SPM1;
    SC_matrix(i_1)=SC1;
    Lshrink_matrix(i_1)=Lshrink1;
    NLshrink_matrix(i_1)=NLshrink1;
    SD_matrix(i_1)=SD1;
    

end


 save('Temp\SC_matrix.mat','SC_matrix');
 save('Temp\Lshrink_matrix.mat','Lshrink_matrix');
 save('Temp\NLshrink_matrix.mat','NLshrink_matrix')
 save('Temp\SD_matrix.mat','SD_matrix')
 save('Temp\SPM_matrix.mat','SPM_matrix')


%%

% 比较SC是否比其他所有列都大（逐行比较）
is_SC_max =(SC_matrix > Lshrink_matrix) & (SC_matrix > NLshrink_matrix) & (SC_matrix > SD_matrix) & (SC_matrix > SPM_matrix);

% 统计SC比其他所有列都大的次数
count_SC_max = sum(is_SC_max);
fprintf('SC比其他所有列都大的次数: %d\n', count_SC_max);


mean_SC=mean(SC_matrix);
mean_Lshrink=mean(Lshrink_matrix);
mean_NLshrink=mean(NLshrink_matrix);
mean_SD=mean(SD_matrix);
mean_SPM=mean(SPM_matrix);
var_SC=var(SC_matrix);
var_Lshrink=var(Lshrink_matrix);
var_NLshrink=var(NLshrink_matrix);
var_SD=var(SD_matrix);
var_SPM=var(SPM_matrix);

mean_matrix = [ mean_Lshrink, mean_NLshrink, mean_SD, mean_SPM,mean_SC]
var_matrix = [var_Lshrink, var_NLshrink, var_SD, var_SPM,var_SC]

