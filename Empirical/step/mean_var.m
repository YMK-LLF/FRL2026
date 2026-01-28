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

% 存储结果
SPM_matrix=zeros(500,1);
SC_matrix=zeros(500,1);
Lshrink_matrix=zeros(500,1);
NLshrink_matrix=zeros(500,1);
SD_matrix=zeros(500,1);

SPM_matrix1=zeros(500,1);
SC_matrix1=zeros(500,1);
Lshrink_matrix1=zeros(500,1);
NLshrink_matrix1=zeros(500,1);
SD_matrix1=zeros(500,1);
train_n=125 %样本量需根据数据进行调整
for i_1 = 1:500
    i_1
 test_sample= data((train_n + i_1):(train_n+29+i_1), :);
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

[Lshrink1, NLshrink1, SPM1, SD1,SC1] = calculate_mean(test_sample, LWEUP, NLWEUP, SPMWEUP, SDWEUP,SCWEUP, L);
    SPM_matrix(i_1)=SPM1;
    SC_matrix(i_1)=SC1;
    Glasso_matrix(i_1)=Glasso1;
    Lshrink_matrix(i_1)=Lshrink1;
    NLshrink_matrix(i_1)=NLshrink1;
    SD_matrix(i_1)=SD1;
    

[Lshrink2, NLshrink2, SPM2, SD2,SC2] = calculate_var(test_sample, LWEUP, NLWEUP, SPMWEUP, SDWEUP,SCWEUP, L);
    SPM_matrix1(i_1)=SPM2;
    SC_matrix1(i_1)=SC2;
    Lshrink_matrix1(i_1)=Lshrink2;
    NLshrink_matrix1(i_1)=NLshrink2;
    SD_matrix1(i_1)=SD2;



end


 

%%



mean_SC=mean(SC_matrix);
mean_Lshrink=mean(Lshrink_matrix);
mean_NLshrink=mean(NLshrink_matrix);
mean_SD=mean(SD_matrix);
mean_SPM=mean(SPM_matrix);
var_SC=mean(SC_matrix1);
var_Lshrink=mean(Lshrink_matrix1);
var_NLshrink=mean(NLshrink_matrix1);
var_SD=mean(SD_matrix1);
var_SPM=mean(SPM_matrix1);

mean_matrix = [ mean_Lshrink, mean_NLshrink, mean_SD, mean_SPM,mean_SC]
var_matrix = [var_Lshrink, var_NLshrink, var_SD, var_SPM,var_SC]



std_SC=std(SC_matrix);
std_Lshrink=std(Lshrink_matrix);
std_NLshrink=std(NLshrink_matrix);
std_SD=std(SD_matrix);
std_SPM=std(SPM_matrix);
std_sig_SC=std(SC_matrix1);
std_sig_Lshrink=std(Lshrink_matrix1);
std_sig_NLshrink=std(NLshrink_matrix1);
std_sig_SD=std(SD_matrix1);
std_sig_SPM=std(SPM_matrix1);


std_matrix = [ std_Lshrink, std_NLshrink, std_SD, std_SPM,std_SC]
std_sig_matrix = [std_sig_Lshrink, std_sig_NLshrink, std_sig_SD, std_sig_SPM,std_sig_SC]
