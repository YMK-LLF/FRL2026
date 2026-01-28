clc;
close all;
clear all;

w=load('//D1/w/J_0.5/wp300.mat');
WEUP=load('Temp/J_0.5/WEUP_p300.mat');%调参后
pWEUP=load('Temp/J_0.5/WEUP_pp300.mat');

LWEUP=w.LWEUP;
NLWEUP=w.NLWEUP;
pWEUP=pWEUP.WEUP;
SWEUP=w.SWEUP;
WEUP=WEUP.WEUP;
% 使用 data.x, data.y 等访问导入的变量

    
    p=300;
    L=5;
    data=load('Temp/step1.mat');
    Sigma=data.Sigma;
    miu=data.miu;
    one=ones(p,1);
    test_n=10000;
    

    SP_matrix=zeros(500,1);
    Spm_matrix=zeros(500,1);
    Lshrink_matrix=zeros(500,1);
    NLshrink_matrix=zeros(500,1);
    SD_matrix=zeros(500,1);

    for i_1=1:500
        i_1
        %***********************生成测试样本集*************************
        rng(100)
        test_sample=mvnrnd(miu,Sigma,test_n);
       
%%

        [Spm1,Lshrink1, NLshrink1,SP1,SD1] = calculate_EUP(test_sample, WEUP, LWEUP, NLWEUP,pWEUP,SWEUP,L,i_1);
  
%%

        SP_matrix(i_1)=SP1;
        Spm_matrix(i_1)=Spm1;
        Lshrink_matrix(i_1)=Lshrink1;
        NLshrink_matrix(i_1)=NLshrink1;
        SD_matrix(i_1)=SD1;

    end
    mean_SP = mean(SP_matrix);
    mean_Spm = mean(Spm_matrix);
    mean_Lshrink = mean(Lshrink_matrix);
    mean_NLshrink = mean(NLshrink_matrix);
    mean_SD = mean(SD_matrix);

    sd_SP = std(SP_matrix);
    sd_Spm = std(Spm_matrix);
    sd_Lshrink = std(Lshrink_matrix);
    sd_NLshrink = std(NLshrink_matrix);
    sd_SD = std(SD_matrix);

aResult = [{'Lshrink','NLshrink','SD','SEM', 'SCME'}; num2cell([mean_Lshrink, mean_NLshrink,mean_SD,mean_SP,mean_Spm])];
  