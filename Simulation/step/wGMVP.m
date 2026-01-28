clc;
close all;
clear all;


    p=300;
    J=0.5;
    train_n=p/J;
   

    run('step1');
    save('temp/step1.mat','L','p','train_n','J','I_p','mu','Sigma','lambda1','lambda2','r1','r2','one_p','sigma_0','vector_1','vector_2','miu');
    [w1_opt, w2_opt, w3_opt, w4_opt, f_val] = optimize_weights_GMVP();
    w_1=[w1_opt, w2_opt, w3_opt];w_2=[w4_opt];
   
    train_n=fix(train_n);



    data=load('Temp/step1.mat');
    I_p=diag(ones(1,p));
    WEUP=zeros(p,500);
        %***************总体协方差矩阵的设计*****************
    sigma_0 = 1;
    lambda1 = [20;10;5]; r1=length(lambda1);  %大spike特征根的个数
    lambda2 = [-0.99]; r2=length(lambda2);        %小spike特征根的个数：只设定一个小的spike特征根

    lambda = [lambda1;lambda2]; r = r1+r2;


    Sigma=data.Sigma;
    miu=data.miu;
    one=ones(p,1);
    for i_1=1:500
        i_1
        rng(1000+i_1);
        %***********************生成训练样本集*************************

        train_sample=mvnrnd(miu,Sigma,train_n);
%%
        %*********对训练样本协方差阵进行特征分解***********
        Sp=cov(train_sample);

        %********对整合的协方差阵进行特征分解**********
        [vector_Sp0,lammda_Sp0]=eig(Sp);
        %**********对特征根降序排列**************************
        [lammda_Sp,index]=sort(diag(lammda_Sp0),'descend');
        %*********对特征向量按照特征根的顺序排列********
        vector_Sp = vector_Sp0(:,index);

        %************大小spike特征值和相应的特征向量*************
        vector_S1=vector_Sp(:,1:r1);
        vector_S2=vector_Sp(:,p-r2+1:p);

 
        Q1=zeros(p,1);Q2=zeros(p,1);
            for i_4=1:r1
                Q1=Q1+w_1(i_4)*vector_S1(:,i_4)*vector_S1(:,i_4)';
            end
            for i_5=1:r2
                Q2=Q2+w_2(i_5)*vector_S2(:,i_5)*vector_S2(:,i_5)';
            end 

        C=1/(sigma_0)^2*(I_p+Q1+Q2);
%%
    % ********** 分开计算分子和分母 **********
    numerator = C * one;              % 分子部分
    denominator = one' * C * one;    % 分母部分

    WGMVP1 = numerator / denominator;

 
%%

    WGMVP(:,i_1)=WGMVP1;

    end




   save('Temp\wGMVP\J_0.5\WGMVP_p300.mat','WGMVP');

