clc;
close all;
clear all;

a=[];  
b=[];
c=[];
d=[];
D_values=[];
R=[];

column_means = zeros(6, 2);
column_stds=zeros(6,2);
    
% 同比循环增大 p 和 n
for i_2=1:6
    i_2
    p=50*i_2
    train_n = p / 0.5; % 保持 p/n 
    J = p / train_n;  %需要变化
    train_n=fix(train_n)
    I_p=diag(ones(1,p));
    %***************总体协方差矩阵的设计*****************
    sigma_0 = 1;
    lambda1 = [20;10;5]; r1=length(lambda1);  %大spike特征根的个数
    lambda2 = [-0.99]; r2=length(lambda2);        %小spike特征根的个数：只设定一个小的spike特征根
    
    lambda = [lambda1;lambda2]; r = r1+r2;
    v_1= zeros(1,p); % 初始化长度为50的行向量，所有元素都设置为零
    v_1(1) = 1; % 将第一个元素设置为1
    v_2= zeros(1,p); % 初始化长度为50的行向量，所有元素都设置为零
    v_2(2) = 1; % 将第二个元素设置为1
    v_3= zeros(1,p); % 初始化长度为50的行向量，所有元素都设置为零
    v_3(3) = 1; % 将第三个元素设置为1
    vector_1=[v_1;v_2;v_3]';
    v_p=zeros(1,p);
    v_p(p)=1;
    vector_2=[v_p]';

    Sigma = construct_sigma(lambda1, lambda2, r1, r2, I_p, sigma_0, p);
    miu=1/sqrt(p)*ones(p,1);
    
    for i_1=1:500
        i_1
        rng(i_2*1000+i_1);
        %***********************生成训练样本集*************************

        train_sample=mvnrnd(miu,Sigma,train_n);

        %求解样本均值
        hat_mu=mean(train_sample);
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

      %  w_1=[-0.9524;-0.9091;-0.8333];     w_2=[99];
         w_1=[-0.991026, -0.982908, -0.968812];     w_2=[18.795580];        % 0.8
      %  w_1=[-0.976489, -0.955130, -0.917789];     w_2=[48.742460];        % 0.5
      
        % w_1=[-0.959387, -0.922418, -0.857587];     w_2=[85.539785];% 0.2

        Q1=zeros(p,1);Q2=zeros(p,1);
            for i_4=1:r1
                Q1=Q1+w_1(i_4)*vector_S1(:,i_4)*vector_S1(:,i_4)';
            end
            for i_5=1:r2
                Q2=Q2+w_2(i_5)*vector_S2(:,i_5)*vector_S2(:,i_5)';
            end 

        C_hat=1/(sigma_0)^2*(I_p+Q1+Q2);
        one_p=ones(p,1);
        
        [a1,a2]= calculate_aj(lambda1, lambda2, J);
        [b1, b2] = calculate_bj(vector_1, vector_2, p);
        [c1, c2] = calculate_cj(miu, vector_1, vector_2);
        X1=one_p'/sqrt(p)*C_hat*(one_p/sqrt(p));
        [Xhat_1, Xhat_2, Xhat_3] = calculate_Xhat( w_1, w_2, a1, a2, b1, b2, c1, c2, r1, r2, sigma_0, one_p, miu, p);  
        a(i_1)=X1;   
        b(i_1)=Xhat_1;

        Y1=one_p'/sqrt(p)*C_hat*Sigma*C_hat*(one_p/sqrt(p));
        [Yhat_1, Yhat_2, Yhat_3] = calculate_Yhat(lambda1, lambda2, w_1, w_2, a1, a2, b1, b2, c1, c2, r1, r2, sigma_0, one_p, miu, p);
        c(i_1)=Y1;   
        d(i_1)=Yhat_1;

    end
    nGMVP=c./(a.^2);
    nGMVP_hat=d./(b.^2);
    nGMVP=nGMVP';
    nGMVP_hat=nGMVP_hat';
    GMVP=[nGMVP nGMVP_hat];
    %  计算每一列的均值
    column_mean  =  mean(GMVP);
    column_means(i_2,:)=column_mean
    %  计算每一列的标准差
    column_std  =  std(GMVP);
    column_stds(i_2,:)=column_std
    % 计算差值
    D_value=column_mean(1)-column_mean(2)
    D_values(i_2)=D_value
    
    RISK=1/(one_p'*inv(Sigma)*one_p)*p;
    R(i_2)=RISK
    R=R'
end    
