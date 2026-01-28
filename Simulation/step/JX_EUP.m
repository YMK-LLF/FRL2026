clc;
close all;
clear all;

x_1=[];  
xhat_1=[];
y_1=[];
yhat_1=[];

x_2=[];  
xhat_2=[];
y_2=[];
yhat_2=[];

x_3=[];  
xhat_3=[];
y_3=[];
yhat_3=[];

eup=[];
eup_hat=[];
L=5;
D_values=[];
eeu=[];
column_means = zeros(6, 2);
column_stds=zeros(6,2);
    
for i_2=1:6
    i_2
    p=50*i_2
    train_n = p / 0.5; % 保持 p/n 
    J = p / train_n;   %需要变化
    train_n=fix(train_n);
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
    QQ1=zeros(p,1);QQ2=zeros(p,1);
    for i_3=1:r1
        QQ1=QQ1+lambda1(i_3)*vector_1(:,i_3)*vector_1(:,i_3)';
    end
    for i_4=1:r2
        QQ2=QQ2+lambda2(i_4)*vector_2(:,i_4)*vector_2(:,i_4)';
    end
    Sigma = (sigma_0)^2*(I_p+QQ1+QQ2);
    miu=1/sqrt(p)*ones(p,1);



 
    for i_1=1:500
        i_1
        rng(i_2*1000+i_1);
        %***********************生成训练样本集*************************

        train_sample=mvnrnd(miu,Sigma,train_n);


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

        % w_1=[-0.9636, -0.9302, -0.8690];     w_2=[89.7846]; %0.4

     % w_1=[-0.9524;-0.9091;-0.8333];     w_2=[99];
     
        w_1=[-0.9730,-0.9482,-0.8993];     w_2=[50.0576]; 
        % w_1=[-0.9758;-0.9538;-0.912];     w_2=[49.9911];   %609
       


        % w_1=[-0.9669, -0.9350, -0.8740];     w_2=[47.5848]; %0.8
        % w_1=[-0.9588, -0.9212, -0.8553];     w_2=[89.8983]; %0.2
        
        % w_1=[-0.9999, -0.99999, -0.9999];     w_2=[52.9330]; %1.5
        % w_1=[-0.9999, -0.99999, -0.9999];     w_2=[52.9330]; %2

        
        
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
        X_1=one_p'/sqrt(p)*C_hat*(one_p/sqrt(p));
        X_2=one_p'/sqrt(p)*C_hat*(miu/norm(miu));
        X_3=miu'/norm(miu)*C_hat*(miu/norm(miu));
        [Xhat_1, Xhat_2, Xhat_3] = calculate_Xhat( w_1, w_2, a1, a2, b1, b2, c1, c2, r1, r2, sigma_0, one_p, miu, p);  


        Y_1=one_p'/sqrt(p)*C_hat*Sigma*C_hat*(one_p/sqrt(p));
        Y_2=one_p'/sqrt(p)*C_hat*Sigma*C_hat*(miu/norm(miu));
        Y_3=miu'/norm(miu)*C_hat*Sigma*C_hat*(miu/norm(miu));
        [Yhat_1, Yhat_2, Yhat_3] = calculate_Yhat(lambda1, lambda2, w_1, w_2, a1, a2, b1, b2, c1, c2, r1, r2, sigma_0, one_p, miu, p);

        
        

        y_1(i_1)=Y_1;
        y_2(i_1)=Y_2;
        y_3(i_1)=Y_3;
        yhat_1(i_1)=Yhat_1;
        yhat_2(i_1)=Yhat_2;
        yhat_3(i_1)=Yhat_3;

 

        EUP=norm(miu)*X_2/(sqrt(p)*X_1)+norm(miu)^2*X_3/L-norm(miu)^2*X_2^2/(L*X_1)-L*Y_1/(2*p*X_1^2)-norm(miu)*Y_2/(sqrt(p)*X_1)+norm(miu)*Y_1*X_2/(sqrt(p)*X_1^2)-norm(miu)^2*Y_3/(2*L)+norm(miu)^2*Y_2*X_2/(L*X_1)-norm(miu)^2*Y_1*X_2^2/(2*L*X_1^2);
    
        EUPhat=norm(miu)*Xhat_2/(sqrt(p)*Xhat_1)+norm(miu)^2*Xhat_3/L-norm(miu)^2*Xhat_2^2/(L*Xhat_1)-L*Yhat_1/(2*p*Xhat_1^2)-norm(miu)*Yhat_2/(sqrt(p)*Xhat_1)+norm(miu)*Yhat_1*Xhat_2/(sqrt(p)*Xhat_1^2)-norm(miu)^2*Yhat_3/(2*L)+norm(miu)^2*Yhat_2*Xhat_2/(L*Xhat_1)-norm(miu)^2*Yhat_1*Xhat_2^2/(2*L*Xhat_1^2);
        
        eup(i_1)=EUP;
        eup_hat(i_1)=EUPhat;
        eup=eup';
        eup_hat=eup_hat';
        ZEUP=[eup eup_hat];
    end
    %  计算每一列的均值
    column_mean  =  mean(ZEUP);
    column_means(i_2,:)=column_mean
    %  计算每一列的标准差
    column_std  =  std(ZEUP);
    column_stds(i_2,:)=column_std
    % 计算差值
    D_value=column_mean(1)-column_mean(2);
    D_values(i_2)=D_value
    
    
    XX_1=one_p'/sqrt(p)*inv(Sigma)*(one_p/sqrt(p));
    XX_2=one_p'/sqrt(p)*inv(Sigma)*(miu/norm(miu));
    XX_3=miu'/norm(miu)*inv(Sigma)*(miu/norm(miu));
 
    
    EEU=norm(miu)^2*XX_3/(2*L)-L/(2*p*XX_1)+norm(miu)*XX_2/(sqrt(p)*XX_1)-norm(miu)^2*XX_2^2/(2*L*XX_1);
  
    eeu(i_2)=EEU;
    eeu=eeu'
end

