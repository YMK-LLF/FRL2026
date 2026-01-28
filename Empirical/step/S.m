clc;
close all;
clear all;

data=readtable('data.csv'); %读取数据
p=100; % 维度
train_n=125;  %样本量需根据不同数据调整
J=p/train_n;  %维容比
one_p=ones(p,1); %p维1向量
I_p= diag(ones(p, 1));  % 将对角线填充为1
L=5;

% 预分配cell数组
matrix_C = cell(500, 1); 
matrix_SC = cell(500, 1); 
matrix_miu=cell(500,1);
matrix_lambda_large=cell(500,1);
matrix_lambda_small=cell(500,1);
matrix_Sp = cell(500, 1); 

for i_1 = 1:500
    i_1
    train_sample=data(i_1:(train_n-1+ i_1), :);
    train_sample= table2array(train_sample);
    Sp=cov(train_sample);
    miu=mean(train_sample);
    miu=miu';
    [eta,coeff,k,SW] = Rho(Sp,1,p,1,train_n);

    large_number = k(1); % 提取第一个
    small_number = k(2); 
    sigma_0=sqrt(coeff);
    lambda_large=eta(1:large_number);
    lambda_small=eta(large_number+1:end);

    r1=large_number;r2=small_number;
    %********对整合的协方差阵进行特征分解**********
        [vector_Sp0,lammda_Sp0]=eig(Sp);
        %**********对特征根降序排列**************************
        [lammda_Sp,index]=sort(diag(lammda_Sp0),'descend');
        %*********对特征向量按照特征根的顺序排列********
        vector_Sp = vector_Sp0(:,index);

        %************大小spike特征值和相应的特征向量*************
        vector_S1=vector_Sp(:,1:r1);
        vector_S2=vector_Sp(:,p-r2+1:p);
        lambda1=lambda_large;  lambda2=lambda_small;
        xita1=-lambda1./(lambda1+1);xita2=-lambda2./(lambda2+1);
        xita2 = xita2(:);xita1 = xita1(:);
% w_1=xita1;w_2=xita2;
%%
%调参

    xita1_init=xita1; xita2_init=xita2;

     save('Temp/step1.mat','L','p','train_n','J','I_p','lambda1','lambda2','r1','r2','one_p','sigma_0','vector_S1','vector_S2','miu','xita1','xita2','xita1_init','xita1_init');

    % 调用优化函数并传递所需的参数
    [w_opt, f_val] = optimize_weights(xita1_init, xita2_init);
    
   
   w_1 = w_opt(1:r1); w_2 = w_opt(r1+1:end);  


        %%
        % 构建协方差矩阵
     Q1=zeros(p,p);Q2=zeros(p,p);
     for i_2=1:r1
         Q1=Q1+w_1(i_2)*vector_S1(:,i_2)*vector_S1(:,i_2)';
     end
     for i_3=1:r2
         Q2=Q2+w_2(i_3)*vector_S2(:,i_3)*vector_S2(:,i_3)';
     end 
    C=(1/(sigma_0)^2)*(I_p+Q1+Q2);

     QQ1=zeros(p,p);QQ2=zeros(p,p);
     for i_4=1:r1
         QQ1=QQ1+lambda1(i_4)*vector_S1(:,i_4)*vector_S1(:,i_4)';
     end
     for i_5=1:r2
         QQ2=QQ2+lambda2(i_5)*vector_S2(:,i_5)*vector_S2(:,i_5)';
     end 
    Sc=(sigma_0)^2*(I_p+QQ1+QQ2);



    matrix_C{i_1} = C;  % 替换为你的矩阵生成代码
    matrix_SC{i_1} = Sc;
    matrix_miu{i_1}=miu;
    matrix_lambda_large{i_1}=lambda_large;
    matrix_lambda_small{i_1}=lambda_small;
    matrix_Sp{i_1} = Sp;


end

% 保存为MAT文件

 save('Temp\matrix_C.mat','matrix_C');
 save('Temp\matrix_SC.mat','matrix_SC');
 save('Temp\matrix_miu.mat','matrix_miu');
 save('Temp\matrix_SD.mat','matrix_Sp')
 save('Temp\matrix_lambda_large.mat','matrix_lambda_large');
 save('Temp\matrix_lambda_smalle.mat','matrix_lambda_small');
