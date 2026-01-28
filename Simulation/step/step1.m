

 %   J = 0.5;

    I_p=diag(ones(1,p));
    L=5;
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
    mu=1/sqrt(p)*ones(p,1);
    one_p=ones(p,1);
    sigma_0 = 1;
    miu=1/sqrt(p)*ones(p,1);
    save('Temp/step1.mat','L','p','train_n','J','I_p','mu','Sigma','lambda1','lambda2','r1','r2','one_p','sigma_0','vector_1','vector_2','miu');
