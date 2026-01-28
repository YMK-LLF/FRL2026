clc;
close all;
clear all;

data=readtable('data.csv'); %读取数据
p=100; % 维度
train_n=125;  %样本量需根据不同数据变化
J=p/train_n;  %维容比
one_p=ones(p,1); %p维1向量
I_p= diag(ones(p, 1));  % 将对角线填充为1
L=5;

% 预分配cell数组
matrix_large=zeros(1, 500);
matrix_small=zeros(1, 500);

for i_1 = 1:500
    i_1
    train_sample=data(i_1:(train_n-1 + i_1), :);
    train_sample= table2array(train_sample);
    Sp=cov(train_sample);
    miu=mean(train_sample);
    miu=miu';
    [eta,coeff,k,SW] = Rho(Sp,1,p,1,train_n);

    large_number = k(1); % 提取第一个
    small_number = k(2); 


    matrix_large(i_1)=large_number;
    matrix_small(i_1)=small_number;
end
   