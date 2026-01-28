function [V, D] = Eig_Update(Sw)
% 使用该函数前请确定输入参数Sw是一个对称矩阵
% 若输入参数是对称矩阵，那么其奇异值分解得到的左右奇异矩阵相等；
% 若得到的特征值是复数，那么其对应的特征向量一定是复数域中的向量
    [singularvector_Sw,singularvalue_Sw,~] = svd(Sw,"econ","vector");  % eigen_vector 列是特征根对应特征向量
    
    % 小于10e-5的特征值，认为其使极小值，遂将其设置为0
    temp1_Sw = singularvalue_Sw > 10e-5; % temp1_Sw 为 逻辑类 因此可以直接用
    V = singularvector_Sw(:,temp1_Sw);
    D = singularvalue_Sw(temp1_Sw)';    
  
end