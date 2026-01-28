function [c1, c2] = calculate_cj(miu, vector_1, vector_2)
    % 计算 c1 和 c2 的函数
    % 输入：
    %   mu       - 均值向量
    %   vector_1 - 大 spike 特征向量矩阵（列向量）
    %   vector_2 - 小 spike 特征向量矩阵（列向量）
    % 输出：
    %   c1       - 对应大 spike 的 c 值
    %   c2       - 对应小 spike 的 c 值

    
    % 计算 c1 和 c2
    c1 = (miu' / norm(miu) * vector_1).^2;
    c2 = (miu' / norm(miu)* vector_2).^2;
end
