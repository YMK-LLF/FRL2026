function [w_opt, f_val] = gradient_ascent_max(f, xita1, xita2, df_dw, w_init, w)
    % 获取变量的维度
    r1 = length(xita1);
    r2 = length(xita2);
    
    % 使用传入的初始值 w_init
    w_val = w_init;  % 初始值由外部给定
    
    % 设定梯度上升法参数
    alpha = 0.01;     % 学习率
    tolerance = 1e-5; % 容忍度
    max_iter = 1000;     % 最大迭代次数
    epsilon = 1e-3;   % 用于约束的偏移量

    % 确保 w_val 是列向量
    w_val = w_val(:);
    
    % 梯度上升过程
    for k = 1:max_iter
        k
        % 创建替换对：将符号变量 w 替换为数值 w_val
        subs_pairs = [w(:), num2cell(w_val)];
        
        % 计算梯度值
        grad_vals = double(subs(df_dw, subs_pairs(:,1), subs_pairs(:,2)));
        
        % 确保梯度是列向量
        grad_vals = grad_vals(:);

        % 更新参数
        w_val_new = w_val + alpha * grad_vals;

        % 施加约束
        w_val_new(1:r1) = max(w_val_new(1:r1), -1 + epsilon); % xita1 > -1
        w_val_new(r1+1:end) = max(w_val_new(r1+1:end), 0 + epsilon); % xita2 > 0

        % 检查收敛
        if norm(grad_vals, 'fro') < tolerance
            fprintf('在 %d 次迭代后收敛。\n', k);
            break;
        end

        % 更新权重值
        w_val = w_val_new;
    end

    % 返回最优权重和目标函数值
    w_opt = w_val;  
    f_val = double(subs(f, w(:), num2cell(w_opt)));  
end