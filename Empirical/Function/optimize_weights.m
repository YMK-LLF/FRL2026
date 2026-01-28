function [w_opt, f_val] = optimize_weights(xita1_init, xita2_init)
    load('Temp/step1.mat')

    % 定义符号变量
    syms xita1 [r1,1] real  % 改为列向量形式
    syms xita2 [r2,1] real  % 改为列向量形式
    
    % 合并变量（保持列向量形式）
    w = [xita1; xita2];

    % 计算目标函数
    [a1, a2] = calculate_aj(lambda1, lambda2, J);
    [b1, b2] = calculate_bj(vector_S1, vector_S2, p);
    [c1, c2] = calculate_cj(miu, vector_S1, vector_S2);
    [Xhat_1, Xhat_2, Xhat_3] = calculate_Xhat(xita1, xita2, a1, a2, b1, b2, c1, c2, r1, r2, sigma_0, one_p, miu, p);
    [Yhat_1, Yhat_2, Yhat_3] = calculate_Yhat(lambda1, lambda2, xita1, xita2, a1, a2, b1, b2, c1, c2, r1, r2, sigma_0, one_p, miu, p);

    f = norm(miu)*Xhat_2/(sqrt(p)*Xhat_1) + norm(miu)^2*Xhat_3/L ...
        - norm(miu)^2*Xhat_2^2/(L*Xhat_1) - L*Yhat_1/(2*p*Xhat_1^2) ...
        - norm(miu)*Yhat_2/(sqrt(p)*Xhat_1) + norm(miu)*Yhat_1*Xhat_2/(sqrt(p)*Xhat_1^2) ...
        - norm(miu)^2*Yhat_3/(2*L) + norm(miu)^2*Yhat_2*Xhat_2/(L*Xhat_1) ...
        - norm(miu)^2*Yhat_1*Xhat_2^2/(2*L*Xhat_1^2);

    % 计算梯度
    df_dw = gradient(f, w);  

    % 合并初始值（确保是列向量）
    w_init = [xita1_init(:); xita2_init(:)];

    % 调用优化函数
    [w_opt, f_val] = gradient_ascent_max(f, xita1, xita2, df_dw, w_init, w);

    % 输出优化结果
    fprintf('最优权重为: \n');
    disp(w_opt);
    fprintf('目标函数值为 %.4f\n', f_val);
end