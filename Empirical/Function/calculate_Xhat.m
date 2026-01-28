function [Xhat_1, Xhat_2, Xhat_3] = calculate_Xhat( w_1, w_2, a1, a2, b1, b2, c1, c2, r1, r2, sigma_0, one_p, miu, p)
        


    % 初始化累加器
    X_hat1 = 0;
    X_hat2 = 0;
    X_hat3 = 0;

    % 计算大spike部分
    for i = 1:r1
        X_hat1 = X_hat1 + w_1(i) * a1(i) * b1(i);
        X_hat2 = X_hat2 + w_1(i) * a1(i) * sqrt(b1(i) * c1(i));
        X_hat3 = X_hat3 + w_1(i) * a1(i) * c1(i);
    end

    % 计算小spike部分
    for i = 1:r2
        X_hat1 = X_hat1 + w_2(i) * a2(i) * b2(i);
        X_hat2 = X_hat2 + w_2(i) * a2(i) * sqrt(b2(i) * c2(i));
        X_hat3 = X_hat3 + w_2(i) * a2(i) * c2(i);
    end

    % 最终计算 Xhat_1, Xhat_2, Xhat_3
    Xhat_1 = 1 / (sigma_0.^2) * (1 + X_hat1);
    Xhat_2 = 1 / (sigma_0.^2) * (one_p' * miu / (sqrt(p) * norm(miu)) + X_hat2);
    Xhat_3 = 1 / (sigma_0.^2) * (1 + X_hat3);
end
