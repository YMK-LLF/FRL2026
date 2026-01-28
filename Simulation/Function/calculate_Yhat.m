function [Yhat_1, Yhat_2, Yhat_3] = calculate_Yhat(lambda1, lambda2, w_1, w_2, a1, a2, b1, b2, c1, c2, r1, r2, sigma_0, one_p, miu, p)


    % 初始化累加器
    Y_hat1 = 0;
    Y_hat2 = 0;
    Y_hat3 = 0;

    % 计算大spike部分
    for i = 1:r1
        Y_hat1 = Y_hat1 + lambda1(i) * b1(i) ...
                        + 2 * w_1(i) * (1 + lambda1(i)) * a1(i) * b1(i) ...
                        + w_1(i)^2 * (1 + lambda1(i) * a1(i)) * a1(i) * b1(i);

        Y_hat2 = Y_hat2 + lambda1(i) * sqrt(b1(i) * c1(i)) ...
                        + 2 * w_1(i) * (1 + lambda1(i)) * a1(i) * sqrt(b1(i) * c1(i)) ...
                        + w_1(i)^2 * (1 + lambda1(i) * a1(i)) * a1(i) * sqrt(b1(i) * c1(i));

        Y_hat3 = Y_hat3 + lambda1(i) * c1(i) ...
                        + 2 * w_1(i) * (1 + lambda1(i)) * a1(i) * c1(i) ...
                        + w_1(i)^2 * (1 + lambda1(i) * a1(i)) * a1(i) * c1(i);
    end

    % 计算小spike部分
    for i = 1:r2
        Y_hat1 = Y_hat1 + lambda2(i) * b2(i) ...
                        + 2 * w_2(i) * (1 + lambda2(i)) * a2(i) * b2(i) ...
                        + w_2(i)^2 * (1 + lambda2(i) * a2(i)) * a2(i) * b2(i);

        Y_hat2 = Y_hat2 + lambda2(i) * sqrt(b2(i) * c2(i)) ...
                        + 2 * w_2(i) * (1 + lambda2(i)) * a2(i) * sqrt(b2(i) * c2(i)) ...
                        + w_2(i)^2 * (1 + lambda2(i) * a2(i)) * a2(i) * sqrt(b2(i) * c2(i));

        Y_hat3 = Y_hat3 + lambda2(i) * c2(i) ...
                        + 2 * w_2(i) * (1 + lambda2(i)) * a2(i) * c2(i) ...
                        + w_2(i)^2 * (1 + lambda2(i) * a2(i)) * a2(i) * c2(i);
    end

    % 最终计算 Yhat_1, Yhat_2, Yhat_3
    Yhat_1 = 1 / (sigma_0^2) * (1 + Y_hat1);
    Yhat_2 = 1 / (sigma_0^2) * (one_p' * miu / (sqrt(p) * norm(miu)) + Y_hat2);
    Yhat_3 = 1 / (sigma_0^2) * (1 + Y_hat3);
end
