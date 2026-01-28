function [b1, b2] = calculate_bj(vector_1, vector_2, p)
    one_p=ones(p,1);
    % 计算大spike对应的bj
    b1 = (one_p' / sqrt(p) * vector_1).^2;

    % 计算小spike对应的bj
    b2 = (one_p' / sqrt(p) * vector_2).^2;
end