function [ Lshrink1, NLshrink1, SPM1, SD1,SC1] = calculate_EUP(test_sample, LWEUP, NLWEUP, SPMWEUP, SDWEUP,SCWEUP, L)

    % 通用计算函数
    calcEUP = @(test_sample, W) mean(test_sample * W) - L / 2 * var(test_sample * W);

    % 各项计算

    Lshrink1 = calcEUP(test_sample, LWEUP);
    NLshrink1 = calcEUP(test_sample, NLWEUP);
    SPM1 = calcEUP(test_sample, SPMWEUP);
    SD1 = calcEUP(test_sample, SDWEUP);
    SC1 = calcEUP(test_sample, SCWEUP);

end
