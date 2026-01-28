function [Lshrink1, NLshrink1, SPM1, SD1,SC1] = calculate_mean(test_sample, LWEUP, NLWEUP, SPMWEUP, SDWEUP,SCWEUP, L)

    % 通用计算函数
    calcmean = @(test_sample, W) mean(test_sample * W)

    % 各项计算

    Lshrink1 = calcmean(test_sample, LWEUP);
    NLshrink1 = calcmean(test_sample, NLWEUP);
    SPM1 = calcmean(test_sample, SPMWEUP);
    SD1 = calcmean(test_sample, SDWEUP);
    SC1 = calcmean(test_sample, SCWEUP);

end
