function [Lshrink2, NLshrink2, SPM2, SD2,SC2] = calculate_var(test_sample, LWEUP, NLWEUP, SPMWEUP, SDWEUP,SCWEUP, L)

    % 通用计算函数
    calcvar = @(test_sample, W) var(test_sample * W);

    % 各项计算

    Lshrink2 = calcvar(test_sample, LWEUP);
    NLshrink2 = calcvar(test_sample, NLWEUP);
    SPM2 = calcvar(test_sample, SPMWEUP);
    SD2 = calcvar(test_sample, SDWEUP);
    SC2 = calcvar(test_sample, SCWEUP);

end
