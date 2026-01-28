function  [a1,a2]= calculate_aj(lambda1, lambda2, J)
 % 计算大spike对应的aj
    a1 = (lambda1.^2 - J) ./ (lambda1 .* (lambda1 + J));
     % 计算小spike对应的aj
    a2 = (lambda2.^2 - J) ./ (lambda2 .* (lambda2 + J));
end