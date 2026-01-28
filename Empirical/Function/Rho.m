function [eta,coeff,k,SW] = Rho(Sw,NormType,SZ1,SZ2,SZ3,varargin)
%RHO 此处显示有关此函数的摘要
%  当输入参数大于5个时，函数仅使用大SPiked特征根
%  假设样本是二维数据，共有 N 个样本，且是 n*m 的。将其看作由m个样本组成的样本矩阵
%  因此，n 表示数据的维数，m表示样本的数量（一个二维样本中包含的一维样本数）
%  J：样本维容比 n/N



[n,m,N] = deal(SZ1,SZ2,SZ3); J = n/N;
Num = N * m;
S = (Sw + Sw') ./ 2 ;

[U1,p,l_large, l_small,k_large,k_small,coeff]= SQM(S,m,N);
k = [k_large,k_small];
% 当且仅当，样本维数p 小于 样本数量Num 时，考虑 小Spiked特征值和对应特征向量的影响
if (p<Num) && (k_small ~= 0) && (nargin <= 5)
    l = [l_large, l_small]' + 1;
    Q = [U1(:,1:k_large),U1(:,(end-k_small+1):end)];
    fprintf("number of lambda_large is %d \n",k_large);
    fprintf("number of lambda_small is %d \n",k_small);
else
    l = l_large' + 1;
    Q =U1(:,1:k_large);
    fprintf("number of lambda_large is %d \n",k_large);
    fprintf("number of lambda_small is %d \n",0);
end


% square of the angle between the eigenvector v of Sigma and the eigenvector u of S
c = (m - J/(l.^2)) ./ (m + J/l); c = reshape(c,[],1);
s = 1 - c; s = reshape(s,[],1);

% 计算eta
switch upper(NormType)
    case {'O' ,1}
        eta = l;
    case {'F',2}
        eta = l.*c + s;
    case {'ST',3}
        eta = l./(c + l.*s);
    case {'N',4}
        eta = zeros(size(l));

        eta(1:k_large) = max( 1 + (l(1:k_large) - 1).*(1 - 2 .* s(1:k_large)),ones([k_large,1]));
        if k_large ~= length(l)
            eta((end-k_small+1):end) = min( 1 + (l((end-k_small+1):end) - 1).*(1 - 2 .* s((end-k_small+1):end)),ones([k_small,1]));
        end

end
eta = eta - 1;
%% Spectral Correction for S

Sp = diag(ones(1,n)) + Q * diag(eta) * Q';
SW = coeff .* (Sp + Sp')./2;
end

