function [V1,p,lambda_large, lambda_small,k_large,k_small,coeff] = SQM(SW,SZ2,SZ3)
%% 估计总体协方差矩阵的大小Spiked特征值
% SW:是一个协方差矩阵
% X:是样本。size(X)返回n m N，n表示二维样本的维数，N表示二维样本的数量

[m,N]=deal(SZ2,SZ3) ;
S = (SW + SW') ./ 2; % 保证协方差矩阵是对称阵
[V1,D1] = Eig_Update(S); % 奇异值分解，实对称正定矩阵的奇异值分解与谱分解结果相同
Num = N * m; p = length(D1); alpha=0.05;

%% 估计大小Sipked特征根的数量和sigma的平方，k_large0，k_small0，coeff
if p>Num
    D1 = D1(1:Num);
    L = [];

    parfor i =1:500
        simga = diag(ones(1,p));
        x1 = mvnrnd(zeros(1,p), simga, Num);
        [~,l1] = eig(x1*x1'); l1 = sort(diag(l1/Num), 'descend');L = [L;l1'];
    end

    l1 = (mean(L));

    k=floor(min(p,Num)*alpha):floor(min(p,Num)*(1-alpha));
    y_lab = D1(:,k);
    x_lab = l1(:,k);

    [coeff,~,~,~,stats] = regress(y_lab',x_lab');

    % r2=stats(1);
    loss = sqrt(stats(4));

    k_large0 = sum(D1>(l1*coeff + 2*loss));

    k_small0 =0;

    % function [k_large, k_small, r2, coeff] = SQM_pSmall(alpha, p, n, l)
else
    L = [];
    parfor i =1:500
        simga = diag(ones(1,p));
        x1 = mvnrnd(zeros(1,p), simga, Num);

        [~,l1] = eig(x1'*x1); l1 = sort(diag(l1/Num), 'descend');L = [L;l1'];

    end
    l1 = (mean(L));
    k=floor(min(p,Num)*alpha):floor(min(p,Num)*(1-alpha));
    if min(k)<1
        k = k(2:end);
    end
    y_lab = D1(:,k);
    x_lab = l1(:,k);

    [coeff,~,~,~,stats] = regress(y_lab',x_lab');
    % r2=stats(1);
    loss = sqrt(stats(4));

    k_large0 = sum(D1>(l1*coeff + 2*loss));

    l_low0 = ones(1,length(l1))*coeff * (1-sqrt(p/Num))^2;
    l_low1 = (l1*coeff - 2*loss);
    l_low = [l_low0;l_low1];

    [max_a,~] = max(l_low, [],1);
    k_small0 = sum(D1<max_a);

end


%% 估计k_large k_small和其对应的lambda_large, lambda_small

c = p/Num;
% 估计 Large Spiked eigenvalues
if k_large0 ~= 0
    for i_kl = 1: k_large0
        s_large =D1(i_kl);
        temp_ = (s_large/coeff + 1-c).^2-4.*s_large./coeff;
        k_large = k_large0;
        if temp_<0
            k_large = i_kl-1;
            break;
        end
    end
    s_large = D1(1:k_large);
    lambda_large = (s_large/coeff + 1-c + sqrt( (s_large/coeff + 1-c).^2-4.*s_large./coeff))/2-1;
else
    lambda_large = []; k_large =0;
end

% 估计 Small Spiked eigenvalues
if k_small0 ~= 0
    k_small = k_small0;
    for i_ks = 1:k_small0
        s_small = D1(length(D1)-i_ks+1);
        temp_ = (s_small/coeff + 1-c).^2-4.*s_small./coeff;
        if temp_<0
            k_small = i_ks-1;
            break
        end
    end
    s_small = D1((p-k_small+1):p);%l(1:k_small);
    lambda_small = (s_small/coeff + 1-c - (sqrt((s_small/coeff + 1-c).^2-4.*s_small./coeff) ) )/2-1;
else
    lambda_small = []; k_small = 0;
end
% 识别的特征根数量
% fprintf("number of lambda_large is %d \n",k_large);
% fprintf("number of lambda_small is %d \n",k_small);

end
