
library(nlshrink)

L=5
p=300

J = 0.5
train_n=p/J
train_n=trunc(train_n)
I_p <- diag(1, p, p)  # 创建 p x p 的单位矩阵
# 定义矩阵名字
matrix_names <- c("LWEUP", "LWGMVP", "NLWEUP", "NLWGMVP", "SWEUP", "SWGMVP")

# 使用循环动态创建矩阵并赋予名称
for (name in matrix_names) {
  assign(name, matrix(0, nrow = p, ncol = 500))  # 创建全为零的矩阵
}

mu <- (1 / sqrt(p)) * rep(1, p)


for (i_1 in 1:500){
  print(i_1)
  
  # ***********************生成训练样本集*************************
  train_sample <- data$data.array[[i_1]]
  X <- data.frame(train_sample)
  

  SD=solve(S)

  # shrink
  X <- as.matrix(X)  # 将数据框转换为矩阵
  L_result <- linshrink_cov(X)
  NL_resut<-nlshrink_cov(X)
  
  Lshrink=solve(L_result)
  
  NLshrink=solve(NL_resut)
  #计算w权重
  compute_WEUP_WGMVP <- function(C) {
    
    
    # 验证输入参数
    if (!is.matrix(C)) stop("C 必须是一个矩阵")
    if (!is.numeric(p) || p <= 0) stop("p 必须是正整数")
    if (!is.numeric(mu) || length(mu) != ncol(C)) stop("mu 必须是与 C 列数相等的向量")
    if (!is.numeric(L) || L <= 0) stop("L 必须是正数")
    
    # 计算 WGMVP1
    ones <- matrix(1, nrow = p, ncol = 1)
    numerator <- C %*% ones                            # 分子部分
    denominator <- t(ones) %*% C %*% ones              # 分母部分
    if (is.matrix(denominator) && length(denominator) == 1) {
      denominator <- as.numeric(denominator)
    }
    WGMVP1 <- numerator / denominator
    
    # 计算 w1
    w1 <- (1 / L) * (C %*% mu)
    
    # 计算 w2
    numerator1 <- C %*% ones %*% t(ones) %*% C %*% mu  # 分子部分
    denominator1 <- t(ones) %*% C %*% ones             # 分母部分
    if (is.matrix(denominator1) && length(denominator1) == 1) {
      denominator1 <- as.numeric(denominator1)
    }
    w2 <- (1 / L) * (numerator1 / denominator1)
    
    # 计算 WEUP1
    WEUP1 <- WGMVP1 + w1 - w2
    
    # 返回结果
    return(list(WEUP = WEUP1, WGMVP = WGMVP1))
  }
  
  wLshrink1<- compute_WEUP_WGMVP(Lshrink)
  wNLshrink1<- compute_WEUP_WGMVP(NLshrink)
  wSD1<- compute_WEUP_WGMVP(SD)
  
  # 输出结果

  LWEUP[,i_1]=wLshrink1$WEUP
  LWGMVP[,i_1]=wLshrink1$WGMVP  
  NLWEUP[,i_1]= wNLshrink1$WEUP
  NLWGMVP[,i_1]= wNLshrink1$WGMVP 
  
  SWEUP[,i_1]= wSD1$WEUP
  SWGMVP[,i_1]= wSD1$WGMVP  
  
}



# 假设 .RData 中包含对象 x, y, z，将它们写入 MAT 文件
writeMat("//D1/w/J_0.5/wp300.mat",LWEUP=LWEUP,LWGMVP=LWGMVP,NLWEUP=NLWEUP,NLWGMVP=NLWGMVP,SWEUP=SWEUP,SWGMVP=SWGMVP)


