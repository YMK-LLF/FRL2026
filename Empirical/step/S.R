
library(nlshrink)
library(R.matlab)

# 读取 CSV 文件
data <- read.csv("C:/Users/data.csv", header = TRUE, stringsAsFactors = FALSE)

matrix_SGL <- list()
matrix_SLS <- list()
matrix_SNL <- list()
train_n=125 #样本量需根据不同数据变化
for (i_1 in 1:500){
  print(i_1)
  
  
  train_sample <- data[i_1:(train_n-1 + i_1), ]
  X <- data.frame(train_sample)
  
  
  # shrink
  X <- as.matrix(X)  # 将数据框转换为矩阵
  SLS <- linshrink_cov(X)
  SNL<- nlshrink_cov(X)
  
  

  matrix_SLS[[paste0("SLS_", i_1)]] <- SLS
  matrix_SNL[[paste0("SNL_", i_1)]] <- SNL
  
  
}


# 假设 .RData 中包含对象 x, y, z，将它们写入 MAT 文件
writeMat("C:/Users/S.mat", matrix_SLS=matrix_SLS,matrix_SNL=matrix_SNL)
