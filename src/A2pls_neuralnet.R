# 目标
# 这本记事本介绍了利用过去20年来过去12个台风的数据，开发了菲律宾台风造成的建筑物破坏的统计模型。
# 该破坏模型的建立是为了预测全国每个城市房屋受损的比例，使用一套与风险相关的指标，如风速和社会经济指标，
# 如贫困发生率或墙壁建筑质量，作为预测指标。
# 
# 建立了两种破坏模型：第一种是预测某市房屋损坏比例的模型，第二种是预测城市房屋完全受损比例的模型。
# 为了解释的目的，这两种模型都被选择为多元逻辑模型(估计系数与对数-概率相关)，并且因为一个概率模型
# 具有相似的性能。

setwd("E:\\storm-goni-main-ok\\test")
library(magrittr)
library(ggplot2)
library(plm)
library(MASS) # 基于此包进行box-cox转换
library(carData)
library(car)
library(mice)# 缺失值处理工具包，本节内容会使用它来进行多重插补
library(VIM) # 缺失值可视化工具包，本节内容中将会使用到它的aggr()函数
library(tidyverse) # tidyverse是为数据科学设计的R软件包，它包含(ggplot2、dplyr、tidyr、stringr、magrittr、tibble)等一系列热门软件包
#install.packages("tidyverse")
library("neuralnet")



#Objective
chosen_link <- "logit"
if(chosen_link == "logit"){link_fn <- logit}
if(chosen_link == "probit"){link_fn <- qnorm}

include_transform <- TRUE # Set = TRUE if we want to compute cross-validation errors accounting for variable transformation, scaling, VIF and stepAIC.


# Read CSV with complete dataset into R
raw_data <- read.csv(file="data//ALL.csv", header=TRUE, sep=",")
summary(raw_data)
# Remove rows with NA in skilled.Agriculture.Forestry.Fishermen (21 rows):
raw_data <- raw_data[!is.na(raw_data$X..skilled.Agriculture.Forestry.Fishermen), ]

# Remove row with NA in experience.
raw_data <- raw_data[!is.na(raw_data$Experience), ]

# Remove unnecessary columns
raw_data <- raw_data[, -which(colnames(raw_data) %in% c("Completely.damaged..abs..", "Partly.damaged..abs..", "Total.damaged.houses..abs..", "total_damage_houses_0p25weight", "Partly.damaged..rel..", "total_damage_houses_0p25weight_perc", "ratio_comp_part", "Total...of.houses", "Predicted.damage.class..1.5.", "Bicol.region", "Population", "land_area", "Ruggedness","ruggedness_stdev"))]
names(raw_data)
raw_data %>%
  #keep(is.numeric) %>% 
  keep(is.numeric) %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()
# 陈圣群 https://blog.csdn.net/weixin_43595036/article/details/119820737
# 创建直观看一下具体的缺失情况
aggr(raw_data, cex.axis = .6, oma = c(9, 5, 5, 1)) # cex.axis调整轴字体大小，oma调整外边框大小

# 陈圣群创建密度图
raw_data %>%
  keep(is.numeric) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~key, scales = "free") +
  geom_density(color = "#348498", fill = "#8ac6d1") +
  theme_bw()

nums_data_norm1 <- scale(is.numeric(raw_data))
nums_data_norm1 <- raw_data
for (col in names(nums_data_norm1))
{
  if(is.numeric(nums_data_norm1[col])){
    xmin <- min(nums_data_norm1[col])
    xmax <- max(nums_data_norm1[col])
    nums_data_norm1[col] <- (nums_data_norm1[col] - xmin) / (xmax - xmin)
    }
}
head(nums_data_norm1)

nums_data_norm1 %>%
  keep(is.numeric) %>%
  gather() %>%
  ggplot(aes(value)) +
  facet_wrap(~key, scales = "free") +
  geom_density(color = "#348498", fill = "#8ac6d1") +
  theme_bw()


# 提取数值字段
nums <- unlist(lapply(raw_data, is.numeric))
nums_data <- raw_data[, nums]
# 数据变形
nums_data.new <- nums_data %>%
  as.data.frame() %>%
  mutate(Cell = rownames(.)) %>%
  gather(., key = colname, value = "value", -Cell)
# 用ggplot画出箱线图
ggplot(data = nums_data.new, aes(x = colname, y = value)) +
  geom_boxplot(aes(1)) +
  facet_wrap(~colname, scales = "free") +
  theme_grey() +
  labs(title = "Outlier Detection On Numeric Data By Boxplot", x = "Numeric Columns", y = "") +
  theme(legend.position = "top") +
  theme_bw()



## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
# Remove outlier:
raw_data <- raw_data[-which(raw_data$Wind.speed == max(raw_data$Wind.speed)),]

sum(raw_data$Total.damaged.houses..abs.. > raw_data$Total...of.houses)

max(raw_data$Total.damaged.houses..rel..[raw_data$Total.damaged.houses..rel..<100])

# Replace >100 Total.damaged.houses..rel.. with 99.99
raw_data$Total.damaged.houses..rel..[raw_data$Total.damaged.houses..rel..>100] <- 99.99

# Create proportion of houses which are damaged:
raw_data$Total.damaged.houses..rel.. <- raw_data$Total.damaged.houses..rel../100

# Create proportion of houses which are completely damaged:
raw_data$Completely.damaged..rel.. <- raw_data$Completely.damaged..rel../100


# Link-transformed damage proportions 
link_fn_td <- link_fn(raw_data$Total.damaged.houses..rel..)
link_fn_pc <- link_fn(raw_data$Completely.damaged..rel..)

# Round lambda est to the nearest half integer
power_est <- function(bC){l.est <- bC$x[bC$y == max(bC$y)]
return(round(l.est/0.5)*0.5)}

# Set up dataframe with first row for td and second row for pc
var_names <- colnames(raw_data)
var_names <- var_names[!(var_names %in% c("disaster_type", "disaster_name", "pcode", "Completely.damaged..rel..", "Total.damaged.houses..rel..", "Region", "prov"))]
power_df <- matrix(NA, nrow = 2, ncol = length(var_names))
colnames(power_df) <- var_names

for (i in 1:length(var_names)){
  
  var_i <- var_names[i]
  power_df[1, var_i] <- power_est(boxCox(link_fn_td - min(link_fn_td)*1.01 ~ raw_data[, var_i],
                                         plotit = FALSE))
  power_df[2, var_i] <- power_est(boxCox(link_fn_pc - min(link_fn_pc)*1.01 ~ raw_data[, var_i],
                                         plotit = FALSE))
  
}

power_df


td_data <- raw_data
pc_data <- raw_data



# Standardise data:

# Save to apply to Goni 2020 dataset later

cols.mean.td <- colMeans(td_data[, var_names] ) 
cols.sd.td <- apply(td_data[, var_names], 2, sd )
scaled.data.td <- lapply(td_data[, var_names], scale)

cols.mean.pc <- colMeans(td_data[, var_names] ) 
cols.sd.pc <- apply(td_data[, var_names], 2, sd )
scaled.data.pc <- lapply(pc_data[, var_names], scale)

# 数据缩放所有有数字列都是标准化的，平均值为0，标准偏差为1。
# Switch out scaled numerical columns (apart from response variables)
td_data[, var_names] <- as.data.frame(scaled.data.td)
pc_data[, var_names] <- as.data.frame(scaled.data.pc)


td_data %>%
  keep(is.numeric) %>% 
  gather(-Total.damaged.houses..rel.., -Completely.damaged..rel.., key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = link_fn(Total.damaged.houses..rel..))) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

pc_data %>%
  keep(is.numeric) %>% 
  gather(-Total.damaged.houses..rel.., -Completely.damaged..rel.., key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = link_fn(Completely.damaged..rel..))) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

# 探索性数据分析，在这一部分，我们探索数据，以确定可能的关系之间的预测和建筑物损坏的比例。

# Dataset for Proportion Damaged:
# Omit Population, Land.Area, Region and prov.
omit_var_td <- which(colnames(raw_data) %in% c("disaster_type", "disaster_name", "pcode", "Completely.damaged..rel..", "Region", "prov"))
td_data <- td_data[, -omit_var_td]
# head(td_data)
# print(names(td_data))
# Dataset for Proportion of Completely Damaged:
omit_var_pc <- which(colnames(raw_data) %in% c("disaster_type", "disaster_name", "pcode", "Total.damaged.houses..rel..", "Region", "prov"))
pc_data <- pc_data[, -omit_var_pc]
head(pc_data)


# -----------PLS
library(pls)
#下面使用8个成分
# #td_data.plsr <- plsr(Total.damaged.houses..rel.. ~ .,data=td_data,scale=T,validation="CV")
Y = as.matrix(td_data[,1])
X = as.matrix(td_data[,-1])
td_data.plsr8 <- plsr(Y ~ X, 8,scale=T,validation="LOO")
td_data.plsr8$coefficients
Z1_scores = td_data.plsr8$scores
new_td_data = cbind(td_data[,1] ,Z1_scores)
td_data = data.frame(new_td_data)
colnames(td_data) = c("Total.damaged.houses..rel..", "x1","x2","x3","x4","x5","x6","x7","x8")
validationplot(td_data.plsr8)
predplot(td_data.plsr8)
# jack.test(td_data.plsr8,8)
# Write new dataset in csv
# write.csv(td_data, 'data\\td_data1.csv', row.names=T)
# -----
Y = as.matrix(pc_data[,1])
X = as.matrix(pc_data[,-1])
pc_data.plsr8 <- plsr(Y ~ X, 8,scale=T,validation="LOO")
Z2_scores = pc_data.plsr8$scores
new_pc_data = cbind(pc_data[,1] ,Z2_scores)
pc_data = data.frame(new_pc_data)
colnames(pc_data) = c("Completely.damaged..rel..", "y1","y2","y3","y4","y5","y6","y7","y8")
# -----------PLS

set.seed(12345)
# Fit Multi-Logistic Model for td
# td_full <- glm(Total.damaged.houses..rel.. ~ ., data= td_data, family = binomial(link=chosen_link), trace = FALSE)
td_full<-neuralnet(Total.damaged.houses..rel.. ~ ., data= td_data, hidden=10,linear.output=FALSE)


## Check on variable collinearity
temp_model <- td_full
# temp_vif <- vif(temp_model) # Use a threshold of 10 to remove variables.
temp_var <- "dummy"
temp_data <- td_data

# while(max(temp_vif)>10){
#   
#   var_remove <- which(temp_vif == max(temp_vif))
#   temp_var <- append(temp_var, names(var_remove))
#   temp_data <- temp_data[, -which(colnames(temp_data) == names(var_remove))]
#   #temp_model <- (glm(Total.damaged.houses..rel.. ~ ., data= temp_data, family = binomial(link=chosen_link), trace = FALSE))
#   temp_model<-(neuralnet(Total.damaged.houses..rel.. ~ ., data= temp_data, hidden=10,linear.output=FALSE))
#   temp_vif <- vif(temp_model) # Use a threshold of 10 to remove variables.
# }

td_col_to_use <- colnames(temp_data)
temp_vif
temp_var
set.seed(12345)
# Fit Multi-Logistic Model for pc
# pc_full <- glm(Completely.damaged..rel.. ~ ., data= pc_data, family = binomial(link=chosen_link), trace = FALSE)
pc_full<- neuralnet(Completely.damaged..rel.. ~ ., data= pc_data, hidden=10,linear.output=FALSE)

## Check on variable collinearity
temp_model <- pc_full
# temp_vif <- vif(temp_model) # Use a threshold of 10 to remove variables.
temp_var <- "dummy"
temp_data <- pc_data

# while(max(temp_vif)>10){
#   
#   var_remove <- which(temp_vif == max(temp_vif))
#   temp_var <- append(temp_var, names(var_remove))
#   temp_data <- temp_data[, -which(colnames(temp_data) == names(var_remove))]
#   #temp_model <- (glm(Completely.damaged..rel.. ~ ., data= temp_data, family = binomial(link=chosen_link), trace = FALSE))
#   temp_model<- (neuralnet(Completely.damaged..rel.. ~ ., data= temp_data, hidden=10,linear.output=FALSE))
#   temp_vif <- vif(temp_model) # Use a threshold of 10 to remove variables.
# }

pc_col_to_use <- colnames(temp_data)
temp_vif
temp_var

# Subset variables based on VIF:
td_data <- td_data[, colnames(td_data) %in% td_col_to_use]
pc_data <- pc_data[, colnames(pc_data) %in% pc_col_to_use]

# Fit Multi-Logistic Model
# td_full <- glm(Total.damaged.houses..rel.. ~ ., data= td_data, family = binomial(link=chosen_link), trace = FALSE)
td_full<-neuralnet(Total.damaged.houses..rel.. ~ ., data= td_data, hidden=10,linear.output=FALSE)
td_step <- stepAIC(td_full)

# Check the variable significance
summary(td_step)


# Fit Multi-Logistic Model
# pc_full <- glm(Completely.damaged..rel.. ~ ., data= pc_data, family = binomial(link=chosen_link), trace = FALSE)
pc_full<-neuralnet(Completely.damaged..rel.. ~ ., data= pc_data, hidden=10,linear.output=FALSE)
pc_step <- stepAIC(pc_full)
pc_step = pc_full

# Check resulting model
pc_step

# Check the variable significance
summary(pc_step)

# 模型验证
# Cross-validation errors

typhoons <- unique(raw_data$disaster_name)
table(raw_data$disaster_name)

cv_pc <- data.frame("Split" = 1:length(typhoons))
cv_pc$RMSE <- NA
cv_pc$MAPE <- NA

cv_td <- data.frame("Split" = 1:length(typhoons))
cv_td$RMSE <- NA
cv_td$MAPE <- NA
#length(typhoons)
for (i in 1:length(typhoons)){
  #print(i)
  #i = 2
  #train.i <- which(raw_data$disaster_name != typhoons[i])
  train.i <- sample(1:nrow(raw_data), dim(raw_data)[1]*0.8)
  #print(typhoons[i])
  if(include_transform){
    
    # Link-transformed damage proportions 
    link_fn_td.i <- link_fn(raw_data[train.i, ]$Total.damaged.houses..rel..)
    link_fn_pc.i <- link_fn(raw_data[train.i, ]$Completely.damaged..rel..)
    
    # Set up dataframe with first row for td and second row for pc
    power_df.i <- matrix(NA, nrow = 2, ncol = length(var_names))
    colnames(power_df.i) <- var_names
    
    td_data.i <- raw_data
    pc_data.i <- raw_data
    
    td_data.i <- td_data.i[, -omit_var_td]
    pc_data.i <- pc_data.i[, -omit_var_pc]   
    
    for (j in 1:length(var_names)){
  
      var_i <- var_names[j]
      power_df.i[1, var_i] <- power_est(boxCox(link_fn_td.i - min(link_fn_td.i)*1.01 ~ raw_data[train.i, ][, var_i], plotit = FALSE))
      power_df.i[2, var_i] <- power_est(boxCox(link_fn_pc.i - min(link_fn_pc.i)*1.01 ~ raw_data[train.i, ][, var_i], plotit = FALSE))
      
      if(power_df.i[1, var_i]==0){ 
        td_data.i[, var_i] <- log(td_data.i[, var_i])
      }else{td_data.i[, var_i] <- td_data.i[, var_i]^power_df.i[1, var_i]}
      
      
      cols.mean.td.i <- mean(td_data.i[train.i, var_i]) 
      cols.sd.td.i <- sd(td_data.i[train.i, var_i])
      td_data.i[, var_i] <- (td_data.i[, var_i] - cols.mean.td.i)/cols.sd.td.i
      
      if(power_df.i[2, var_i]==0){
        pc_data.i[, var_i] <- log(pc_data.i[, var_i])
      }else{pc_data.i[, var_i] <- pc_data.i[, var_i]^power_df.i[2, var_i]}
      
      cols.mean.pc.i <- mean(pc_data.i[train.i, var_i]) 
      cols.sd.pc.i <- sd(pc_data.i[train.i, var_i])
      pc_data.i[, var_i] <- (pc_data.i[, var_i] - cols.mean.pc.i)/cols.sd.pc.i
      
    }
    
    #-----------PLS
    library(pls)
    #下面使用8个成分
    # #td_data.plsr <- plsr(Total.damaged.houses..rel.. ~ .,data=td_data,scale=T,validation="CV")
    Y = as.matrix(td_data.i[,1])
    X = as.matrix(td_data.i[,-1])
    td_data.i.plsr8 <- plsr(Y ~ X, 8,scale=T,validation="LOO")
    Z1_scores = td_data.i.plsr8$scores
    new_td_data.i = cbind(td_data.i[,1] ,Z1_scores)
    td_data.i = data.frame(new_td_data.i)
    colnames(td_data.i) = c("Total.damaged.houses..rel..", "x1","x2","x3","x4","x5","x6","x7","x8")
    #-----
    Y = as.matrix(pc_data.i[,1])
    X = as.matrix(pc_data.i[,-1])
    pc_data.i.plsr8 <- plsr(Y ~ X, 8,scale=T,validation="LOO")
    Z2_scores = pc_data.i.plsr8$scores
    new_pc_data.i = cbind(pc_data.i[,1] ,Z2_scores)
    pc_data.i = data.frame(new_pc_data.i)
    colnames(pc_data.i) = c("Completely.damaged..rel..", "y1","y2","y3","y4","y5","y6","y7","y8")
    #-----------PLS

    
    # Fit Multi-Logistic Model
    #td_full.i <- glm(Total.damaged.houses..rel.. ~ ., data= td_data.i, family = binomial(link=chosen_link), trace = FALSE)
    td_full.i<-neuralnet(Total.damaged.houses..rel.. ~ ., data= td_data.i, hidden=10,linear.output=FALSE)
    ## Check on variable collinearity
    temp_model <- td_full.i
    # temp_vif <- vif(temp_model) # Use a threshold of 10 to remove variables.
    temp_var <- "dummy"
    temp_data <- td_data.i
    #print(temp_data)
    
    # while(max(temp_vif)>10){
    #   
    #   var_remove <- which(temp_vif == max(temp_vif))
    #   temp_var <- append(temp_var, names(var_remove))
    #   temp_data <- temp_data[, -which(colnames(temp_data) == names(var_remove))]
    #   #temp_model <- (glm(Total.damaged.houses..rel.. ~ ., data= temp_data, family = binomial(link=chosen_link), trace = FALSE))
    #   temp_model<-neuralnet(Total.damaged.houses..rel.. ~ ., data= temp_data, hidden=10,linear.output=FALSE)
    #   temp_vif <- vif(temp_model) # Use a threshold of 10 to remove variables.
    # }
    
    td_col_to_use.i <- colnames(temp_data)
    
    # Fit Multi-Logistic Model
    #pc_full.i <- glm(Completely.damaged..rel.. ~ ., data= pc_data.i, family = binomial(link=chosen_link), trace = FALSE)
    pc_full.i<-neuralnet(Completely.damaged..rel.. ~ ., data= pc_data.i, hidden=10,linear.output=FALSE)
    ## Check on variable collinearity
    temp_model <- pc_full.i
    # temp_vif <- vif(temp_model) # Use a threshold of 10 to remove variables.
    temp_var <- "dummy"
    temp_data <- pc_data.i
    
    # while(max(temp_vif)>10){
    #   
    #   var_remove <- which(temp_vif == max(temp_vif))
    #   temp_var <- append(temp_var, names(var_remove))
    #   temp_data <- temp_data[, -which(colnames(temp_data) == names(var_remove))]
    #   #temp_model <- (glm(Completely.damaged..rel.. ~ ., data= temp_data, family = binomial(link=chosen_link), trace = FALSE))
    #   temp_model<-(neuralnet(Completely.damaged..rel.. ~ ., data= temp_data, hidden=10,linear.output=FALSE))
    #   temp_vif <- vif(temp_model) # Use a threshold of 10 to remove variables.
    # }
    
    pc_col_to_use.i <- colnames(temp_data)
    
    td_data.i <- td_data.i[, colnames(td_data.i) %in% td_col_to_use.i]
    pc_data.i <- pc_data.i[, colnames(pc_data.i) %in% pc_col_to_use.i]
    
    td.training.data.i <- td_data.i[train.i, ]
    td.test.data.i <- td_data.i[-train.i, ]
    
    pc.training.data.i <- pc_data.i[train.i, ]
    pc.test.data.i <- pc_data.i[-train.i, ]
    
    #pc_full.i <- glm(Completely.damaged..rel.. ~ ., data= pc.training.data.i, family = binomial(link=chosen_link), trace = FALSE)
    #pc_model.i <-  stepAIC(pc_full.i, trace = 0)
    pc_full.i<-neuralnet(Completely.damaged..rel.. ~ ., data= pc.training.data.i, hidden=10,linear.output=FALSE)
    pc_model.i <-  pc_full.i
    #td_full.i <- glm(Total.damaged.houses..rel.. ~ ., data= td.training.data.i, family = binomial(link=chosen_link), trace = FALSE)
    td_full.i<-neuralnet(Total.damaged.houses..rel.. ~ ., data= td.training.data.i, hidden=10,linear.output=FALSE)
    td_model.i <- td_full.i
    #td_model.i <-  stepAIC(td_full.i, trace = 0)
    #print(td_model.i)
  }else{
    
    pc.training.data.i <- pc_data[train.i, ]
    pc.test.data.i <- pc_data[-train.i, ]
    
    pc_model <- pc_step$formula
    #pc_model.i <- glm(pc_model, data= pc.training.data.i, family = binomial(link=chosen_link), trace = FALSE)
    pc_model.i<-neuralnet(pc_model, data= pc.training.data.i, hidden=10,linear.output=FALSE)
    
    td.training.data.i <- td_data[train.i, ]
    td.test.data.i <- td_data[-train.i, ]
    
    td_model <- td_step$formula
    #td_model.i <- glm(td_model, data= td.training.data.i, family = binomial(link=chosen_link), trace = FALSE)
    td_model.i<-neuralnet(td_model, data= td.training.data.i, hidden=10,linear.output=FALSE)
  }
  
  pc.test.pred.i <- predict(pc_model.i, newdata = pc.test.data.i, type = "response")
  pc.test.res.i <- pc.test.pred.i - pc.test.data.i$Completely.damaged..rel..
  pc.test.rmse.i <- sqrt(mean(pc.test.res.i^2))
  pc.test.mape.i <- mean(abs(pc.test.res.i/pc.test.data.i$Completely.damaged..rel..))
  #cv_pc[i, c("RMSE", "MAPE")] <- c(pc.test.rmse.i, pc.test.mape.i)
  
  pc.test.meangap.i <- mean(pc.test.data.i$Completely.damaged..rel..)-pc.test.data.i$Completely.damaged..rel..
  RR <- 1-sum(pc.test.res.i^2)/sum(pc.test.meangap.i^2)
  cv_pc[i, c("RMSE","MAPE", "RR")] <- c(pc.test.rmse.i, pc.test.mape.i, RR)
  

  # pcd1 <- paste(pc.test.pred.i, pc.test.data.i$Completely.damaged..rel..)
  
  # write.csv(pc.test.pred.i, 'data\\predictdata.csv', row.names=T)
  # write.csv(pc.test.data.i$Completely.damaged..rel.., 'data\\realdata.csv', row.names=T)
  
  # If the estimate of houses damaged is less than the estimate of houses completely damaged, set it to the latter.
  td.test.pred.i <- predict(td_model.i, newdata = td.test.data.i, type = "response")
  td.test.pred.i[td.test.pred.i<pc.test.pred.i] <- pc.test.pred.i
  td.test.res.i <- td.test.pred.i - td.test.data.i$Total.damaged.houses..rel..
  td.test.rmse.i <- sqrt(mean(td.test.res.i^2))
  td.test.mape.i <- mean(abs(td.test.res.i/td.test.data.i$Total.damaged.houses..rel..))
  td.test.meangap.i = mean(td.test.data.i$Total.damaged.houses..rel..)-td.test.data.i$Total.damaged.houses..rel..
  RR2 = 1-sum(td.test.res.i^2)/sum(td.test.meangap.i^2)
  cv_td[i, c("RMSE", "MAPE","RR2")] <- c(td.test.rmse.i, td.test.mape.i, RR2)
  
}

colMeans(cv_pc[, c("RMSE", "MAPE","RR")])
cv_pc

colMeans(cv_td[, c("RMSE", "MAPE","RR2")])
cv_td
cv_pc[length(typhoons)+1, c("RMSE", "MAPE","RR")]<- colMeans(cv_pc[, c("RMSE", "MAPE","RR")])
cv_td[length(typhoons)+1, c("RMSE", "MAPE","RR2")]<- colMeans(cv_td[, c("RMSE", "MAPE","RR2")])
# Write new dataset in csv
write.csv(cv_pc, 'data\\A2pls_neuralnet_pc.csv', row.names=F)
write.csv(cv_td, 'data\\A2pls_neuralnet_td.csv', row.names=F)
# #残差分析
# pc.pred <- predict(pc_step, type = "response")
# pc.res <- pc.pred - pc_data$Completely.damaged..rel..
# 
# hist(pc.res)
# plot(link_fn(pc_data$Completely.damaged..rel..), link_fn(pc.pred), asp = 1)
# abline(a = 0, b = 1, add = TRUE)
# 
# td.pred <- predict(td_step, type = "response")
# # If the estimate of houses damaged is less than the estimate of houses completely damaged, set it to the latter.
# td.pred[td.pred<pc.pred] <- pc.pred
# 
# td.res <- td.pred - td_data$Total.damaged.houses..rel..
# 
# hist(td.res)
# 
# plot(link_fn(td_data$Total.damaged.houses..rel..), link_fn(td.pred), asp = 1)
# abline(a = 0, b = 1, add = TRUE)

#Prediction for Goni Typhoon

# Read CSV with new dataset
new_data <- read.csv(file="data//All_Goni_pred.csv", header=TRUE, sep=",")

# Choose which wind speed data to use
new_data$Wind.speed <- new_data$Wind.speed.max
# new_data$Wind.speed <- new_data$Wind.speed.mean

# Set a minimum distance from typhoon and distance from first impact
min_dist <- 0.01*min(new_data$distance_first_impact[new_data$distance_first_impact>0], new_data$Distance.to.typhoon[new_data$Distance.to.typhoon>0])
new_data$distance_first_impact[new_data$distance_first_impact==0] <- min_dist
new_data$Distance.to.typhoon[new_data$Distance.to.typhoon==0] <- min_dist

head(new_data)

#数据转换和缩放
new_td_data <- new_data
new_pc_data <- new_data

for (i in 1:length(var_names)){
  
  if(power_df[1, var_names[i]]==0){
    new_td_data[, var_names[i]] <- (log(new_td_data[, var_names[i]]) - cols.mean.td[i])/cols.sd.td[i]
  }else{
    new_td_data[, var_names[i]] <- ((new_td_data[, var_names[i]]^power_df[1, var_names[i]]) - cols.mean.td[i])/cols.sd.td[i]
  }
  
  if(power_df[2, var_names[i]]==0){
    new_pc_data[, var_names[i]] <- (log(new_pc_data[, var_names[i]]) - cols.mean.pc[i])/cols.sd.pc[i]
  }else{
    new_pc_data[, var_names[i]] <- (new_pc_data[, var_names[i]]^power_df[2, var_names[i]] - cols.mean.pc[i])/cols.sd.pc[i]
  }
  
}

td_step

# Predict proportion of all damaged houses and standard errors
pred.td_prop_bldgs = predict(td_step, newdata = new_td_data, type = "response", se.fit = TRUE)
# Predict number of all damaged houses
pred.td_num_bldgs = pred.td_prop_bldgs$fit * new_data$Total...of.houses
# Add predicted values and standard errors to data set
new_data$Total.damaged.houses..rel.. = pred.td_prop_bldgs$fit
new_data$Total.damaged.houses..abs.. = pred.td_num_bldgs
new_data$Total.damaged.houses..rel..se = pred.td_prop_bldgs$se.fit
new_data$Total.damaged.houses..abs..se = pred.td_prop_bldgs$se.fit * new_data$Total...of.houses

pc_step



# Predict proportion of all damaged houses
pred.pc_prop_bldgs = predict(pc_step, newdata =new_pc_data, type = "response", se.fit = TRUE)
# Predict number of all damaged houses
pred.pc_num_bldgs = pred.pc_prop_bldgs$fit * new_data$Total...of.houses
# Add predicted values and standard errors to data set
new_data$Completely.damaged..rel.. = pred.pc_prop_bldgs$fit
new_data$Completely.damaged..abs.. = pred.pc_num_bldgs
new_data$Completely.damaged..rel..se = pred.pc_prop_bldgs$se.fit
new_data$Completely.damaged..abs..se = pred.pc_prop_bldgs$se.fit * new_data$Total...of.houses


head(new_data)

# Write new dataset in csv
write.csv(new_data, 'data\\All_Goni_pred.csv', row.names=F)


#读取并绘制Shapefile
#install.packages("sp")
library(sp)

library(rgdal)

library(sf)

library(ggplot2)
dsn <- system.file("vectors", package = "rgdal")
# Read shapefiles
ph_mun = st_read('data\\Shapefiles\\PH_municipality.shp')

goni_trace = st_read('data\\Shapefiles\\goni-path.shp')


# Create new shapefile
ph_mun_new = ph_mun

tdprop.list = replicate(length(ph_mun$Mun_Code),NA)
tdnum.list = replicate(length(ph_mun$Mun_Code),NA)
pcprop.list = replicate(length(ph_mun$Mun_Code),NA)
pcnum.list = replicate(length(ph_mun$Mun_Code),NA)
tdprop.se.list = replicate(length(ph_mun$Mun_Code),NA)
tdnum.se.list = replicate(length(ph_mun$Mun_Code),NA)
pcprop.se.list = replicate(length(ph_mun$Mun_Code),NA)
pcnum.se.list = replicate(length(ph_mun$Mun_Code),NA)
windspeed.list = replicate(length(ph_mun$Mun_Code),NA)
dist_first_impact.list = replicate(length(ph_mun$Mun_Code),NA)
Slope.list  = replicate(length(ph_mun$Mun_Code),NA)
Elevation.list  = replicate(length(ph_mun$Mun_Code),NA)
Population.density.list  = replicate(length(ph_mun$Mun_Code),NA)
X..strong.roof.type.list = replicate(length(ph_mun$Mun_Code),NA)
Experience.list  = replicate(length(ph_mun$Mun_Code),NA)

for (i in 1:length(new_data$pcode)){
  pcode = new_data$pcode[i]
  windspeed = new_data$Wind.speed[i]
  dist_first_impact = new_data$distance_first_impact[i]
  tdprop = new_data$Total.damaged.houses..rel..[i]
  tdnum = new_data$Total.damaged.houses..abs..[i]
  pcprop = new_data$Completely.damaged..rel..[i]
  pcnum = new_data$Completely.damaged..abs..[i]
  tdprop.se = new_data$Total.damaged.houses..rel..se[i]
  tdnum.se = new_data$Total.damaged.houses..abs..se[i]
  pcprop.se = new_data$Completely.damaged..rel..se[i]
  pcnum.se = new_data$Completely.damaged..abs..se[i]
  Slope  = new_data$Slope[i]
  Elevation  = new_data$Elevation[i]
  Population.density  = new_data$Population.density[i]
  X..strong.roof.type = new_data$X..strong.roof.type[i]
  Experience  = new_data$Experience[i]
  
  
  # Look for pcode in shapefile attributes
  index = which(ph_mun$Mun_Code == as.character(pcode))
  tdprop.list[index] = tdprop
  tdnum.list[index] = tdnum
  pcprop.list[index] = pcprop
  pcnum.list[index] = pcnum
  windspeed.list[index] = windspeed
  dist_first_impact.list[index] = dist_first_impact
  tdprop.se.list[index] = tdprop.se
  tdnum.se.list[index] = tdnum.se
  pcprop.se.list[index] = pcprop.se
  pcnum.se.list[index] = pcnum.se
  Slope.list[index]  = Slope
  Elevation.list[index]  = Elevation
  Population.density.list[index]  = Population.density
  X..strong.roof.type.list[index] = X..strong.roof.type
  Experience.list[index]  = Experience
}
# Save new fields in output shapefile
ph_mun_new$Wind_Speed = windspeed.list
ph_mun_new$dist_first_impact = dist_first_impact.list
ph_mun_new$tdprop = tdprop.list
ph_mun_new$tdnum = tdnum.list
ph_mun_new$pcprop = pcprop.list
ph_mun_new$pcnum = pcnum.list
ph_mun_new$tdprop.se = tdprop.se.list
ph_mun_new$tdnum.se = tdnum.se.list
ph_mun_new$pcprop.se = pcprop.se.list
ph_mun_new$pcnum.se = pcnum.se.list
ph_mun_new$Slope = Slope.list 
ph_mun_new$Elevation = Elevation.list 
ph_mun_new$Population.density = Population.density.list 
ph_mun_new$X..strong.roof.type = X..strong.roof.type.list 
ph_mun_new$Experience = Experience.list 

# st_write(ph_mun_new, 'data\\Shapefiles\\PH_municipality_results.shp')

cap.df <- data.frame(x = 126.8, y = 15.75, text = "Developed by the Disaster Analytics for Society Lab (DASL) at the \n Earth Observatory of Singapore and Nanyang Technological University.")

# Plot tdprop
ggplot() + theme(
  panel.background = element_rect(fill = "lightblue",
                                  colour = "lightblue",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.2, linetype = 'dashed',
                                  colour = "white")
) +
  geom_sf(data = ph_mun_new, aes(fill=tdprop), size = 0.1, color = "black") + 
  scale_fill_fermenter(breaks = c(0,0.2,0.4,0.6,0.8,1), type = "seq", palette = 7, direction=1,
                       na.value='grey50') + 
  geom_sf(data = goni_trace, size = 0.2, color = "red", fill = "cyan1") + 
  coord_sf(xlim=c(118,127), ylim=c(11,16), expand=F) + labs(fill = "Damage rate", x = "", y = "") + geom_text(data = cap.df, aes(x = x, y = y, label = text), size = 2.5, hjust = 1) +
  theme(legend.position =  c(0.025, 0.15), legend.justification = c(0, 0), legend.title = element_text(size = 7.5), legend.text = element_text(size = 7.5))

ggsave("graphics/goni-damagerate.png", dpi = 300, height = 6, width = 10, units = "in")

# Plot tdprop standard errors:

ggplot() + theme(
  panel.background = element_rect(fill = "lightblue",
                                  colour = "lightblue",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.2, linetype = 'dashed',
                                  colour = "white")
) +
  geom_sf(data = ph_mun_new, aes(fill=tdprop.se), size = 0.1, color = "black") + 
  scale_fill_gradient(low="white", high="red") + 
  geom_sf(data = goni_trace, size = 0.2, color = "red", fill = "cyan1") + 
  coord_sf(xlim=c(118,127), ylim=c(11,16), expand=F) + labs(fill = "Damage rate \n standard error", x = "", y = "") + geom_text(data = cap.df, aes(x = x, y = y, label = text), size = 2.5, hjust = 1) +
  theme(legend.position =  c(0.025, 0.15), legend.justification = c(0, 0), legend.title = element_text(size = 7.5), legend.text = element_text(size = 7.5))


# Plot tdnum
ggplot() + theme(
  panel.background = element_rect(fill = "lightblue",
                                  colour = "lightblue",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.2, linetype = 'dashed',
                                  colour = "white")
) +
  geom_sf(data = ph_mun_new, aes(fill=tdnum), size = 0.1, color = "black") + 
  scale_fill_fermenter(breaks = c(0,0.2,0.4,0.6,0.8)*10000, type = "seq", palette = 7, direction=1,
                       na.value='grey50') + 
  geom_sf(data = goni_trace, size = 0.2, color = "red", fill = "cyan1") + 
  coord_sf(xlim=c(118,127), ylim=c(11,16), expand=F) + labs(fill = "# Damaged houses", x = "", y = "") + geom_text(data = cap.df, aes(x = x, y = y, label = text), size = 2.5, hjust = 1) +
  theme(legend.position =  c(0.025, 0.15), legend.justification = c(0, 0), legend.title = element_text(size = 7.5), legend.text = element_text(size = 7.5))

ggsave("graphics/goni-damagedhouses.png", dpi = 300, height = 6, width = 10, units = "in")

# Plot pcprop
ggplot() + theme(
  panel.background = element_rect(fill = "lightblue",
                                  colour = "lightblue",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.2, linetype = 'dashed',
                                  colour = "white")
) +
  geom_sf(data = ph_mun_new, aes(fill=pcprop), size = 0.1, color = "black") + 
  scale_fill_gradient(low="white", high="red")  +
  geom_sf(data = goni_trace, size = 0.2, color = "red", fill = "cyan1") + 
  coord_sf(xlim=c(118,127), ylim=c(11,16), expand=F) + labs(fill = "Complete damage rate", x = "", y = "") + geom_text(data = cap.df, aes(x = x, y = y, label = text), size = 2.5, hjust = 1) +
  theme(legend.position =  c(0.025, 0.15), legend.justification = c(0, 0), legend.title = element_text(size = 7.5), legend.text = element_text(size = 7.5))

# Plot pcprop standard errors
ggplot() + theme(
  panel.background = element_rect(fill = "lightblue",
                                  colour = "lightblue",
                                  size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.2, linetype = 'dashed',
                                  colour = "white")
) +
  geom_sf(data = ph_mun_new, aes(fill=pcprop.se), size = 0.1, color = "black") + 
  scale_fill_gradient(low="white", high="red") + 
  geom_sf(data = goni_trace, size = 0.2, color = "red", fill = "cyan1") + 
  coord_sf(xlim=c(118,127), ylim=c(11,16), expand=F) + labs(fill = "Complete damage rate \n standard error", x = "", y = "") + geom_text(data = cap.df, aes(x = x, y = y, label = text), size = 2.5, hjust = 1) +
  theme(legend.position =  c(0.025, 0.15), legend.justification = c(0, 0), legend.title = element_text(size = 7.5), legend.text = element_text(size = 7.5))

#鍦ㄦ崯澶辩巼鏂归潰锛屽彈鎵撳嚮鏈�澶氱殑10涓煄甯傛槸锛?
top10.dr <- sort(ph_mun_new$tdprop, decreasing = TRUE)[1:10]

top10.dr.df <- data.frame(ph_mun_new[ph_mun_new$tdprop %in% top10.dr, c("Mun_Name", "tdprop")])

top10.dr.df <- top10.dr.df[order(top10.dr.df$tdprop, decreasing = TRUE), ]

top10.dr.df[, c("Mun_Name", "tdprop")]

#鎸夊彈鎹熸埧灞嬫暟閲忚绠楋紝鍙楃伨鏈�澶氱殑10涓煄甯傛槸锛?
top10.db <- sort(ph_mun_new$tdnum, decreasing = TRUE)[1:10]

top10.db.df <- data.frame(ph_mun_new[ph_mun_new$tdnum %in% top10.db, c("Mun_Name", "tdnum")])

top10.db.df <- top10.db.df[order(top10.db.df$tdnum, decreasing = TRUE), ]

top10.db.df[, c("Mun_Name", "tdnum")]

#鏈�鍚庯紝鎴戜滑鎻愪緵浜嗗湴鍥剧殑KMZ鏂囦欢锛屼互渚块�氳繃GoogleEarth瀹炵幇鍙鍖栵細
spdf <- as(ph_mun_new, "Spatial")
spdf@data <- spdf@data[, c("Mun_Name", "tdprop", "tdprop.se", "tdnum", "tdnum.se", "pcprop", "pcprop.se", "pcnum", "pcnum.se", "Wind_Speed", "dist_first_impact", "Slope", "Elevation", "Population.density", "X..strong.roof.type", "Experience")]
colnames(spdf@data) <- c("Municipality", "Damage_rate", "Damage_rate_standard_error", "No_of_damaged_houses", "No_of_damaged_houses_standard_error", "Complete_damage_rate", "Complete_damage_rate_standard_error", "No_of_completely_damaged_houses", "No_of_completely_damaged_houses_standard_error", "Wind_speed", "Distance_to_first_impact", "Slope", "Elevation", "Population_density", "Proportion_of_strong_roof_type", "Experience")
spdf@data$Municipality <- as.character(spdf@data$Municipality)
# Replace pattern in string:
spdf@data$Municipality[spdf@data$Municipality == "GABALDON (BITULOK & SABANI)"] <- "GABALDON (BITULOK AND SABANI)"

shape = "http://maps.google.com/mapfiles/kml/pal2/icon18.png"

# Round variables

# To 2 decimal places
spdf@data[, which(!(colnames(spdf@data)%in% c("Municipality", "Damage_rate", "Complete_damage_rate", "No_of_damaged_houses", "No_of_completely_damaged_houses")))] <-  round(spdf@data[, which(!(colnames(spdf@data)%in% c("Municipality", "Damage_rate", "Complete_damage_rate", "No_of_damaged_houses", "No_of_completely_damaged_houses")))], digits = 2)

# To the nearest 5%
spdf@data[, "Damage_rate"] <- round(spdf@data[, "Damage_rate"]/0.05)*0.05
spdf@data[, "Complete_damage_rate"] <- round(spdf@data[, "Complete_damage_rate"]/0.05)*0.05

# To the nearest 100
spdf@data[, "No_of_damaged_houses"] <- round(spdf@data[, "No_of_damaged_houses"]/100)*100
spdf@data[, "No_of_completely_damaged_houses"] <- round(spdf@data[, "No_of_completely_damaged_houses"]/100)*100

head(spdf@data)

#install.packages("plotKML")
library(plotKML)
plotKML(spdf, 'graphics/goni_damage', colour = "Damage_rate", plot.labpt = TRUE, shape = shape, alpha = 0.75, size = 0.2, balloon = TRUE)


