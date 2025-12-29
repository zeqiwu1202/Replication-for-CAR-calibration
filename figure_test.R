# 加载必要的包
library(ggplot2)
library(dplyr)

# 设置随机种子，保证结果可复现
set.seed(123)

# 样本大小
sample_sizes <- c(500, 1000, 1500, 2000)

# 方法
methods <- c("DELL", "PDELL", "AIPW", "IPW")

# 生成数据框
data <- expand.grid(SampleSize = sample_sizes, Method = methods, Replicate = 1:30)

# 为每个方法生成不同的值
data <- data %>%
  mutate(Value = case_when(
    Method == "DELL"  ~ rnorm(n(), mean = 0.45, sd = 0.05),
    Method == "PDELL" ~ rnorm(n(), mean = 0.5, sd = 0.1),
    Method == "AIPW"  ~ rnorm(n(), mean = 0.8, sd = 0.3),
    Method == "IPW"   ~ rnorm(n(), mean = 1.5 - 0.5 * (SampleSize / 1000), sd = 0.6)
  ))

# 真实值（基准线）
real_value <- 0.5

# 计算每个方法在每个样本大小下的均值
method_means <- data %>%
  group_by(SampleSize, Method) %>%
  summarise(MeanValue = mean(Value), .groups = 'drop')

# 绘图
ggplot(data, aes(x = factor(SampleSize), y = Value, fill = Method)) +
  geom_boxplot(outlier.size = 1, position = position_dodge(width = 0.8), width = 0.7) +  # 箱线图
  geom_line(data = method_means, aes(x = factor(SampleSize), y = MeanValue, color = Method, group = Method), 
            size = 1) +  # 方法的均值线
  geom_point(data = method_means, aes(x = factor(SampleSize), y = MeanValue, color = Method), 
             size = 2, shape = 3) +  # 添加均值点
  geom_hline(yintercept = real_value, linetype = "dashed", color = "purple", size = 1) +  # 真实值的虚线
  scale_fill_manual(values = c("DELL" = "red", "PDELL" = "blue", "AIPW" = "orange", "IPW" = "green")) + 
  scale_color_manual(values = c("DELL" = "red", "PDELL" = "blue", "AIPW" = "orange", "IPW" = "green")) +
  labs(
    title = "Boxplots of Estimated Average Treatment Effects",
    x = "Sample Size",
    y = "Average Treatment Effects"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.85, 0.85),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  )
