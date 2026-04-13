setwd("/Users/MycologyLab/yanqi/DADA2/OIL_SPIL_SAMPLE/new_0410/")
# =========================
# 0. 环境
# =========================
library(vegan)
library(ggplot2)
library(minpack.lm)

# =========================
# 1. 读取数据
# =========================
dat <- read.csv("asv_t2o.csv", row.names = 1, check.names = FALSE)

# 去掉Sequence列
comm <- dat[, !(colnames(dat) %in% "Sequence")]

# 转置 → 行=样本
comm <- t(comm)

# 去掉全0 ASV
comm <- comm[, colSums(comm) > 0]

# =========================
# 2. 计算输入变量
# =========================

# 相对丰度
relab <- decostand(comm, method = "total")

# 平均相对丰度
p <- colMeans(relab)

# 出现频率（presence frequency）
freq <- colSums(comm > 0) / nrow(comm)

df <- data.frame(p = p, freq = freq)

# 去掉0（避免log问题）
df <- df[df$p > 0, ]

df <- df[is.finite(df$p) & is.finite(df$freq), ]

# 去掉极端值（非常重要）
df <- df[df$p > 1e-6 & df$freq > 0 & df$freq < 1, ]

# log转换
df$logp <- log10(df$p)

df <- df[is.finite(df$p) & is.finite(df$freq), ]

# 去掉极端值（非常重要）
df <- df[df$p > 1e-6 & df$freq > 0 & df$freq < 1, ]

# =========================
# 3. 拟合 Neutral model
# =========================

# Sloan model（简化版）
fit <- nlsLM(
  freq ~ 1 - exp(-m * p),
  data = df,
  start = list(m = 0.1),
  control = nls.lm.control(maxiter = 1000)
)

m <- coef(fit)[["m"]]

# 预测
df$pred <- predict(fit)

# =========================
# 4. 计算R²
# =========================
R2 <- cor(df$freq, df$pred)^2

# =========================
# 5. 置信区间（经验方法）
# =========================
df$upper <- df$pred + 1.96 * sqrt(df$pred * (1 - df$pred) / nrow(comm))
df$lower <- df$pred - 1.96 * sqrt(df$pred * (1 - df$pred) / nrow(comm))

# 限制范围
df$upper[df$upper > 1] <- 1
df$lower[df$lower < 0] <- 0

# =========================
# 6. 分类点（above / neutral / below）
# =========================
df$class <- "neutral"
df$class[df$freq > df$upper] <- "above"
df$class[df$freq < df$lower] <- "below"

# =========================
# 7. 作图（论文级）
# =========================
ggplot(df, aes(x = logp, y = freq)) +
  geom_point(aes(color = class), size = 2) +
  
  geom_line(aes(y = pred), color = "blue", size = 1) +
  geom_line(aes(y = upper), color = "blue", linetype = "dashed") +
  geom_line(aes(y = lower), color = "blue", linetype = "dashed") +
  
  scale_color_manual(values = c(
    "above" = "#d73027",
    "neutral" = "black",
    "below" = "#4575b4"
  )) +
  
  theme_classic() +
  
  labs(
    x = "log(Mean relative abundance)",
    y = "Occurrence frequency"
  ) +
  
  annotate("text",
           x = min(df$logp),
           y = 0.9,
           hjust = 0,
           label = paste0(
             "R² = ", round(R2, 4),
             "\nm = ", signif(m, 3)
           ))
p <- ggplot(df, aes(x = logp, y = freq)) +
  geom_point(aes(color = class), size = 2) +
  
  geom_line(aes(y = pred), color = "blue", size = 1) +
  geom_line(aes(y = upper), color = "blue", linetype = "dashed") +
  geom_line(aes(y = lower), color = "blue", linetype = "dashed") +
  
  scale_color_manual(values = c(
    "above" = "#d73027",
    "neutral" = "black",
    "below" = "#4575b4"
  )) +
  
  theme_classic() +
  
  labs(
    x = "log(Mean relative abundance)",
    y = "Occurrence frequency"
  ) +
  
  annotate("text",
           x = min(df$logp),
           y = 0.9,
           hjust = 0,
           label = paste0(
             "R² = ", round(R2, 4),
             "\nm = ", signif(m, 3)
           ))

ggsave("neutral_model_plot_t2o.pdf", p, width = 4, height = 3)
