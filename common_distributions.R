# Load packages
library(ggplot2)
library(tidyverse)
library(tidyr)
library(ggsci)
library(cowplot)
library(ggpubr)

# distributions

set.seed(12345)

# Generate representative samples
n = 5000

df = tibble(
  Normal = rnorm(n, mean = 0, sd = 1),
  t = rt(n, df = 5),
  ChiSq = rchisq(n, df = 4),
  Binomial = rbinom(n, size = 20, prob = 0.5),
  NegBinomial = rnbinom(n, size = 5, mu = 10),
  Poisson = rpois(n, lambda = 3)
)

# Convert to long form for ggplot2
df_long = df %>%
  pivot_longer(cols = everything(), names_to = "distribution", values_to = "value")

# Plot
distros = ggplot(df_long, aes(x = value, fill = distribution)) +
  geom_histogram(bins = 50) +
  facet_wrap(~ distribution, scales = "free") +
  theme_bw(base_size = 14, base_family = "Arial") +
  scale_fill_nejm() +
  theme(legend.position = "none") +
  labs(x = "", y = "")

distros

png("~/OneDrive - Northwestern University/Teaching/Choosing the Right Test/representative_distros.png",
    width = 16,
    height = 9,
    units = "in",
    res = 300)
distros
dev.off()
   
# Central limit theorem
set.seed(12345)
ns = c(5, 10, 50, 100, 500, 1000, 5000, 100000)
distros = lapply(ns, function(n){
  rnorm(n, mean = 0, sd = 1) })
plots_normal = lapply(distros, function(d){
  ggplot(NULL, aes(x = d)) +
    geom_histogram(fill = "dodgerblue4") +
    labs(x = "", y = "Count", title = paste("n = ", length(d))) +
    theme_bw(base_family = "Arial",
             base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5)) })

clt_normal = plot_grid(plotlist = plots_normal, ncol = 4)
png("~/OneDrive - Northwestern University/Teaching/Choosing the Right Test/central_limit_normal.png",
    width = 16,
    height = 9,
    units = "in",
    res = 300)
clt_normal
dev.off()

set.seed(12345)
ns = c(5, 10, 50, 100, 500, 1000, 5000, 100000)
distros = lapply(ns, function(n){
  rnbinom(n, size = 5, mu = 10) })
plots_negbinom = lapply(distros, function(d){
  ggplot(NULL, aes(x = d)) +
    geom_histogram(fill = "dodgerblue4") +
    labs(x = "", y = "Count", title = paste("n = ", length(d))) +
    theme_bw(base_family = "Arial",
             base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5)) })

clt_negbinom = plot_grid(plotlist = plots_negbinom, ncol = 4)
png("~/OneDrive - Northwestern University/Teaching/Choosing the Right Test/central_limit_negbinom.png",
    width = 16,
    height = 9,
    units = "in",
    res = 300)
clt_negbinom
dev.off()   

# IL-6 ranges (nomimal vs biological significance)
il6 = data.frame(Condition = factor(c(rep("Healthy", 20), 
                               rep("Alzheimer's Disease", 20),
                               rep("Sepsis", 20)),
                               levels = c("Healthy", "Alzheimer's Disease", "Sepsis")),
                 pg_ml = c(rnorm(n = 20, mean = 2, sd = 1),
                                  rnorm(n = 20, mean = 10, sd = 1),
                                  rnorm(n = 20, mean = 800, sd = 100)))

il6_ad = il6 %>% 
  dplyr::filter(Condition != "Sepsis") %>% 
  ggplot(aes(x = Condition, y = pg_ml, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  stat_compare_means(comparisons = list(c("Healthy", "Alzheimer's Disease")),
                     method = "wilcox.test",
                     p.adjust.method = "fdr") +
  theme_bw(base_family = "Arial", base_size = 20) +
  theme(legend.position = "none") +
  labs(x = "Diagnosis", y = "Plasma [IL-6] (pg/mL)") +
  scale_fill_nejm()

png("~/OneDrive - Northwestern University/Teaching/Choosing the Right Test/il6_ad.png",
    width = 8,
    height = 9,
    units = "in",
    res = 300)
il6_ad
dev.off()   

il6_ad

il6_all = il6 %>% 
  ggplot(aes(x = Condition, y = pg_ml, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  stat_compare_means(comparisons = list(c("Healthy", "Alzheimer's Disease"),
                                        c("Healthy", "Sepsis"),
                                        c("Alzheimer's Disease", "Sepsis")),
                     method = "wilcox.test",
                     p.adjust.method = "fdr") +
  theme_bw(base_family = "Arial", base_size = 20) +
  theme(legend.position = "none") +
  labs(x = "Diagnosis", y = "Plasma [IL-6] (pg/mL)") +
  scale_fill_nejm()

png("~/OneDrive - Northwestern University/Teaching/Choosing the Right Test/il6_all.png",
    width = 8,
    height = 9,
    units = "in",
    res = 300)
il6_all
dev.off()   


il6_all

cxcl10_two_main = data.frame(Condition = factor(c(rep("Healthy, Untreated", 20), 
                                      rep("Healthy, Dex", 20),
                                      rep("Pneumonia, Untreated", 20),
                                      rep("Pneumonia, Dex", 20)),
                                    levels = c("Healthy, Untreated", "Healthy, Dex", 
                                               "Pneumonia, Untreated","Pneumonia, Dex")),
                 pg_ml = c(rnorm(n = 20, mean = 2, sd = 1),
                           rnorm(n = 20, mean = 0, sd = 1),
                           rnorm(n = 20, mean = 50, sd = 10),
                           rnorm(n = 20, mean = 25, sd = 10)))

two_main_plot = cxcl10_two_main %>% 
  ggplot(aes(x = Condition, y = pg_ml, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  stat_compare_means(comparisons = list(c("Healthy, Untreated", "Healthy, Dex"),
                                        c("Pneumonia, Untreated","Pneumonia, Dex"),
                                        c("Healthy, Untreated", "Pneumonia, Untreated"),
                                        c("Healthy, Dex", "Pneumonia, Dex")),
                     method = "t.test",
                     p.adjust.method = "fdr") +
  theme_bw(base_family = "Arial", base_size = 20) +
  theme(legend.position = "none") +
  labs(x = "", y = "Plasma [CXCL10] (pg/mL)") +
  scale_fill_nejm()
png("~/OneDrive - Northwestern University/Teaching/Choosing the Right Test/two_mains.png",
    width = 10,
    height = 9,
    units = "in",
    res = 300)
two_main_plot
dev.off()   

set.seed(53)
cxcl10_no_main = data.frame(Condition = factor(c(rep("Healthy, Untreated", 5), 
                                                  rep("Pneumonia, Untreated", 5),
                                                  rep("Pneumonia, Dex", 5)),
                                                levels = c("Healthy, Untreated",
                                                           "Pneumonia, Untreated",
                                                           "Pneumonia, Dex")),
                             pg_ml = c(rnorm(n = 5, mean = 2, sd = 1),
                                       rnorm(n = 5, mean = 50, sd = 10),
                                       rnorm(n = 5, mean = 50, sd = 50)))

no_main_plot = cxcl10_no_main %>% 
  ggplot(aes(x = Condition, y = pg_ml, fill = Condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  stat_compare_means(comparisons = list(c("Healthy, Untreated", "Pneumonia, Untreated"),
                                        c("Healthy, Untreated", "Pneumonia, Dex")),
                     method = "wilcox.test",
                     p.adjust.method = "fdr") +
  theme_bw(base_family = "Arial", base_size = 20) +
  theme(legend.position = "none") +
  labs(x = "Diagnosis", y = "Plasma [CXCL10] (pg/mL)") +
  scale_fill_nejm()
png("~/OneDrive - Northwestern University/Teaching/Choosing the Right Test/no_main.png",
    width = 10,
    height = 9,
    units = "in",
    res = 300)
no_main_plot
dev.off()   