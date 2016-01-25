args <- commandArgs(trailingOnly = TRUE)

table_file_name <- args[1]
plot_file_name <- args[2]

table <- read.table(table_file_name,header=T)
#cor.test(table$log_length,table$proportion_consistent,method="kendall",alternative = "less")
pdf(plot_file_name, height=7, width=7)
plot(table$log_length,table$proportion_consistent,xlab="Log branch duration",ylab="Proportion of indels with RI=1.0")
fit <- lm(table$proportion_consistent ~ table$log_length)
summary(fit)
r2 <- summary(fit)$r.squared
coef_string <- fit$coefficients["(Intercept)"]
coef <- as.numeric(substr(coef_string, 0, 10))
text(0.7,coef-0.1,"r2 = ")
text(1.05,coef-0.1,round(r2,digits=3))
dev.off()