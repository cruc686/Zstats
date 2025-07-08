#' My Function
#'
#' This function calculate mediation effection.
#'
#' @param beta1, se1 beta2,se2,beta_total,n.
#' @return data.frame.
#' @export
#' 
#' 
#beta1、se1：XM的回归系数和标准误，beta2、se2：MY的回归系数和标准误，beta_total:XY总效应回归系数
.med_cal <- function(beta1,se1,beta2,se2,beta_total,n=1000){
  mediate_effect <- beta1 * beta2
  mediate_se <- sqrt((beta1^2 * se2^2) + (beta2^2 * se1^2))
  m_low <- mediate_effect - 1.96* mediate_se
  m_upper <- mediate_effect + 1.96* mediate_se
  Z <- mediate_effect/mediate_se
  p <- 2*pnorm(q=abs(Z), lower.tail=FALSE)
  # 计算总效应
  total_effect <- beta_total
  # 计算直接效应（总效应减去中介效应）
  direct_effect <- total_effect - mediate_effect
  # 计算中介效应百分比
  mediate_percentage <- (mediate_effect / total_effect) * 100
  # Bootstrap设置
  mediate_percentages <- numeric(n)
  # 执行Bootstrap
  set.seed(123)  # 确保结果可重复
  for (i in 1:n) {
    beta1_boot <- rnorm(1, beta1, se1)
    beta2_boot <- rnorm(1, beta2, se2)
    mediate_effect_boot <- beta1_boot * beta2_boot
    direct_effect_boot <- total_effect - mediate_effect_boot
    mediate_percentages[i] <- (mediate_effect_boot / total_effect) * 100
  }
  # 计算中介效应百分比的可信区间
  pro <- (beta1*beta2 / total_effect) * 100
  se <- sd(mediate_percentages)
  pro_low <- pro - 1.96*se 
  pro_upper <- pro + 1.96*se 
  res <- data.frame(mediate_effect,m_low,m_upper,p,pro,pro_low,pro_upper)
  return(res)
}
