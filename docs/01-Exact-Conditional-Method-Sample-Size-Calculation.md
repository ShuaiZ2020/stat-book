# Exact Conditional Method疫苗有效性(Vaccine Efficacy)的样本量计算

我们首先导入后面会用到的包，后续会用到dplyr的%>%管道函数和kableExtra输出表格。

``` r
library(dplyr)
library(kableExtra)
```


## 疫苗有效性(Vaccine Efficacy)

- 当我们检验一款新的疫苗是否有效时，通常会使用疫苗有效性$\pi$(Vaccine Efficacy)这个指标。

- 让$P_1$, $P_2$分别为安慰剂组和接种疫苗组的发病率；$N_1$, $N_2$分别为安慰剂组和接种疫苗组的群体总数；$X, Y$分别为安慰剂和疫苗组的病例人数。则$X \overset{i.i.d}{\sim}\text{Binomial}(N_1, P_1)$, $Y \overset{i.i.d}{\sim}\text{Binomial}(N_2, P_2)$。 

- 疫苗有效性可以可以用可以被写成：
\begin{equation} 
  \pi = 1-(P_2/P_1)
  (\#eq:pi)
\end{equation} 
- 原假设和备择假设可以写成

\begin{equation} 
  H_0: \pi \leq \pi_0 \text{ versus } \ H_1:\pi>\pi_0   
  (\#eq:hypo1)
\end{equation} 



## 基于泊松假设的大样本

- 如果一个疾病发病率较低，我们就需要更多受试者参加实验，在这种情况下$X, Y$可以被近似成为独立的泊松分布. 

- $X \overset{i.i.d}{\sim}\text{Poisson}(\lambda_1)$, $Y \overset{i.i.d}{\sim}\text{Poisson}(\lambda_2)$ with $\lambda_1 = N_1\cdot P_1, \lambda_2 = N_2\cdot P_2$. 
- $Y$的条件概率分布可以被写成 
$Y|T \sim \text{Binomial}(T, \theta)=\binom{T}{k} \theta^k(1-\theta)^{T-k}$ with $T = X + Y, \theta = \frac{\lambda_1}{\lambda_1=\lambda_2}$

- 此时原假设和备择假设可以写成

\begin{equation} 
  H_0: \theta \geq \theta_0 \text{ versus } \ H_1:\theta<\theta_0 \text{ where } \theta_0 = \frac{1-\pi_0}{2-\pi_0}
  (\#eq:hypo2)
\end{equation} 


- P value 可以被计算为
\begin{equation} 
  p = \Pr[Y \leq Y_{\text{obs}} \mid Y \sim \text{Binomial}(T, \theta_0)] = \sum_{k=0}^{Y_{\text{obs}}} \binom{T}{k} \theta_0^k (1 - \theta_0)^{T - k}
  (\#eq:pvalue)
\end{equation} 
 
- Statistical power 可以被计算为
\begin{equation} 
  1 - \beta = \Pr[Y \leq Y_c \mid Y \sim \text{Binomial}(T, \theta_1)] = \sum_{k=0}^{Y_c} \binom{T}{k} \theta_1^k (1 - \theta_1)^{T - k}
  (\#eq:power)
\end{equation} 

- 临床实验所需样本量可以被计算为
\begin{equation} 
  N_2 = T/[(2-\pi_1)/P_1]
  (\#eq:n2)
\end{equation} 

- 计算样本量所需的变量为$\pi_0, \pi_1, \alpha$, 期望统计功效(1-$\beta$)和安慰剂组发病率($P_1$)。

## R代码示例


``` r
exact_conditonal_test <- function(T, alpha, power, theta0, theta1) {
  Y_c <- qbinom(alpha, size = T, prob = theta0)-1
  p_value <- pbinom(Y_c, size = T, prob = theta0)
  power <- pbinom(Y_c, size = T, prob = theta1)
  return(c(T, Y_c, power, p_value))
}

ect_sample_size <- function(T, pi0, pi1, incidence, alpha, power){
  theta0 <- (1-pi0) / (2-pi0)
  theta1 <- (1-pi1) / (2-pi1)
  table_t <- as.data.frame(do.call(rbind, 
                               lapply(T, exact_conditonal_test, alpha = alpha, theta0= theta0, theta1 = theta1)), .name_repair = "unique")
  colnames(table_t) <- c('T', 'Y_c', 'power', 'p-value')
  
  for (n in 1:nrow(table_t)){
    if (table_t$power[n] >= power && all(table_t$power[n:nrow(table_t)] >= power)) {
      min_n <- table_t$T[n]
      break
    }
  }
  result <- list(
    T_table = table_t,
    T = min_n,
    text = paste0('The min value of T achieve power of ', power, ' is ', min_n)
  )
  class(result) <- 'result'
  result
}

print.myresult <- function(x, ...) {
  cat(x$text, "\n")
}


T2N <- function(T_value, incidence, pi1, dropout_rate){
  N2 <- T_value/((2-pi1)*incidence)/(1-dropout_rate)
  cat('The sample of vaccine group considering the drop out rate:', N2)
}

```

详细推导过程可以看 @chan_exact_1998 和 @loiacono_sample_nodate

### 例1
@chan_exact_1998 这篇论文介绍了exact conditional method的公式推导。他还给了一个不同样本量对应的power和significance level。 我们可以运行上面的函数来验证结果的准确性。设置病例范围从33到40，$\pi_0=0.2,\pi_1=0.8, P_1=0.006, \alpha=0.025$, 期望统计功效为95%。输入以上参数到ect_sample_size这个函数。

<img src="image/Screenshot 2025-07-21 155606.png" width="334" />


``` r
T <- 33:40 # the number of T
pi0 <- 0.2 # null hypothesis efficacy
pi1 <- 0.8 # True efficacy under alternative hypothesis
alpha <- 0.025 # type I error
incidence <- 0.006 # placebo incidence rate
power <- 0.95
res1 <- ect_sample_size(T, pi0, pi1, incidence, alpha, power)
kbl(res1$T_table)%>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left")
```

<table class="table table-striped" style="width: auto !important; ">
 <thead>
  <tr>
   <th style="text-align:right;"> T </th>
   <th style="text-align:right;"> Y_c </th>
   <th style="text-align:right;"> power </th>
   <th style="text-align:right;"> p-value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 0.9139690 </td>
   <td style="text-align:right;"> 0.0136117 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.9540856 </td>
   <td style="text-align:right;"> 0.0244451 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 35 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.9449925 </td>
   <td style="text-align:right;"> 0.0178969 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 36 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.9347919 </td>
   <td style="text-align:right;"> 0.0129998 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.9653937 </td>
   <td style="text-align:right;"> 0.0227940 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 38 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.9584044 </td>
   <td style="text-align:right;"> 0.0168288 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.9504998 </td>
   <td style="text-align:right;"> 0.0123313 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.9738542 </td>
   <td style="text-align:right;"> 0.0211901 </td>
  </tr>
</tbody>
</table>


代码中的incidence 就是$P_1$， 其他参数与上面提及的保持一致。可以看到输出的表格中样本量，critical value， statistical power和p value与论文中的表格完全一致。我们打印出能够达到95% statistical power的病例数

``` r
res1$T
#> [1] 37
```
打印出T=37与论文一致，再用公式\@ref(eq:n2)计算出疫苗组样本量，为了方便计算我也写成了函数t2n。此示例计算样本量时没有考虑脱落所以dropout_rate 参数可以设置为0。得到结果样本量至少为5138.889。这个结果乘2就是疫苗组和安慰剂组所需的总样本量，总样本量至少为10277.78，进一10278。此结果与论文一致，证明我们的算法无误。

``` r
T2N(37, incidence, pi1, dropout_rate = 0)
#> The sample of vaccine group considering the drop out rate: 5138.889
```

### HRV-三期
<img src="image/Screenshot 2025-07-22 173413.png" width="284" />

@wu_efficacy_2022 这篇论文原为我们想要复现的样本量。但是过程中发现其中并未明确提及所需参数 $\pi_1, \alpha$和期望统计功效。无法计算出样本量所以后续我们使用这篇论文作为参考。下一个示例为HRV-三期这篇论文的参考文献，也是rotavirus vaccine的疫苗有效性临床实验。

### RV5 
@mo_efficacy_2017 这篇论文会作为我们后续疫苗有效性检验样本量计算的参考文章。提供参数：15%脱落率， $\pi_0=0, \pi_1=0.6, P_1 = 0.02, \alpha = 0.025$，期望统计功效为80%。带入以上参数可以得到病例总数至少为47，我们取双数48。再用T2N函数计算疫苗组样本量总数至少2016.807，我们取整为2020，与论文一致。


<img src="image/Screenshot 2025-07-21 162616.png" width="328" />


``` r
T <- 40:50 # the number of T
pi0 <- 0 # null hypothesis efficacy
pi1 <- 0.6 # True efficacy under alternative hypothesis
alpha <- 0.025 # type I error
incidence <- 0.02 # placebo incidence rate
power <- 0.8
res_rv5 <- ect_sample_size(T, pi0, pi1, incidence, alpha, power)
kbl(res_rv5$T_table)%>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left")
```

<table class="table table-striped" style="width: auto !important; ">
 <thead>
  <tr>
   <th style="text-align:right;"> T </th>
   <th style="text-align:right;"> Y_c </th>
   <th style="text-align:right;"> power </th>
   <th style="text-align:right;"> p-value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.7692914 </td>
   <td style="text-align:right;"> 0.0192387 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 41 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 0.7363326 </td>
   <td style="text-align:right;"> 0.0137666 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 42 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.8052771 </td>
   <td style="text-align:right;"> 0.0217793 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 43 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 0.7757295 </td>
   <td style="text-align:right;"> 0.0157697 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.8362319 </td>
   <td style="text-align:right;"> 0.0243834 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 45 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.8100042 </td>
   <td style="text-align:right;"> 0.0178489 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 46 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 0.7819032 </td>
   <td style="text-align:right;"> 0.0129480 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 47 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.8396107 </td>
   <td style="text-align:right;"> 0.0199930 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 0.8146130 </td>
   <td style="text-align:right;"> 0.0146525 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 49 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.8650285 </td>
   <td style="text-align:right;"> 0.0221921 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 50 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 0.8429717 </td>
   <td style="text-align:right;"> 0.0164196 </td>
  </tr>
</tbody>
</table>



``` r
res_rv5$T
#> [1] 47
```


``` r
T2N(48, incidence, pi1, dropout_rate = 0.15)
#> The sample of vaccine group considering the drop out rate: 2016.807
```

## Reference{-}




