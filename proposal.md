# Stats506 Final Project

## Topic
Performance of BH-Adjust in an atypical situation

## Introduction
In the multiple hypothesis testing problem, usually we care about the False Discovery Rate(FDR) instead of Type 1 Error. Benjamini-Hochberg procedure is commonly used to control the FDR in multiple hypothesis testing.  
**Benjamini-Hochberg procedure:**  
Begin by ordering the p-values in ascending order.
$$
p_{(1)} \le p_{(2)} \le \cdots \le p_{(n)},
$$
Fix a significant level $q \in [0,1]$. Let $i_0$ be the largest $i$ for which
$$
p_{(i)} \le \frac{i}{n} q.
$$
Reject all $H_{(i)}$ with $i \le i_0$. This procedure is proved in [Benjamini & Yekutieli (2001)](https://projecteuclid.org/euclid.aos/1013699998). If the joint distribution of the statistics (or joint dist. of the p-values) is PRDS(positive regression dependency on each one from a subset) on the set of true nulls $H_0$, then the Benjamini-Hochberg procedure BH(q) controls the FDR at level q. However, in practise, we do not always know the precise correlation structure. In [Benjamini & Yekutieli (2001)](https://projecteuclid.org/euclid.aos/1013699998), it also introduced a BY-Adjust procedure, which use rejection threshold $q=q/\sum_{i=1}^n \frac{1}{i}$ to control the FDR at level $q$. However, this procedure is quite conservative. In this project, I want to demostrate the performance of original BH procedure on various correlation structure.

## Experiment Setting
Consider $n$ i.i.d. samples $x_t$, $t = 1,...,n$, $x_t = [x_{t1},...,x_{tm}] \in R^m$. Suppose $x_{ti} \sim N(\mu_i,1)$. We want to test the null hypothesis $H_{0i} : \mu_i \le 0$ for $i = 1,...,m$. Consider the following setting:  
(a) **i.i.d.**: $x_{ti}$ for $i = 1, ... , m$ are i.i.d.  
(b) **Positive Correlated.**: $Cov(x_{ti}, x_{tj}) = 0.8$ for $i\ne j$.  
(c) **Negative Correlated.**: $Cov(x_{ti}, x_{tj}) = 0.8$ for $1 \le i <j \le m/2$, $x_{t, i+m/2} = -x_{t,i}$ for $1 \le i \le m/2$.  
 
Here we apply the BH-Adjust and BY-Adjust at level $q=0.1$ by conducting $B=1000$ simulations with $n=200$ and $m=100$ in each setting, and compare the FDR and power of these adjustment procedure. We will try the following cases:
1. Global Null: $\mu_i=0$
2. $\mu_1=\mu_2=1, \mu_3=\mu_0.6, \mu_5=\mu_6=0.2$, and $\mu_i=0$ for $i>6$.