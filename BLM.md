# Black Litterman Model

## 1. Intro: Markowitz Mean-variance Model

find w to maximize $w^T\mu - \lambda w' \Sigma w/2 $

Weaknesses:
1. hard to give $\mu$
2. historical return is changing based on the period
3. no constraints on the trade

## 2. Expected Returns
### 1. Reverse Optimization
Based on the market information, we can gain *implied market equilibrium return*, using reverse optimization;

$\Pi = \lambda \Sigma w_{mkt}$

This return is identical to the CAPM implied return, and also is very different from the historical return.

With implied market equilibrium return, we can combine our views to change the distribution of the return, and we call it *new combined return*;

### 2. Black-Litterman formula
The new combined return has a new distribution, and this expectation is 

$E[R]=[(\tau\Sigma)^{-1}+P'\Omega^{-1}P]^{-1}[(\tau\Sigma)^{-1}\Pi+P'\Omega Q]$

- Views

    View 1: Asset A will have an absolute excess return of 5.25% (Confidence of View = 25%).
    
    View 2: B will outperform C by 0.25% (Confidence of View = 50%).
    
    View 3: D and E will outperform F and G value by 2% (Confidence of View = 65%).

- Related matrix：

    $$P:=\quad\left[ \begin{matrix}0&0&0&0&0&0&1&0\\
    -1&1&0&0&0&0&0&0\\
    0&0&.9&-.9&.1&-.1&0&0\end{matrix}\right]
    $$
    

$$Q:=\quad\left[ \begin{matrix}0.0525\\0.0025\\0.02\end{matrix}\right]$$

$$\Omega:=diag(0.25*0.25,0.5*0.5,0.65*0.65)$$

- $\tau$:scalar,but guidance in the literature for setting the scalar’s value is scarce.
    - 1.$\tau$ is close to zero,from 0.01 to 0.05(black, litterman,1992)
    - 2.$\tau=1$,Satchell and Scowcroft,2000
    - 3.$\tau=0.25$
    - 4.$\tau\Sigma$ is the standard of the estimate of the Implied Equilibrium Return Vector ($\Pi$) ,Blamont and Firoozye (2003),He and Litterman (1999) calibrate the confidence of a view so that the ratio of $\omega/\tau$ is equal to the variance of the view portfolio ( $p_k'\Sigma p_k$ ).

    if the last estimate is used, then the $\Omega$ matrix can be specified autonomously.

### 3. Optimization
using the *new combined return*, substitute it to the Markowitz Optimazation model,

i.e. maximize $w^T\mu_p - \lambda w^T \Sigma w/2 $

the w is the optimal weight of the portfolio.

In reality, we have the constraints:

$$\Sigma w_i = 1;\\w_i>=0$$




```python

```
