# Statistical Process Monitoring

Statistical Process Monitoring (SPM) involves using various tools to assess process stability. Here's an overview of the terminology, types of control charts, and methodologies used in SPM. 

## Basic terminology

Control charts are fundamental tools for assessing process stability under the SPM framework. 
In the online (Phase II) setting, the process is monitored sequentially as new data is collected. The goal is to detect deviations from the in-control (IC) state as soon as possible.

- **Monitoring Statistic**: Sequential calculation of a statistic $C_{t}$ for $t = 1, 2, \ldots$. An alarm is raised when $C_{t}$ falls outside the control limits $(\text{LCL}_t, \text{UCL}_t)$.
- **Run Length (RL)**: represents the number of time points required for the monitoring procedure to signal an alarm, $\text{RL} = \inf\left\{ t > 0: C\\_{t} > \text{UCL}\\_{t} \text{ or } C\\_{t} < \text{LCL}\\_{t} \right\}$.
    
- **Control Limits**: Selected to constrain some IC properties of the chart's RL to a nominal value. For Phase II control charts, a common design is $\text{ARL}_\text{IC} := \mathbb{E}_{0}[\text{RL}] = A_0$
    where $A_0 > 1$ and $\mathbb{E}_{0}[\cdot ]$ represents the expectation assuming the process always remains IC.
    Other designs use the median of the IC run length $\text{MRL}_\text{IC}$ or the run length's quantiles.

## Taxonomy of control charts

Control charts can be classified into three main categories:

1. Shewhart-type: memoryless, reactive to large changes. Uses only the information about $\bm{X}_t$ at each time $t > 0$ to compute the monitoring statistic, $C_{t} = f(\boldsymbol{X}_t)$.

2. CUSUM-type: chart with memory, dampens the historical information with an update mechanism of the form $C_{t} = \max\left\{ 0, C_{t-1} + f(\bm{X}_t) \right\}.$
3. EWMA-type: chart with memory, the historical information is weighted using exponentially-decaying weights such as $C_{t} = (1 - \lambda)C_{t-1} + \lambda X_{t}.$

## Nonparametric control charts

Traditional control charts rely on i.i.d. continuous quality variables following a parametric distribution. When these assumptions are violated, control charts designed under these assumptions can display poor performance. 

Various nonparametric methods have been developed to relax the parametric assumptions:
- **Rank-Based Charts**: these use a rank transformation of the data.
- **Data Categorization**: numerical data is categorized and then a log-linear model is subsequently monitored over time.

These methods help when parametric assumptions are infeasible but might lose effectiveness in information compared to parametric charts.

## Selecting hyperparameters

Choosing the appropriate values of tuning parameters $\bm{\zeta} \in \mathcal{Z} \subseteq \mathbb{R}^{d}$ (e.g., smoothing constant $\lambda$ in EWMA, allowance constant $k$ in CUSUM) is crucial.

The choice of tuning parameters aims to detect specific magnitudes of parameter change efficiently, typically formulated as:

$\begin{aligned}
  &\frac{\partial \mathbb{E}_1[\text{RL}]}{\partial \bm{\zeta}}\Big|_{\bm{\zeta}=\bm{\zeta}^*} = \bm{0},\\
    & \text{s.t. } \mathbb{E}_0[\text{RL}] = A_0,
\end{aligned}$

When analitycal solutions are unavailable, numerical optimization methods like stochastic approximations or Monte-Carlo simulations are used.

## Multi-chart monitoring schemes

In complex monitoring scenarios, multiple control charts may be run simultaneously. These schemes are useful for monitoring multiple parameters jointly, such as the mean and variance of a distribution.

The control limits $\bm{LCL}_t = (\text{LCL}_{t1}, \ldots, \text{LCL}_{tJ})$ and $$\bm{UCL}_t = (\text{UCL}_{t1}, \ldots, \text{UCL}_{tJ})$$ are specific to each control chart, and are usually determined so that:

$\begin{cases}
\mathbb{E}_{0}\left[ \min( \text{RL}_1, \text{RL}_2, \ldots, \text{RL}_J) \right] = A_0,\\
\mathbb{E}_{0}[\text{RL}_1] = \mathbb{E}_{0}[\text{RL}_2] = \ldots = \mathbb{E}_{0}[\text{RL}_J].
\end{cases}$

This design assumes equal importance of each control chart, but weighting schemes might be used to emphasize the relative importance of each control chart. Algorithms like stochastic approximations are commonly used to solve these designs.
Furthermore, the expectation can be replaced by the median or quantiles of the RL in the above scheme.