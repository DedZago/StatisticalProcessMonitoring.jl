# Statistical Process Monitoring

Statistical Process Monitoring (SPM) involves using various tools to assess process stability. Here's an overview of the terminology, types of control charts, and methodologies used in SPM. 

## Basic Terminology

Control charts are pivotal tools for assessing process stability under SPM. According to [woodall2000](@citet), control charts can be broadly categorized into Phase I and Phase II control charts.

### Phase II Control Charts
- **Objective**: Monitor deviations from the IC state as new data is collected.
- **Monitoring Statistic**: Sequential calculation of a statistic $C_{t}$ for $t = 1, 2, \ldots$. An alarm is raised when $C_{t}$ falls outside the control limits $(\text{LCL}_t, \text{UCL}_t)$.
- **Run Length (RL)**:
    $$
    \text{RL} = \inf\left\{ t > 0: C_{t} > \text{UCL}_{t} \text{ or } C_{t} < \text{LCL}_{t} \right\}
    $$
    RL represents the number of time points required for the monitoring procedure to signal an alarm.
- **Control Limits**: Selected to constrain some IC properties of the chart's RL to a nominal value. For Phase II control charts, a common design is:
    $$
    \text{ARL}_\text{IC} := \mathbb{E}_{0}[\text{RL}] = A_0,
    $$
    where $A_0 > 1$ and $\mathbb{E}_{0}[\cdot ]$ represents the expectation assuming the process always remains IC.
    Other designs use the median of the IC run length $\text{MRL}_\text{IC}$ or the run length's quantiles.

## Taxonomy of Control Charts

Control charts can be classified into three main categories as per [chakraborti2001](@citet):

1. Shewhart-type: memoryless, reactive to large changes. Uses only the information about $\bm{X}_t$ at each time $t > 0$.

2. CUSUM-type: chart with memory, dampens the historical information with an update mechanism of the form
    $$
    C_{t} = \max\left\{ 0, C_{t-1} + f(\bm{X}_t) \right\}.
    $$
3. EWMA-type: chart with memory, the historical information is weighted using exponentially-decaying weights such as
    $$
    C_{t} = (1 - \lambda)C_{t-1} + \lambda X_{t}.
    $$

## Nonparametric Control Charts

Traditional control charts rely on i.i.d. continuous quality variables following a parametric distribution. When these assumptions are violated, control charts designed under these assumptions are limited ([crosier1988](@citep), [lowry1992](@citep), [chakraborti2001](@citet), [chakraborti2019](@citet)).

Various nonparametric methods have been developed:
- **Rank-Based Charts**: [qiu2001](@citet), [zou2010](@citet), [li2017a](@citet).
- **Data Categorization**: [qiu2008](@citet), [li2012](@citet), [wang2017](@citet).

These methods help when parametric assumptions are infeasible but might lose effectiveness in information compared to parametric charts ([xie2022](@citep)).

## Selecting Hyperparameters

Choosing the appropriate values of tuning parameters $\bm{\zeta} \in \mathcal{Z} \subseteq \mathbb{R}^{d}$ (e.g., smoothing constant $\lambda$ in EWMA, allowance constant $k$ in CUSUM) is crucial.

### Optimization Problem
The choice of tuning parameters aims to detect specific magnitudes of parameter change efficiently, typically formulated as:
$$
\begin{aligned}
  &\frac{\partial \mathbb{E}_1[\text{RL}]}{\partial \bm{\zeta}}\Big|_{\bm{\zeta}=\bm{\zeta}^*} = \bm{0},\\
    & \text{s.t. } \mathbb{E}_0[\text{RL}] = A_0,
\end{aligned}
$$

Methods like Monte-Carlo simulations might be used when analytical solutions are unavailable, and are the focus of this package.

## Multi-Chart Monitoring Schemes

In complex monitoring scenarios, multiple control charts may be run simultaneously. These schemes are useful for monitoring multiple parameters jointly, such as the mean and variance of a distribution ([gan1995](@citet)).

### Design
The control limits $\bm{LCL}_t = (\text{LCL}_{t1}, \ldots, \text{LCL}_{tJ})$ and $$\bm{UCL}_t = (\text{UCL}_{t1}, \ldots, \text{UCL}_{tJ})$$ are determined so that:
$$
\mathbb{E}_{0}\left[ \min( \text{RL}_1, \text{RL}_2, \ldots, \text{RL}_J) \right] = A_0,
\mathbb{E}_{0}[\text{RL}_1] = \mathbb{E}_{0}[\text{RL}_2] = \ldots = \mathbb{E}_{0}[\text{RL}_J].
$$
This design assumes equal importance of each control chart, but weighting schemes might be used to emphasize the relative importance of each control chart. Algorithms like stochastic approximations are commonly used to solve these designs.
