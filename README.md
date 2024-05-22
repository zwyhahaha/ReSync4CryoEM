# Orientation Determination of Cryo-EM Images Using Block Stochastic Riemannian Subgradient Methods

## Disclaimer

The experiment frameowrk is downloaded from

-  https://github.com/wenyouwei/2022cryoem which is associated with the paper "Pan, Huan, et al. "Orientation estimation of cryo-EM images using projected gradient descent method." *Inverse Problems* 39.4 (2023): 045002."
- The majority code in Utils is obtained in ASPIRE Toolkit: https://spr.math.princeton.edu/
- Implementation of ReSync algorithm is obtained from https://github.com/Huikang2019/ReSync , which is associated with the paper "Liu, Huikang, Xiao Li, and Anthony Man-Cho So. "ReSync: Riemannian Subgradient-based Robust Rotation Synchronization." *Advances in Neural Information Processing Systems* 36 (2024)."

We upload these dependencies in order to be self-contained. In this experiment framework, we add ReSync-norm as well as stochastic ReSync algorithms (SGD, BCD, BSGD) in the `Algorithm` folder.

To accelerate ReSync algorithm, we provide C-implementation for gradient estimation.

## Usage

(1) run the following command in the command line of Matlab, and you can find `Utils` folder in our supplementary material

```matlab
initpath
```

(2) The code for reproduce our results can be found in `Runners` folder.

- `run_synthetic.m` is the code for Table 1.
- `run_real.m` is the code for Table 2/3.
- `run_plot_benchmark.m` is the code for Figure 2.
- `run_plot_BSGD.m` is the code for Figure 3.

## Issues

If you have encountered with some `function or file not found` error, you can first check `initpath` command. If the path setting is ok, then it is probably due to you have not complied some code in C. To address this issue, you can find that file (e.g. `Q_theta.c`), and run the following in Matlab command line:

```matlab
mex Q_theta.c
```

