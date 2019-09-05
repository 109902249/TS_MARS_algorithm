# MATLAB implementation of the TS-MARS algorithm

0. Corresponding author: Qi Zhang, Department of Applied Mathematics and Statistics, Stony Brook University, Stony Brook, NY 11794-3600

1. The TS-MARS algorithm by Qi Zhang and Jiaqiao Hu [1] is implemented for solving single-objective box-constrained expensive deterministic optimization problems.

2. TS-MARS a random search algorithm for solving deterministic optimization problems in a black-box scenario. The algorithm has a model-based nature and iteratively finds improved solutions by modifying and sampling from a probability distribution over the solution space. In contrast to existing algorithms in the class, which are mostly population-based, our approach employs a two-time-scale stochastic approximation idea and uses only a single candidate solution per iteration. We prove global convergence of the algorithm and carry out numerical experiments to illustrate its performance.

3. In this implementation, the algorithm samples candidate solutions from a sequence of independent multivariate normal distributions that recursively  approximiates the corresponding Boltzmann distributions [2].

### Reference:
1. Qi Zhang and Jiaqiao Hu (2017): A Two-time-scale Adaptive Search Algorithm for Global Optimization. The Proceedings of the 2017 Winter Simulation Conference, pp. 2069-2079.
2. Jiaqiao Hu and Ping Hu (2011): Annealing adaptive search, cross-entropy, and stochastic approximation in global optimization. Naval Research Logistics 58(5):457-477.
