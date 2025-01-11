[![Open in Codespaces](https://classroom.github.com/assets/launch-codespace-7f7980b617ed060a017424585567c406b6ee15c891e84e1186181d67ecf80aa0.svg)](https://classroom.github.com/open-in-codespaces?assignment_repo_id=13065509)

<h1>Implementation, Parallelisation and Evaluation of N Body Methods</h1>
<p>The N-Body problem is a fundamental astrophysics problem concerning how to model the
motion of a number, <b>N</b>, different objects that interact with each other gravitationally, where <b>N</b>
is greater than or equal to three. An analytical solution for the Two-Body Problem (<b>N</b>=2) was first solved by Isaac Newton; however, it is currently
impossible to produce an analytical solution to the N-Body problem, due to the vast amount
of unknown variables introduced by modelling a system with <b>N</b> bodies. As a result of this,
finding a solution to the N-Body problem is an area of active research and numerical
methods have been created to approximate the solution. Numerical solutions to model the
N-Body problem are used in many areas such as the modelling of galaxy formations, stars,
and orbits. This project aims to explore various algorithmic solutions to the N-Body problem,
followed by the development of software that implements and compares these different
methods. Parallel computing and optimisation techniques will also be explored in order to
gauge the improvements that can gained over the serial/unoptimised versions of these
methods.</p>

<h2>Summary</h2>
<p>This project aims to investigate the N-Body problem by developing, parallelising, and
evaluating different computational methods that model the problem. Initially, the project
began by implementing the serial versions of the following N-Body methods: the Brute Force
method, the ‘Barnes-Hut’ algorithm and the Fast Multipole Method. Following the completion
of the implementation of these methods, the parallel versions were produced.</p>

<p>The first stage of development focused on implementing the Brute Force method as this is
the intuitive solution. This was subsequently followed by the ‘Barnes-Hut’ and the Fast
Multipole Method. In order to facilitate the implementation of the latter two methods, a
quadtree data structure was also implemented, which is responsible for handling the
decomposition of the simulation space. Parallel and optimised versions of the algorithms
were then implemented to evaluate the performance of these methods against their
unoptimized counterparts.</p>

<p>The software used to model these methods was rigorously tested using a testing strategy
that comprised of unit testing and integration testing. Further testing was then performed in
order to produce quantitative data that could be used to assess the execution time and
accuracy of each method with one another.</p>

<h2>Usage</h2>
<b>
<p>make all</p>
<p>Usage: ./simulation brute-force | barnes-hut | fast-multipole | parallel-brute-force | parallelbarnes-
hut | parallel-fast-multipole [numbers-of-bodies-simulated]</p>
</b>
