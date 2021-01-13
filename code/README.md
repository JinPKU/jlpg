# Basic Framework of JLPG
We first define a **problem**, and set up the **options**, and then use
```C++
pgm(problem,initialvalue,option);
```
to solve the problem.

<!-- $$\newcommand{\prox}{\textbf{prox}}$$ -->
## Problem Construction
The problem consists of the following three part:
$$ p(x) = f(x) + \mu h(x)$$
The proximal gradient method updates by
$$ x^+ = \textbf{prox}_{\lambda\mu h}(x - \lambda\nabla f(x))$$
each step. 
Hence we need the value of $f$, and its gradient function; we need the value of $h$, and its proximal operator. 

### `grad_pair`
The pair $(f,\nabla f)$ is represented by the template class `grad_pair<T>`. It contains the following member, `f` and `proxf`. Here we require `f(x)` returns $f(x)$, and `gradf(x)` returns $\nabla f(x)$. 

### `prox_pair`
The pair $(h, \textbf{prox}_{th}(x))$ is represented by the template class `prox_pair<T>`. It contains the following member, `h` and `gradh`. Here we require `h(x)` returns $h(x)$ and $proxh(x,t)$ returns $\textbf{prox}_{th}(x)$.

