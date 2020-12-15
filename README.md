# A Software Package of Proximal Gradient Descent Method


Implemented by \
**Zeyu Jin**, PKU\
**Ting Lin**, PKU




The input might be
```
Problem<Vec> prob;
prob.f = f;
prob.gradf = gradf;
prob.h = h;
prob.proxh = h;
prob.mu = mu;
prob.x0 = x0;
prob.setup();
prob.solve();
```