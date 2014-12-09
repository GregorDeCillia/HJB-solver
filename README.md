## HJB-solver

Numerical tool to solve linear [Hamilton Jacobi Bellman Equations](https://en.wikipedia.org/wiki/Hamilton%E2%80%93Jacobi%E2%80%93Bellman_equation). I.e. an equation of the form

![\dot{v}(t,x)=\nabla_xv(t,x)(f_0(t,x)+F(t,x)u), \ v(T,x)=v_T(x)](pictures/HJB.jpg)


It is assumed that the space and the control space are one dimenional.

![f_0,F:\mathbb{R}\times \mathbb{R} \rightarrow \mathbb{R},\ U\subseteq \mathbb{R} ](pictures/f0-f.jpg)