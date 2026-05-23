# Adaptive Rejection Sampling

**By Matthew Seguin**

## Summary

Project Repository URL: [github.com/MMSeguin2003/ars](https://github.com/MMSeguin2003/ars)

The goal of this project is to implement adaptive rejection sampling as described below on a differentiable, log-concave density where $$\frac{dh}{dx} =\frac{d}{dx}\log f(x)$$ is monotonically decreasing.

Adaptive rejection sampling uses lower hull and upper hulls to squeeze a (possibly un-normalized) log-concave density $$f(x)$$ defined in a connected domain $$D$$ (i.e. with support over $$D$$) and sample from it.

<!--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-->

## Background

First we must define what a concave function is.

---

**Definition:**

A function $$\zeta(x)$$ is said to be concave if it satisfies the following condition for all $x,y$ in its domain and $$0<\lambda<1$$: $$\zeta(\lambda x + (1-\lambda) y)\geq\lambda\zeta(x)+(1-\lambda)\zeta(y)$$. It is said to be strictly concave if this inequality is strict for $$x\neq y$$.

---

<!--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-->

Then we can define what a log-concave function is.

---

**Definition:**

A function $$f(x)$$ is said to be log-concave if $$h(x) =\log f(x)$$ is concave. It is said to be strictly log-concave if $$h(x)$$ is strictly concave.

---

This text is only concerned with twice differentiable densities so this is equivalent to saying $$h''(x)\leq 0$$ for all $$x$$ (or using the strict inequality in the case of strict log-concavity). The proof of this statement is in [Appendix I](#i).

Note that this condition is sufficient to say $$h'(x)$$ is strictly decreasing for all $$x$$.

<!--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-->

Now a few important results about log-concave densities.

---

**Theorem:** ([Appendix II](#ii))

Every log-concave density is unimodal. Meaning there exists some $$M\in ℝ$$ such that it is non-decreasing on $$(-\infty,M]$$ and non-increasing on $$[M,\infty)$$.

---

**Theorem:** ([Appendix III](#iii))

Every log-concave density has a finite mean and variance. They also have finite moments of every order.

---

**Theorem:** ([Appendix IV](#iv))

For every unimodal density with finite mean and variance the absolute distance between the mode and mean is at most $$\sqrt{3}$$ standard deviations. In other words $$|\mu-M|\leq\sigma\sqrt{3}$$ where $$\mu$$ is the mean, $$M$$ is the mode, and $$\sigma$$ is the standard deviation.

---

Note that this means every log-concave density is unimodal, has finite moments of every order, and the distance between the mode and mean is at most $$\sigma\sqrt{3}$$.

<!--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-->

## Algorithm

- Find $$h(x) =\log f(x)$$ and $$h'(x) =\frac{d}{dx}\log f(x)$$

- Generate a $$k$$ abscissae $$T_k =$${ $$x_j\in D : j\in$${ $$1, 2, ..., k$$} $$, x_1\leq x_2\leq\dots\leq x_k$$}. If $$D$$ is unbounded on the left choose $$x_1$$ such that $$h'(x_1) > 0$$ and similarly if $$D$$ is unbounded on the right choose $$x_k$$ such that $$h'(x_k) < 0$$

- Find the intersection points for the upper hull

$$z_j =\frac{h(x_{j+1}) - h(x_j) - x_{j+1} h'(x_{j+1}) + x_j h'(x_j)}{h'(x_j) - h'(x_{j+1})} =\frac{(h(x_{j+1}) - x_{j+1}h'(x_{j+1})) - (h(x_j) - x_j h'(x_j))}{h'(x_j) - h'(x_{j+1})}$$

The final expression can be easily calculated in a vectorized manner in python.

- Find the upper hull $$u_k (x)$$ where $$u_k (x) = h(x_j) + (x - x_j)h'(x_j)$$ for $$x\in [z_{j-1}, z_j]$$ for each $$j\in$${ $$1, ..., k$$} and $$z_0 =\inf D$$ (or $$-\infty$$ if $$D$$ is unbounded on the left) and $$z_k =\sup D$$ (or $$\infty$$ if $$D$$ is unbounded on the right)

- Find $$s_k (x) =\frac{\exp u_k (x)}{\int_D\exp u_k (t)\:dt}$$

- Find the lower hull $$l_k (x)$$ where $$l_k (x) =\frac{(x_{j+1} - x)h(x_j) + (x - x_j)h(x_{j+1})}{x_{j+1} - x_j}$$ for $$x\in [x_j, x_{j+1}]$$ for each $$j\in$${ $$1, ..., k-1$$}. Then for $$x < x_1$$ or $$x > x_k$$ we define $$l_k (x) =-\infty$$

- Sample $$x^{(cand)}$$ from $$s_k (x)$$ and independently sample $$w\sim\mbox{Uniform}(0, 1)$$

- If $$w\leq\exp (l_k(x^{(cand)}) - u_k(x^{(cand)}))$$ then accept $$x^{(cand)}$$ and continue from the step where we sample from $$s_k (x)$$ if we don't already have $$n$$ points, otherwise continue to the following step.

- Evaluate $$h(x^{(cand)})$$ and $$h'(x^{(cand)})$$. Then if $$w\leq\exp (h(x^{(cand)}) - u_k(x^{(cand)}))$$ accept $$x^{(cand)}$$, otherwise we reject $$x^{(cand)}$$.

- If $$h(x^{(cand)})$$ and $$h'(x^{(cand)})$$ were evaluated in the sampling process add $$x^{(cand)}$$ to $$T_k$$ to form $$T_{k+1}$$ and repeat the steps starting from after we found the abscissae until $$n$$ points have been sampled.

<!--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-->

## Methodology

### Generating the initial abscissae $$T_k$$

The first thing to note is that generating the initial abscissae is not trivial based on the provided domain:

- If $$D$$ is bounded we can simply take $$k$$ linearly spaced points in $D$ to create the abscissae.
- If $$D$$ is unbounded we want to generate our upper and lower hulls in order to encompass a decent part of the density otherwise they will not serve much purpose in squeezing around the density.

This procedure entails generating the initial abscissae for an unbounded domain $$D$$.

We firstly note that since the (possibly un-normalized) density is log concave it is also unimodal. Therefore there exists some $$x_0$$ such that $$f(x)$$ is monotonically increasing for $$x < x_0$$ and is monotonically decreasing for $$x > x_0$$.

This means in order to encompass a decent part of the density we want to generate $$T_k$$ around the mode. There are infinitely many places the mode could be so we need to have some kind of maximization algorithm, all of which require a starting point. Luckily we know that $$|\mu-M|\leq\sigma\sqrt{3}$$ or equivalently $$\frac{|\mu-M|}{\sigma}\leq\sqrt{3}$$ so we can take our starting point to be $$\mu$$. The mean can easily be found by integration. Even numerical integration which is not exact will suffice since we simply need to be close enough to the mode.

Note that if the domain is unbounded we can not simply start at a point and increment it until we find the function's maximum because we might skip over the whole distribution (even with increments decreasing in size).

Since the function is unimodal (and by assumption the derivative exists) we know the derivative will be positive up until the mode (where it equals 0) then will be negative so finding the mode amounts to finding the zero of the derivative.

Another property that is useful is that

$$h'(x_0) =\left.\frac{d}{dx} h(x)\right|_{x=x_0} =\left.\frac{d}{dx}\log f(x)\right|_{x=x_0} =\left.\frac{f'(x)}{f(x)}\right|_{x=x_0} =\frac{f'(x_0)}{f(x_0)} = 0$$

when $$x_0$$ is the mode. This is useful because the rate of change of the density might be very small at some points so if we try to directly apply an algorithm on $$f'(x)$$ it might think we started at the root so working on the log scale could be important.

After we have the mode we can just take $$k$$ linearly spaced points centered around it for our abscissae $$T_k$$.

### Sampling From $$s_k(x)$$

The next thing to note is that we need to sample $$x^{(cand)}$$ from $$s_k (x)$$ so we need to develop a method of doing that. Such a method is described below:

Clearly $$s_k (x) =\frac{\exp u_k (x)}{\int_D\exp u_k (t)\:dt}$$ is a piecewise function. So we can first sample a categorical variable that tells us which range to choose from. Then based on that we can sample from inside that range.

#### Sampling The Range

Define range $$j$$ to be $$R_j = [z_{j-1}, z_j]$$ for each $$j\in$${ $$1, ..., k$$}. Then

$$\mathbb{P}[x^{(cand)}\in R_j] =\int_{R_j}\frac{\exp u_k (x)}{\int_D\exp u_k (t)\:dt}\:dx$$

First note that $$\bigcup_j R_j = D$$. Now we will evaluate simpler integrals.

$$\int_{R_j}\exp u_k (x)\:dx =\int_{z_{j-1}}^{z_j}\exp (h(x_j) + (x-x_j) h'(x_j))\:dx =\exp (h(x_j) - x_j h'(x_j))\int_{z_{j-1}}^{z_j}\exp (xh'(x_j))\:dx$$

Now there are two cases $$h'(x_j) = 0$$ vs $$h'(x_j)\neq 0$$.

- **First if** $$h'(x_j) = 0$$:

$$\int_{R_j}\exp u_k (x)\:dx =\exp (h(x_j) - x_j h'(x_j))\int_{z_{j-1}}^{z_j}\exp (xh'(x_j))\:dx =\exp h(x_j)\int_{z_{j-1}}^{z_j} 1\:dx = (z_j - z_{j-1})\exp h(x_j)$$

- **Now if** $$h'(x_j)\neq 0$$:

$$\int_{R_j}\exp u_k (x)\:dx =\exp (h(x_j) - x_j h'(x_j))\int_{z_{j-1}}^{z_j}\exp (xh'(x_j))\:dx =\exp (h(x_j) - x_j h'(x_j))\left(\frac{1}{h'(x_j)}\exp (xh'(x_j))\Big{|}_{z_{j-1}}^{z_j}\right)$$

$$=\exp (h(x_j) - x_j h'(x_j))\left(\frac{\exp (z_j h'(x_j)) -\exp (z_{j-1} h'(x_j))}{h'(x_j)}\right)$$

Quickly note that $$\int_D\exp u_k (t)\:dt =\sum_{j=1}^k\int_{R_j}\exp u_k (x)\:dx$$

Therefore if we let $$w_j =\int_{R_j}\exp u_k (x)\:dx$$ we see that:

$$\mathbb{P}[x^{(cand)}\in R_j] =\int_{R_j}\frac{\exp u_k (x)}{\int_D\exp u_k (t)\:dt}\:dx =\frac{1}{\int_D\exp u_k (t)\:dt}\int_{R_j}\exp u_k (x)\:dx$$

$$=\frac{1}{\sum_{m=1}^k\int_{R_m}\exp u_k (x)\:dx}\int_{R_j}\exp u_k (x)\:dx =\frac{w_j}{\sum_{m=1}^k w_m} = p_j$$

Now we know how to sample from a probability vector like this, there are $$k$$ possible outcomes we need to choose one. So we have successfully sampled from the ranges.

#### Sampling Within The Range

Now assuming we are inside of $$R_j$$ we need to sample a value from in that range. The denisty is easily seen to be $$\frac{\exp (h(x_j) + (x-x_j) h'(x_j))}{w_j}$$ since we are just normalizing for already being in $$R_j$$.

Again we have two cases $$h'(x_j) = 0$$ vs $$h'(x_j)\neq 0$$.

- **First if** $$h'(x_j) = 0$$:

If $$h'(x_j) = 0$$ this is clearly a uniform random variable over $$[z_{j-1}, z_j]$$ as the density does not depend on $$x$$ so we can just sample a uniform.

- **Now if** $$h'(x_j)\neq 0$$:

If $$h'(x_j)\neq 0$$ we get the following CDF for $$x^{(cand)} | x^{(cand)}\in R_j$$:

$$F(x) =\int_{z_{j-1}}^x \frac{\exp (h(x_j) + (x-x_j) h'(x_j))}{w_j}\:dx =\frac{\exp (h(x_j) - x_j h'(x_j))}{w_j}\int_{z_{j-1}}^x\exp (xh'(x_j))\:dx$$

$$=\frac{\exp (h(x_j) - x_j h'(x_j))}{w_j}\left(\frac{1}{h'(x_j)}\exp (xh'(x_j))\Big{|}_{z_{j-1}}^x\right)$$

$$=\frac{1}{w_j}\left(\frac{\exp (h(x_j) - x_j h'(x_j))}{h'(x_j)}\right)\left(\exp (xh'(x_j)) -\exp (z_{j-1} h'(x_j))\right)$$

Quickly note that

$$\frac{1}{w_j} =\left(\frac{1}{\exp (z_j h'(x_j)) -\exp (z_{j-1} h'(x_j))}\right)\left(\frac{h'(x_j)}{\exp (h(x_j) - x_j h'(x_j))}\right)$$

Continuing we get

$$F(x)=\frac{1}{w_j} =\left(\frac{1}{\exp (z_j h'(x_j)) -\exp (z_{j-1} h'(x_j))}\right)\left(\frac{h'(x_j)}{\exp (h(x_j) - x_j h'(x_j))}\right)$$

$$=\left(\frac{1}{\exp (z_j h'(x_j)) -\exp (z_{j-1} h'(x_j))}\right)\left(\exp (xh'(x_j)) -\exp (z_{j-1} h'(x_j))\right) =\frac{\exp (xh'(x_j)) -\exp (z_{j-1} h'(x_j))}{\exp (z_j h'(x_j)) -\exp (z_{j-1} h'(x_j))}$$

Since the term in the $$w_j$$ cancels out to simplify to the above. We can then use the inverse CDF method to sample from this. Namely sample $$u\sim\mbox{Uniform}(0, 1)$$ then our sample for $$x^{(cand)}$$ is $$F^{-1} (u)$$. This is solved below:

If $$F^{-1} (u) = x^{(cand)}$$ then $$u = F(x^{(cand)})$$ so

$$u = F(x^{(cand)}) =\frac{\exp (x^{(cand)}h'(x_j)) -\exp (z_{j-1} h'(x_j))}{\exp (z_j h'(x_j)) -\exp (z_{j-1} h'(x_j))}$$

$$u\left(\exp (z_j h'(x_j)) -\exp (z_{j-1} h'(x_j))\right) =\exp (x^{(cand)}h'(x_j)) -\exp (z_{j-1} h'(x_j))$$

$$u\left(\exp (z_j h'(x_j)) -\exp (z_{j-1} h'(x_j))\right) +\exp (z_{j-1} h'(x_j)) =\exp (x^{(cand)}h'(x_j))$$

$$\log\left(u\left(\exp (z_j h'(x_j)) -\exp (z_{j-1} h'(x_j))\right) +\exp (z_{j-1} h'(x_j))\right) = x^{(cand)}h'(x_j)$$

$$x^{(cand)} =\frac{\log\left(u\left(\exp (z_j h'(x_j)) -\exp (z_{j-1} h'(x_j))\right) +\exp (z_{j-1} h'(x_j))\right)}{h'(x_j)}$$

$$=\frac{\log\left(u\exp (z_j h'(x_j)) + (1-u)\exp (z_{j-1} h'(x_j))\right)}{h'(x_j)}$$

Quickly note that the cases of $$z_0 = -\infty$$ and $$z_k =\infty$$ are not an issue for these integrals converging or the probabilities being positive. Since in these cases we chose $$x_1$$ such that $$h'(x_1) > 0$$ and similarly $$x_k$$ such that $$h'(x_k) < 0$$. This implies that $$z_0 h'(x_1) = -\infty$$ and $$z_k h'(x_k) = -\infty$$ which tells us:

$$w_1 =\exp (h(x_1) - x_1 h'(x_1))\left(\frac{\exp (z_1 h'(x_1)) -\exp (z_0 h'(x_1))}{h'(x_1)}\right)$$

$$=\frac{\exp (h(x_1) - x_1 h'(x_1))}{h'(x_1)}\left(\exp (z_1 h'(x_1)) -\exp (z_0 h'(x_1))\right) =\frac{\exp (h(x_1) - x_1 h'(x_1))}{h'(x_1)}\left(\exp (z_1 h'(x_1)) -\exp (-\infty)\right)$$

$$=\frac{\exp (h(x_1) - x_1 h'(x_1))}{h'(x_1)}\left(\exp (z_1 h'(x_1))\right) > 0$$

Because $$h'(x_1) > 0$$. Similarly:

$$w_k =\exp (h(x_k) - x_k h'(x_k))\left(\frac{\exp (z_k h'(x_k)) -\exp (z_{k-1} h'(x_k))}{h'(x_k)}\right)$$

$$=\frac{\exp (h(x_k) - x_k h'(x_k))}{h'(x_k)}\left(\exp (z_k h'(x_k)) -\exp (z_{k-1} h'(x_k))\right)$$

$$=\frac{\exp (h(x_k) - x_k h'(x_k))}{h'(x_k)}\left(\exp (-\infty) -\exp (z_{k-1} h'(x_k))\right)=\frac{\exp (h(x_k) - x_k h'(x_k))}{h'(x_k)}\left(-\exp (z_{k-1} h'(x_k))\right) > 0$$

Because $$h'(x_k) < 0$$

<!--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-->

## Appendix

<!--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-->

### I

**Theorem:**

A twice differentiable function $h$ is concave if and only if $$h''(x)\leq 0$$ for all $x$ in its domain. It is strictly concave if and only if this inequality is strict.

**Proof:**

*Case 1* - *Non-strict case*

- First we show that $$h''(x)\leq 0$$ for all $$x$$ implies that $$h$$ is concave.

Suppose $$h''(x)\leq 0$$ for all $$x$$ and let $$0<\lambda<1$$ be arbitrary. This means that if $x<y$ then $$h'(x)\geq h'(y)$$. Firstly, if $$x = y$$ then we know trivially that

$$h(\lambda x+(1-\lambda)y)=h(\lambda x+(1-\lambda)x)=h(x)\geq h(x)=\lambda h(x)+(1-\lambda)h(x)=\lambda h(x)+(1-\lambda)h(y)$$

Otherwise we can say without loss of generality that $$x<y$$ and hence $$h'(x)\geq h'(y)$$. Let $$z=\lambda x+(1-\lambda)y$$. Then $$z$$ is a weighted average of $$x$$ and $$y$$ so $$x < z < y$$. We want to show that $$h(z)\geq\lambda h(x)+(1-\lambda)h(y)$$ which is equivalent to showing that

$$0\leq h(z)-\lambda h(x)-(1-\lambda)h(y)=\lambda(h(z)-h(x))+(1-\lambda)(h(z)-h(y))$$

Now by the mean value theorem there exists a $$\xi_1\in(x,z)$$ and $$\xi_2\in(z,y)$$ such that $$h(z)-h(x)=h'(\xi_1)(z-x)$$ and $$h(z)-h(y)=h'(\xi_2)(z-y)$$. Therefore we may rewrite the previous condition as

$$0\leq\lambda(h(z)-h(x))+(1-\lambda)(h(z)-h(y))=\lambda h'(\xi_1)(z-x)+(1-\lambda)h'(\xi_2)(z-y)$$

Now we note that $$z-x=\lambda x+(1-\lambda)y-x=(1-\lambda)(y-x)$$ and $$z-y=\lambda x+(1-\lambda)y-y=-\lambda(y-x)$$. Therefore the condition becomes

$$0\leq\lambda h'(\xi_1)(z-x)+(1-\lambda)h'(\xi_2)(z-y)=\lambda(1-\lambda)(y-x)h'(\xi_1)-\lambda(1-\lambda)h'(\xi_2)(y-x)$$

$$=\lambda(1-\lambda)(y-x)(h'(\xi_1)-h'(\xi_2))$$

Since $$0<\lambda<1$$ we know $$\lambda(1-\lambda)>0$$. Then we already know $$y-x>0$$. And finally since $$x<\xi_1<z<\xi_2<y$$ we know that $$h'(\xi_1)\geq h'(\xi_2)$$ which shows that the desired inequality is valid. Since $$\lambda$$, $$x$$, and $$y$$ were arbitrary this shows that $$h$$ is concave by definition.

- Now we show that if $$h$$ is concave and twice differentiable then $$h''(x)\leq 0$$ for all $$x$$.

Let $$\lambda =\frac{1}{2}$$, $$t\neq 0$$, $$x=w-t$$, and $$y=w+t$$. Firstly we know $$\lambda x+(1-\lambda)y=\frac{1}{2}(w-t)+\frac{1}{2}(w+t)=w$$. Then using the definition of concavity we know that

$$h(w)=h(\lambda x+(1-\lambda)y)\geq\lambda h(x)+(1-\lambda)h(y)=\frac{1}{2}h(w-t)+\frac{1}{2}h(w+t)$$

Which is equivalent to $$h(w-t)-2h(w)+h(w+t)\leq 0$$. Dividing both sides by $$t^2$$ we get that

$$\frac{h(w-t)-2h(w)+h(w+t)}{t^2}\leq 0$$

Therefore the limit is also true

$$\lim_{t\rightarrow 0}\frac{h(w-t)-2h(w)+h(w+t)}{t^2}\leq 0$$

Now the limit of the numerator is

$$\lim_{t\rightarrow 0} h(w-t)-2h(w)+h(w+t)=h(w)-2h(w)+h(w)=0$$

Since $$h$$ is continuous. Clearly the limit of the denominator is $$\lim_{t\rightarrow 0}t^2=0$$. So we may use L'Hopital's rule

$$\lim_{t\rightarrow 0}\frac{h(w-t)-2h(w)+h(w+t)}{t^2}=\lim_{t\rightarrow 0}\frac{-h'(w-t)+h'(w+t)}{2t}=\lim_{t\rightarrow 0}\frac{(h'(w)-h'(w-t))+(h'(w+t)-h'(w))}{2t}$$

$$=\frac{1}{2}\left(\lim_{t\rightarrow 0}\frac{h'(w)-h'(w-t)}{t}+\frac{h'(w+t)-h'(w)}{t}\right)=\frac{1}{2}\left(\lim_{t\rightarrow 0}\frac{h'(w)-h'(w-t)}{t}+\lim_{t\rightarrow 0}\frac{h'(w+t)-h'(w)}{t}\right)$$

$$=\frac{1}{2}\left(\lim_{t\rightarrow 0}\frac{h'(w-t)-h'(w)}{-t}+h''(w)\right)=\frac{1}{2}\left(\lim_{s\rightarrow 0}\frac{h'(w+s)-h'(w)}{s}+h''(w)\right)=\frac{1}{2}\left(h''(w)+h''(w)\right)=h''(w)$$

By using $$s=-t$$ and the definition of the derivative. Therefore we have shown that $$h''(w)\leq 0$$ for an arbitrary $$w$$.

Therefore if $$h$$ is a twice differentiable function then $$h$$ is concave if and only if $$h''(x)\leq 0$$ in its domain $$\blacksquare$$

*Case 2* - *Strict case*

- First we show that $$h''(x)<0$$ for all $x$ implies that $$h$$ is strictly concave.

Using the same path of logic as the non-strict case we now get $$h'(\xi_1)>h'(\xi_2)$$ and hence

$$0<\lambda(1-\lambda)(y-x)(h'(\xi_1)-h'(\xi_2))$$

Which implies $$h(\lambda x+(1-\lambda)y)>\lambda h(x)+(1-\lambda)h(y)$$ for $$x\neq y$$, proving strict concavity.

- Now we show that if $$h$$ is strictly concave and twice differentiable then $$h''(x)<0$$ for all $$x$$.

Using the same path of logic as the non-strict case we now get

$$h(w)=h(\lambda x+(1-\lambda)y)>\lambda h(x)+(1-\lambda)h(y)=\frac{1}{2}h(w-t)+\frac{1}{2}h(w+t)$$

Which implies $$h(w-t)-2h(w)+h(w+t)<0$$ and $$\frac{h(w-t)-2h(w)+h(w+t)}{t^2}<0$$. Again this inequality holds in the limit and since the limit is just $$h''(w)$$ as before we have shown that $$h''(w)<0$$ for an arbitrary $w$ in its domain.

Therefore if $$h$$ is a twice differentiable function then $$h$$ is strictly concave if and only if $$h''(x)<0$$ in its domain $$\blacksquare$$

<!--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-->

### II

**Theorem:**

Every log-concave density is unimodal. Meaning there exists some $$M\in ℝ$$ such that it is non-decreasing on $$(-\infty,M]$$ and non-increasing on $$[M,\infty)$$.

**Proof:**

Let $$f(x)$$ be log concave and $$h(x)=\log f(x)$$ then by definition $$h(x)$$ is concave. First let $$x<z<y$$. Then there exists some $$\lambda>0$$ such that $$z=\lambda x+(1-\lambda)y$$. Therefore we know that $$h(z)=h(\lambda x+(1-\lambda)y)\geq\lambda h(x)+(1-\lambda)h(y)$$ by the concavity of $$h$$. Then since $$y-x>0$$ we can rewrite this as

$$(y-x)h(z)\geq(y-x)(\lambda h(x)+(1-\lambda)h(y))$$

Then we can rewrite this as

$$(y-x)(\lambda h(z)+(1-\lambda)h(z))\geq(y-x)(\lambda h(x)+(1-\lambda)h(y))$$

$$\lambda(y-x)(h(z)-h(x))\geq(1-\lambda)(y-x)(h(y)-h(z))$$

Then note that $$y-z=y-\lambda x-(1-\lambda)y=\lambda(y-x)$$ and $$z-x=\lambda x+(1-\lambda)y-x=(1-\lambda)(y-x)$$ so

$$(y-z)(h(z)-h(x))=\lambda(y-x)(h(z)-h(x))\geq(1-\lambda)(y-x)(h(y)-h(z))=(z-x)(h(y)-h(z))$$

Which tells us that

$$\frac{h(z)-h(x)}{z-x}\geq\frac{h(y)-h(z)}{y-z}$$

In other words, $$h$$ has decreasing slopes. This is trivial if we assumed $$h$$ was twice differentiable but this proof is more general.

Now assume for the sake of contradiction that there exists some $$x_1<x_2<x_3$$ such that $$h(x_2)<h(x_1)$$ and $$h(x_3)>h(x_2)$$, meaning $$h$$ decreases from $$x_1$$ to $$x_2$$ then increases from $$x_2$$ to $$x_3$$. Then

$$\frac{h(x_2)-h(x_1)}{x_2-x_1}<0\hspace{1cm}\text{and}\hspace{1cm}\frac{h(x_3)-h(x_2)}{x_3-x_2}>0$$

Which is a contradictions since $$h$$ must have decreasing slopes. Therefore once $$h$$ begins decreasing it must never increase again.

Now let $$M=\sup$${ $$x\in ℝ:h\text{ is non-decreasing on }(-\infty,x]$$}. Such a value exists since $$\int_{ℝ}e^{h(x)}dx=1$$. If no such $M$ exists then $$h(x)$$ and hence $$f(x)=e^{h(x)}$$ is always non-decreasing on ℝ which violates the fact that the integral is finite. Note that this implies that $$h(x)$$ decreases after $$M$$ and is then non-increasing on $$[M,\infty)$$ by the previous result.

Finally since $$e^x$$ is an increasing function in $$x$$ we know that $$f(x)=e^{h(x)}$$ increases when $$h(x)$$ increases and decreases when $$h(x)$$ decreases.

Therefore there exists some $$M\in ℝ$$ such that $$f(x)$$ is non-decreasing on $$(-\infty,M]$$ and non-increasing on $$[M,\infty)$$. Showing that $$f$$ is unimodal by definition $$\blacksquare$$

<!--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-->

### III

**Theorem:**

Every log-concave density has a finite mean and variance. They also have finite moments of every order.

**Proof:**

Let $$f$$ be a log-concave density. It suffices to show that the $$k\text{th}$$ moment is finite for an arbitrary $$k\in$${ $$0,1,2,\dots$$}.

*Case 1* - *Finite domain* $$D$$

If its domain $$D$$ is finite then clearly

$$\int_D |x|^kf(x)dx=\int_{\inf D}^{\sup D} |x|^kf(x)dx<\infty$$

Since it is an integral of finite values on a finite domain.

*Case 2* - *Unbounded domain* $$D$$ *on right side*

Let $$m=\sup$${ $$x\in ℝ:f\text{ is non-decreasing on }(-\infty,x]$$}. Define $$h(x)=\log f(x)$$, then by definition $$h$$ is concave. Therefore the slope

$$\frac{h(x)-h(m)}{x-m}$$

is non-increasing and negative for $$x>m$$. So there exists a $$c>0$$ such that for all $$x>m$$

$$\frac{h(x)-h(m)}{x-m}\leq -c\implies h(x)\leq h(m)-c(x-m)$$

Therefore we see

$$\int_D|x|^kf(x)dx=\int_{\inf D}^m|x|^kf(x)dx+\int_m^\infty|x|^ke^{h(x)}dx\leq\int_{\inf D}^m|x|^kf(x)dx+\int_m^\infty|x|^ke^{h(m)-c(x-m)}dx$$

$$=\int_{\inf D}^m|x|^kf(x)dx+e^{h(m)+cm}\int_m^\infty|x|^ke^{-cx}dx<\infty$$

Since the exponential decay guarantees convergence for a fixed $$k$$.

*Case 3* - *Unbounded domain* $$D$$ *on left side*

Let $$m=\sup$${ $$x\in ℝ:f\text{ is non-decreasing on }(-\infty,x]$$}. Define $$h(x)=\log f(x)$$, then by definition $$h$$ is concave. Therefore the slope

$$\frac{h(m)-h(x)}{m-x}$$

is non-decreasing and positive for $$x<m$$. So there exists a $$c'>0$$ such that for all $$x<m$$

$$\frac{h(m)-h(x)}{m-x}\leq c'\implies h(x)\leq h(m)-c'(m-x)$$

Therefore we see

$$\int_D|x|^kf(x)dx=\int_{-\infty}^m|x|^ke^{h(x)}dx+\int_m^{\sup D}|x|^kf(x)dx\leq\int_{-\infty}^m|x|^ke^{h(m)-c'(m-x)}dx+\int_m^{\sup D}|x|^kf(x)dx$$

$$=e^{h(m)-c'm}\int_{-\infty}^m|x|^ke^{c'x}dx+\int_m^{\sup D}|x|^kf(x)dx<\infty$$

Since the exponential decay guarantees convergence for a fixed $$k$$.

*Case 4* - *Unbounded domain* $$D$$ *on both side*

Using the same definitions from Case 2 and Case 3 we see that

$$\int_D|x|^kf(x)dx=\int_{-\infty}^m|x|^ke^{h(x)}dx+\int_m^\infty|x|^ke^{h(x)}dx\leq\int_{-\infty}^m|x|^ke^{h(m)-c'(m-x)}dx+\int_m^\infty|x|^ke^{h(m)-c(x-m)}dx<\infty$$

Since the exponential decay guarantees convergence for a fixed $$k$$.

Therefore every log-concave density has finite mean, variance, and moments of every order $$\blacksquare$$

<!--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-->

### IV

**Theorem:**

For every unimodal density with finite mean and variance the absolute distance between the mode and mean is at most $$\sqrt{3}$$ standard deviations. In other words $$|\mu-M|\leq\sigma\sqrt{3}$$ where $$\mu$$ is the mean, $$M$$ is the mode, and $$\sigma$$ is the standard deviation.

**Proof:**

Without loss of generality consider a unimodal density $$f$$ with mode $$M$$, mean $$\mu\geq M$$, and variance $$\sigma^2$$. If $$\mu < M$$ we may simply consider the density $$f(-x)$$ which has mode $$-M$$, mean $$-\mu >-M$$, and variance $$\sigma^2$$ since the quantity $$|(-\mu)-(-M)|=|\mu-M|$$ remains unchanged.

*Case 1* - $$\mu=M$$

Trivially $$|\mu-M|=0\leq\sigma\sqrt{3}$$

*Case 2* - $$\mu >M$$ *and mass to the left of* $$M$$

If $$f(x)>0$$ for some $$x<M$$ then by redistributing the mass to $$x>M$$ we can increase the mean while keeping the mode unchanged so we have increased $$|\mu-M|$$. Formally, define the density $$\phi$$ for $$x\geq M$$ to be

$$\phi(x)=f(x)\left(\int_M^\infty f(x)dx\right)^{-1}$$

Then

$$\int_{-\infty}^\infty\phi(x)dx=\int_M^\infty\phi(x)dx=\int_M^\infty f(x)\left(\int_M^\infty f(x)dx\right)^{-1}dx=\left(\int_M^\infty f(x)dx\right)^{-1}\left(\int_M^\infty f(x)dx\right)=1$$

So $$\phi$$ is indeed a density. It is also still unimodal, with the same mode, since we are simply rescaling the density starting at the mode. Additionally we see that

$$\mu=\int_{-\infty}^\infty xf(x)dx=\int_{-\infty}^Mxf(x)dx+\int_M^\infty xf(x)dx\leq M\int_{-\infty}^Mf(x)dx+\int_M^\infty xf(x)dx$$

$$=M\left(1-\int_M^\infty f(x)dx\right)+\int_M^\infty xf(x)dx$$

Now as a quick note

$$M\int_M^\infty f(x)dx\leq\int_M^\infty xf(x)dx\hspace{0.25cm}\implies\hspace{0.25cm}M\left(1-\int_M^\infty f(x)dx\right)\left(\int_M^\infty f(x)dx\right)\leq\left(1-\int_M^\infty f(x)dx\right)\left(\int_M^\infty xf(x)\right)$$

Which implies

$$M\left(1-\int_M^\infty f(x)dx\right)\leq\left(\int_M^\infty xf(x)\left(\int_M^\infty f(x)dx\right)^{-1}dx\right)-\int_M^\infty xf(x)dx$$

Giving us the result

$$\mu\leq M\left(1-\int_M^\infty f(x)dx\right)+\int_M^\infty xf(x)dx\leq\int_M^\infty xf(x)\left(\int_M^\infty f(x)dx\right)^{-1}dx=\int_{-\infty}^\infty x\phi(x)dx=\mu'$$

Therefore the mean has increased and since $$\mu>M$$ we see that $$|\mu-M|=\mu-M\leq\mu'-M=|\mu'-M|$$, i.e. the distance between the mean and the mode has increased. There are other ways to redistribute the mass to increase the mean but this shows that the extreme case is one where $$\mu>M$$ and $$M$$ is the minimum of the domain. Therefore the distance between the mean and the mode in this case is at most the upper bound derived in Case 3.

*Case 3* - $$M$$ *is the minimum of the domain*

Let $$\zeta(x)=f(x+M)$$ for $$x\geq 0$$ then

$$\int_{-\infty}^\infty\zeta(x)dx=\int_0^\infty\zeta(x)dx=\int_0^\infty f(x+M)dx\int_M^\infty f(w)dw=1$$

So $$\zeta$$ is a density. Additionally, $$\arg\max_x\zeta(x)=\arg\max_xf(x+M)=0$$ is the mode of $$\zeta$$ since $$M=\arg\max_xf(x)$$ is the mode of $f$. Then

$$\int_{-\infty}^\infty x\zeta(x)dx=\int_0^\infty x\zeta(x)dx=\int_0^\infty xf(x+M)dx=\int_M^\infty (w-M)f(w)dw$$

$$=\int_M^\infty wf(w)dw-\int_M^\infty Mf(w)dw=\mu-M=\delta$$

is the mean of $$\zeta$$. And finally

$$\int_{-\infty}^\infty(x-\delta)\zeta(x)dx=\int_0^\infty(x-(\mu-M))\zeta(x)dx=\int_0^\infty(x-(\mu-M))f(x+M)dx$$

$$=\int_M^\infty((w-M)-(\mu-M))f(w)dw=\int_M^\infty(w-\mu)f(w)dw=\sigma^2$$

Is the variance of $$\zeta$$.

Let $$g(y)=\sup$${ $$x\geq 0:\zeta(x)\geq y$$} then clearly $$g$$ is non-increasing and by Fubini's Theorem

$$1=\int_0^\infty\zeta(x)dx=\int_0^\infty\int_0^{\zeta(x)}dy\:dx=\int_0^{\zeta(0)}\int_0^{g(y)}dx\:dy=\int_0^{\zeta(0)}g(y)dy$$

$$\delta=\int_0^\infty x\zeta(x)dx=\int_0^\infty\int_0^{\zeta(x)}x\:dy\:dx=\int_0^{\zeta(0)}\int_0^{g(y)}x\:dx\:dy=\int_0^{\zeta(0)}\frac{(g(y))^2}{2}dy$$

$$\delta^2+\sigma^2=\int_0^\infty x^2\zeta(x)dx=\int_0^\infty\int_0^{\zeta(x)}x^2\:dy\:dx=\int_0^{\zeta(0)}\int_0^{g(y)}x^2\:dx\:dy=\int_0^{\zeta(0)}\frac{(g(y))^3}{3}dy$$

So we want to maximize

$$|\mu-M|=\mu-M=\delta=\int_0^{\zeta(0)}\frac{(g(y))^2}{2}dy$$

Subject to the constraints

$$1=\int_0^{\zeta(0)}g(y)\:dy\hspace{3cm}\delta^2+\sigma^2=\int_0^{\zeta(0)}\frac{(g(y))^3}{3}dy$$

We know by Hölder's Inequality with $$p=q=2$$ that

$$\int_0^{\zeta(0)}(g(y))^2dy=\int_0^{\zeta(0)}\left|(g(y))^{1/2}(g(y))^{3/2}\right|dy=\lVert g^{1/2}g^{3/2}\rVert_1\leq\lVert g^{1/2}\rVert_p\lVert g^{3/2}\rVert_q$$

$$=\left(\int_0^{\zeta(0)}\left|(g(y))^{1/2}\right|^pdy\right)^{1/p}\left(\int_0^{\zeta(0)}\left|(g(y))^{3/2}\right|^qdy\right)^{1/q}=\left(\int_0^{\zeta(0)}g(y)dy\right)^{1/2}\left(\int_0^{\zeta(0)}(g(y))^3dy\right)^{1/2}$$

Then squaring both sides we get

$$\left(\int_0^{\zeta(0)}(g(y))^2dy\right)^2\leq\left(\int_0^{\zeta(0)}g(y)dy\right)\left(\int_0^{\zeta(0)}(g(y))^3dy\right)=3(\delta^2+\sigma^2)$$

Therefore by dividing both sides by 4 we see

$$\delta^2=\left(\int_0^{\zeta(0)}\frac{(g(y))^2}{2}dy\right)^2\leq\frac{3}{4}(\delta^2+\sigma^2)$$

Simplifying we get $$\delta^2\leq 3\sigma^2$$ which implies $$|\mu-M|=\mu-M=\delta\leq\sigma\sqrt{3}$$.

This proof covers all possible unimodal densities with finite mean, and variance. Therefore we have shown that $$|\mu-M|\leq\sigma\sqrt{3}$$ for any unimodal density $$\blacksquare$$

---

As an extra note

The inequality becomes an equality if and only if $$g=\left|g^{1/2}\right|^2=c\left|g^{3/2}\right|^2=g^3$$. This would mean $$g(y)(cg(y)^2-1)=0$$. We know $$g(y)$$ is not identically 0 since we have assumed the mean of our distribution $$\mu'>0$$ so there is some mass greater than 0. Therefore the equality holds if and only if $$g(y)=\frac{1}{\sqrt{c}}$$ i.e. $$g(y)$$ is a constant. Note that this represents a uniform distribution since it means $$\zeta(x)$$ stays at the same level on its domain.

The variance of a Uniform($$a$$, $$b$$) distribution is $$\frac{1}{12}(b-a)^2$$. Here we are setting $$a=0$$ and constraining variance to $$\sigma^2$$ so $$\frac{1}{12}b^2=\sigma^2$$ and hence $$b=\sigma\sqrt{12}$$. The density is then $$\zeta(x)=\frac{1}{\sigma\sqrt{12}}$$ for $$x\in[0,\sigma\sqrt{12}]$$ and 0 elsewhere. Therefore $$g(y)=\sup$${ $$x\geq 0:\zeta(x)\geq y$$} $$=\sigma\sqrt{12}$$ for $$y\leq\frac{1}{\sigma\sqrt{12}}=\zeta(0)=f(M)$$. Checking that the conditions are met we see

$$\int_0^\frac{1}{\sigma\sqrt{12}}g(y)dy=\int_0^\frac{1}{\sigma\sqrt{12}}\sigma\sqrt{12}\:dy=1$$

$$\int_0^\frac{1}{\sigma\sqrt{12}}\frac{(g(y))^2}{2}dy=\int_0^\frac{1}{\sigma\sqrt{12}}6\sigma^2\:dy=\frac{6\sigma^2}{\sigma\sqrt{12}}=\frac{3\sigma}{\sqrt{3}}=\sigma\sqrt{3}=\mu'=\mu-M=|\mu-M|$$

$$\int_0^\frac{1}{\sigma\sqrt{12}}\frac{(g(y))^3}{3}dy=\int_0^\frac{1}{\sigma\sqrt{12}}4\sigma^3\sqrt{12}\:dy=4\sigma^2=(\sigma\sqrt{3})^2+\sigma^2=(\mu')^2+\sigma^2$$

The same result holds for any uniform distribution as well, for simplicity this case assumed the left bound was 0.

<!--
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
-->