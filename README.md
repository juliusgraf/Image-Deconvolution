# Image Deconvolution

## Overview of the project

The goal is to estimate a sharp image $\widehat{x}[m,n]$ from a blurred image $y[m,n]$, which is obtained by convolving the sharp image $x[m,n]$ with the impulse response $h[p,q]$, followed by corruption with additive noise $\epsilon[m,n]$. This is achieved by minimizing a least squares criterion $\mathcal{F}(x) = \lVert y - h \star x \rVert_2^2$ with added regularization:

- **Quadratic $(\ell_2)$ Regularization**: 
$$\mathcal{F}_\alpha(x) = \lVert y - h \star x \rVert_2^2 + \alpha \lVert d_1 \star x \rVert_2^2$$ 
where $d_1$ is a 2D Laplacian filter.

- **Convex Edge-Preserving Regularization**: 
$$\mathcal{F}_\alpha(x) = \lVert y - h \star x \rVert_2^2 + \alpha \sum _{p\sim q} \varphi(x_p - x_q)$$ 
where $\varphi(t) = \sqrt{t^2 + T^2} - T$.

A brief presentation of the main results can be found in the file _presentation_image_deconvolution.pdf_.

## Acknowledgements

This project was conducted in the context of the elective course "Numerical Image Processing" taught by Charles Soussen and Élisabeth Lahalle from L2S CentraleSupélec. The project team included Arthur Vogels.
