# Image-Deconvolution

L'objectif est d'estimer une image nette $\widehat{x}[m,n]$ à partir d'une image floue $y[m,n]$ obtenue par convolution de l'image nette $x[m,n]$ avec la réponse impulsionnelle $h[p,q]$, suivie d'une corruption par bruit additif $\epsilon[m,n]$. Pour cela, on cherche à minimiser un critère des moindres carrés $\mathcal{F}(x) = \lVert y - h \star x \rVert_2^2$, auquel on ajoute une régularisation :

- Par régularisation quadratique ($\ell_2$) : 
$$\mathcal{F}_\alpha(x) = \lVert y - h \star x \rVert_2^2 + \alpha \lVert d_1 \star x \rVert_2^2$$ 
où $d_1$ est un filtre laplacien 2D.

- Par régularisation convexe préservant les contours : 
$$\mathcal{F}_\alpha(x) = \lVert y - h \star x \rVert_2^2 + \alpha \sum _{p\sim q} \varphi(x_p - x_q)$$ 
où $\varphi(t) = \sqrt{t^2 + T^2} - T$.

Une rapide présentation des principaux résultats obtenus se trouve dans le fichier _presentation_image_deconvolution.pdf_.
