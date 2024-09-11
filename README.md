# L’orbite du télescope spatial James Webb

## Mise en équation pour la simulation numérique

La position $\vec{r_J}$ et la vitesse $\vec{v_J}$ de JWST sont définies dans le référentiel tournant à la vitesse $\Omega_T$, dont l'origine est le point de Lagrange $L_2$.

Le point $S$ désigne le Soleil. Le point $T$ désigne le barycentre Terre-Lune. L'axe $X$ est dans la direction Soleil-Terre.

$$
\vec{r_J}=
\left(\begin{matrix}
X\\
Y\\
Z
\end{matrix}\right)\quad
\vec{v_J}=
\left(\begin{matrix}
V_X\\
V_Y\\
V_Z
\end{matrix}\right)
$$

L'équation dynamique est la suivante ($\mu_S$ est le paramètre gravitationnel standard du Soleil, $\mu_T$ est celui du couple Terre-Lune) :

\begin{align}
\frac{d\vec{r_J}}{dt}&=\vec{v_J}\\
\frac{d\vec{v_J}}{dt}&=-\mu_S\frac{\vec{SL_2}+\vec{r_J}}{\left\lVert\vec{SL_2}+\vec{r_J}\right\lVert^3}-\mu_T\frac{\vec{TL_2}+\vec{r_J}}{\left\lVert\vec{TL_2}+\vec{r_J}\right\lVert^3}+\Omega_T^2\left(\begin{matrix}\left\lVert\vec{SL_2}\right\lVert+X\\Y\\0\end{matrix}\right)+2\Omega_T\left(\begin{matrix}V_Y\\-V_X\\0\end{matrix}\right)
\end{align}
