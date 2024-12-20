@doc raw"""

    wigner_eckart(j, j', m, m', k, q)

return the geometrical coefficient of the Wigner-Eckart Theorem, equal to `0` unless `q = m - m'`
```math
\left\langle\alpha jm\left|\mathbf{T}_{\mathrm{q}}^{(\mathbf{k})}\right| \alpha^{\prime} j' m^{\prime}\right\rangle=(-1)^{j-m}\left(\begin{array}{ccc}
j & k & j^{\prime} \\
-m & q & m^{\prime}
\end{array}\right)\left\langle\alpha j\left\|\mathbf{T}^{(\mathbf{k})}\right\| \alpha^{\prime} j^{\prime}\right\rangle
```
see *The theory of atomic structure and spectra* (11.15).
"""
wigner_eckart(j1, j2, m1, m2, k, q) = (-1)^Int(j1 - m1) * wigner3j(j1, k, j2, -m1, q, m2)


@doc raw"""

    reducedME(j, j')

reduced matrix element $\left\langle\alpha j\left\|\mathbf{J}^{(1)}\right\| \alpha^{\prime} j^{\prime}\right\rangle$
```math
\left\langle\alpha j\left\|\mathbf{J}^{(1)}\right\| \alpha^{\prime} j^{\prime}\right\rangle=\delta_{\alpha j, \alpha^{\prime} j^{\prime}}[j(j+1)(2 j+1)]^{1 / 2}
```
see *The theory of atomic structure and spectra* (11.20).
"""
reducedME(j1, j2) = j1 == j2 ? signedroot(j1 * (j1 + 1) * (2j1 + 1)) : RationalRoot(0)


@doc raw"""

    uncoup_T1(j1, j2, j, j1', j2', j', k)

Uncoupling formular with **T** operates only on $|\alpha_1j_1m_1\rangle$.
```math
\left\langle\alpha_{1} j_{1} \alpha_{2} j_{2} j\left\|\mathbf{T}^{(\mathbf{k})}\right\| \alpha_{1}^{\prime} j_{1}^{\prime} \alpha_{2}^{\prime} j_{2}^{\prime} j^{\prime}\right\rangle = \delta_{\alpha_{2} j_{2}, \alpha_{2}^{\prime} j_{2}^{\prime}}(-1)^{j_{1}+j_{2}+j^{\prime}+k}\left[j, j^{\prime}\right]^{1 / 2}\left\{\begin{array}{ccc}
j_{1} & j_{2} & j \\
j^{\prime} & k & j_{1}^{\prime}
\end{array}\right\}\left\langle\alpha_{1} j_{1}\left\|\mathbf{T}^{(\mathbf{k})}\right\| \alpha_{1}^{\prime} j_{1}^{\prime}\right\rangle
```
see *The theory of atomic structure and spectra* (11.38)
"""
function uncoup_T1(j11, j21, j1, j12, j22, j2, k)
    if j21 == j22 
        c = (-1)^Int(j11 + j21 + j2 + k) * signedroot((2j1 + 1) * (2j2 + 1)) * 
            wigner6j(j11, j21, j1, j2, k, j12)
    else
        c = RationalRoot(0)
    end
    return c
end


@doc raw"""

    uncoup_T2(j1, j2, j, j1', j2', j', k)

Uncoupling formular with **W** operates only on $|\alpha_2j_2m_2\rangle$.
```math
\left\langle\alpha_{1} j_{1} \alpha_{2} j_{2} j\left\|\mathbf{W}^{(\mathbf{k})}\right\| \alpha_{1}^{\prime} j_{1}^{\prime} \alpha_{2}^{\prime} j_{2}^{\prime} j^{\prime}\right\rangle = \delta_{\alpha_{1} j_{1}, \alpha_{1}^{\prime} j_{1}^{\prime}}(-1)^{j_{1}+j_2^{\prime}+j+k}\left[j, j^{\prime}\right]^{1 / 2}\left\{\begin{array}{ccc}
j_{1} & j_{2} & j \\
k & j^{\prime} & j_{2}^{\prime}
\end{array}\right\}\left\langle\alpha_{2} j_{2}\left\|\mathbf{W}^{(\mathbf{k})}\right\| \alpha_{2}^{\prime} j_{2}^{\prime}\right\rangle
```
see *The theory of atomic structure and spectra* (11.39)
"""
function uncoup_T2(j11, j21, j1, j12, j22, j2, k)
    if j11 == j12
        c = (-1)^Int(j11 + j22 + j1 + k) * signedroot((2j1 + 1) * (2j2 + 1)) *
            wigner6j(j11, j21, j1, k, j2, j22)
    else
        c = RationalRoot(0)
    end
    return c
end