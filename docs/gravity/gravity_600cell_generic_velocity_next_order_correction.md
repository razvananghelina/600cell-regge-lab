# Correction protocol: next-order direct-limit timeout

Date: 2026-08-21

Registered implementation commit: `bd68c52`.

Preserved first-execution record:

```text
reproducible/gravity_600cell_generic_velocity_next_order_first_timeout.json
SHA-256 0aa9d9fca5d6f0dc4dd5f6a5ac93e04ceb5885ab277cede2d4a4d75df919e1db
```

The first execution passed only the frozen-provenance gate and then spent
more than ten minutes inside

```text
sympy.limit(2*F/h,h,0+)
```

after the complete exponential path had been substituted into the full
unreduced derivative.  It was interrupted with exit code `130`.  No
next-order coefficient, common root or physical outcome had been evaluated.

## Exact mechanical replacement

Introduce, before taking a limit,

```text
tau=sqrt(rho),
q=(L_plus-L_minus)/tau,
w=600*sqrt(3)*(L_minus^2-L_plus^2)/tau.
```

The action is identically

```text
S=tau[A(L_minus,L_plus,q,M)+w*eta(q)],

A=360(L_minus+L_plus)*sqrt(1+q^2/4)*epsilon(q)-8*pi*M.
```

At fixed endpoints, exact differentiation gives

```text
2F/tau=A-q*A_q-q*w*eta_q.
```

The exact boundary derivatives are

```text
P_minus=(L_minus/2)
 [tau*A_Lminus-A_q+1200*sqrt(3)*L_minus*eta-w*eta_q],

P_plus=(L_plus/2)
 [tau*A_Lplus+A_q-1200*sqrt(3)*L_plus*eta+w*eta_q].
```

For a path

```text
L_minus=exp(xm1*h+xm2*h^2),
L_plus =exp(xp1*h+xp2*h^2),
tau=c*h,
```

freeze the exact first jets

```text
q0=(xp1-xm1)/c,
q1=[xp2+xp1^2/2-xm2-xm1^2/2]/c,

w0=-1200*sqrt(3)*q0,
w1=(1200*sqrt(3)/c)
   [xm2-xp2+xm1^2-xp1^2].
```

Extract every registered first coefficient by the multivariate chain rule at
`h=0`.  This is algebraically identical to the direct limit but prevents
SymPy from rediscovering the scaling through the complete transcendental
expression.

## Scope of the correction

Do not change the action, state, endpoint ansatz, coefficient definitions,
root census, sampled controls, thresholds, composition equations or outcome
hierarchy.  The existing direct 100-decimal controls against the original
unexpanded `F` and momenta remain mandatory and are the independent check of
the scaled chain-rule extraction.
