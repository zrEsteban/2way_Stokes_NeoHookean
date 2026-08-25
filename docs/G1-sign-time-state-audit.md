# G1 — signos, integración temporal, estados y persistencia Robin

Estado: **FAIL**. Commit de entrada: `9897ceba2de926397c261f1f8c9a7d22a5b739e8`.
G0 permanece aprobado; no se inició G2 ni se modificaron tolerancias, datos
físicos, transferencia, IQN-ILS o formulación FSI.

## Convención de signos

La convención verificada es

```text
t_f = sigma_f n_f = tau_f n_f - p n_f
f_s = -T_{f->s} t_f
```

`robinRobinCouplingInterface::updateForce()` forma
`faceZoneViscousForce - faceZonePressureForce*n_f`; el adaptador deal.II
escribe su negativo como carga estructural. El residual débil del sólido suma
esa carga externa. En una interfaz plana superior, `n_f=+e_z` y `p>0` producen
`t_f=-p e_z` y `f_s=+p e_z`, compresión hacia el PDMS exterior. La prueba
unitaria y el smoke real obtuvieron desplazamiento estructural `+z`.

La forma reservada para gates futuros queda documentada, no implementada:

```text
A_f = df_s/dd_Gamma
R_stab(d) = R_s(d) - f_s^k - A_f(d-d^k)
J_stab = J_s - A_f
```

## BE, BDF2 y estado aceptado

El código usa, mediante `TimeIntegration.H`:

```text
BE:   v=(d-dn)/dt; a=(v-vn)/dt; da/dd=1/dt^2
BDF2: v=(3d-4dn+dnm1)/(2dt)
      a=(3v-4vn+vnm1)/(2dt); da/dd=9/(4dt^2)
```

La velocidad y aceleración exportadas usan las mismas funciones que residual
y Jacobiano. Las pruebas polinomiales y de derivada por diferencias finitas
pasan. `accepted-state.bin` permanece inmutable durante correctores; cada solve
escribe `trial-state.bin`, y la promoción ocurre una vez después de todos los
gates. `stateAudit` (default `false`) contó 10 aceptaciones en continuo y 5+5
en partido. El primer paso retomado informó BDF2. Un fallo anterior a la
promoción no modifica el aceptado.

## Persistencia Robin

Clasificación:

- `constantHs`: parámetro físico persistente.
- `value`, `coeff0`, `coeff1`, `rhs`, `prevPressure`, `prevAcceleration`:
  estado/coeficientes aceptados persistentes.
- `rhoSolidHs`: cache derivado que debe persistir porque reconstruirlo sobre
  una geometría movida no es restart-idempotente.
- datos transferidos, residuales e historias IQN: iterativos; no se escriben.

Se corrigieron: omisión de `constantHs`; lectura de `prevPressure` que antes
sobrescribía `value`; ausencia de `prevAcceleration`; pérdida de coeficientes
Robin al leer; y pérdida de estados en constructores de copia, mapping,
`autoMap` y `rmap`. El cache usa semántica por valor. La condición actualiza
presión/aceleración finales antes de promover el estado deal.II.

Los diccionarios iniciales sin `coeff0/coeff1/rhs` siguen aceptados con los
defaults históricos; los restart que contienen las tres entradas las
restauran. `decomposePar -no-libs` permanece obligatorio como defensa.

## Pruebas y resultados

| Prueba | Resultado |
|---|---|
| build `librobinRobinCoupling.so` | PASS |
| CTest material/tiempo/signo (3/3) | PASS |
| smoke deal.II: BE, BDF2, signo, topología y límites Robin | PASS |
| round-trip e idempotencia Robin | PASS |
| `decomposePar` con plugin y `-no-libs` | PASS |
| validador negativo G0 | PASS |
| mínimo serial continuo, 10 pasos, `End`, sin NaN/Inf | PASS |
| aceptación exactamente una vez | PASS (10; partido 5+5) |
| restart serial continuo/partido, `rtol=5e-5` | **FAIL** |
| mínimo/restart MPI de G1 | NOT RUN: bloqueado por el fallo serial previo |

Comandos principales:

```bash
wmake libso src/robinRobinCoupling
src/dealiiPdmsSolid/tests/run-smoke.sh
scripts/test-robin-serialization.sh
scripts/test-decomposition-preserves-constantHs.sh
G1_RESTART_MODES=serial G1_RESTART_WORK=/tmp/g1-restart-20260825-h \
  scripts/test-g1-restart.sh
```

El continuo serial alcanzó `t=1e-6`, 10 pasos y `End`. Métrica final:
residuo geométrico relativo `7.8076143124498678e-05`, ratio físico
`0.50885219891990696`, `max|r_u|=1.541974818387503e-07`, sin NaN/Inf.
SHA-256: `fsiResiduals.dat`
`8bc7ad2dd28176344230a6bef559688b6784aab19935a2b415effeea0aa080e0`;
`robin-out.csv`
`cceaabe1995fc8bf45b148a7375f0b632e5afb2c4baf6312c3ad1243c0a14cc1`.

El partido también alcanzó `End` y usó BDF2, pero su `robin-out.csv` fue
`cfd677ac4f397391696e22b0492c1f2cc5a2ce6d77241e41d29ece8980b240dd`.
Los máximos relativos continuo/partido fueron aproximadamente 4.0–6.3 % en
desplazamiento, 0.8–7.9 % en velocidad, 1.4–4.8 % en tracción y 11.9–42.3 %
en aceleración: exceden `5e-5`. La primera divergencia significativa aparece
en el primer solve de presión a `t=6e-7`: residual inicial
`0.059691256396702426` continuo frente a `0.059624512273709143` restart,
aunque el solve de velocidad anterior coincide a los dígitos impresos.

## Decisión y riesgos

G1 no puede aprobarse: falta identificar y persistir el estado del participante
fluido que causa la divergencia de presión posterior al restart. No se ajustó
ninguna tolerancia para ocultarlo. MPI restart no debe considerarse validado.

El siguiente trabajo recomendado sigue siendo **G1**, no G2: instrumentar de
forma opt-in los términos de la ecuación de presión y sus campos `oldTime`,
localizar el primer operando distinto en `t=6e-7`, persistirlo y repetir la
matriz serial/MPI. No iniciar G2 hasta cerrar esta deuda.
