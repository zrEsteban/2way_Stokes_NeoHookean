# Auditoría de convergencia, estabilidad y condicionamiento FSI

Fecha: 2026-08-22; actualización BDF2: 2026-08-23. Caso auditado principalmente:
`cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS`.

## Resultado ejecutivo

**Conclusión demostrada:** la implementación seleccionada **no corresponde exactamente** a un
“Robin–Neumann semi-implicit kinematically coupled β/AMP scheme with IQN-ILS correction” publicado.
Es un acoplamiento particionado **AMP-inspired pressure-Robin / Neumann**, iterado externamente,
con un participante estructural deal.II, una corrección IQN-ILS sólo sobre desplazamientos y una
inercia superficial Robin ajustada empíricamente. `semiImplicitBeta=1` sólo se valida; no multiplica
ningún término de las ecuaciones. Por ello la clasificación rigurosa es **AMP-inspired + IQN-ILS**,
no β-scheme ni AMP exacto.

El caso no es estable/convergente durante todo el intervalo configurado. Una corrida archivada de
10 pasos, hasta `1e-6 s`, termina; la corrida vigente avanza once pasos y en `1.2e-6 s` agota las
100 correcciones: el residuo de desplazamiento cumple `4.3876e-5 < 1e-4`, pero el cierre cinemático
falla con `max|u_f-v_s|=1.2540e-6 m/s` y razón física `4.138`. Newton de deal.II converge en dos
iteraciones y `minJ=1`; por tanto ese fallo es de la iteración de interfaz, no de Newton ni de la
definitud geométrica del sólido reducido.

No existe evidencia que permita declarar estabilidad incondicional. La mezcla temporal BDF2–Euler,
el valor parcial de masa virtual, el mapeo nearest-neighbour no adjunto y la ausencia de una
identidad energética discreta impiden aplicar directamente el teorema clásico del β-scheme.

Las etiquetas usadas abajo son **demostrada**, **inferida** y **no verificada**.

## Alcance e inventario

- **Demostrada:** no existe `AGENTS.md` bajo el árbol visible y el directorio no es un repositorio
  Git (`git status` devuelve “no es un repositorio”). No fue posible atribuir cambios ni distinguir
  modificaciones sin seguimiento; se conservaron resultados existentes.
- La biblioteca runtime se construye con los tres archivos del acoplador y se instala como
  `librobinRobinCoupling.so` ([Make/files](../src/robinRobinCoupling/Make/files)); el caso la carga en
  [controlDict](../cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/system/controlDict).
- El `wmake libso` posterior al cambio diagnóstico terminó correctamente. `wmkdepend` emitió avisos
  de headers opcionales y el enlazador avisó que no encontraba el `rpath` de PETSc, pero produjo la
  biblioteca. La ejecución requiere el `LD_LIBRARY_PATH` que ya establece `runCase.sh`.

## Flujo efectivo de llamadas

```text
controlDict carga librobinRobinCoupling.so
        |
fsiProperties: fluidSolidInterface robinRobin
        |
robinRobinCouplingInterface::evolve()
        |
        +-- updateDisplacement(): relajación o IQN-ILS de dΓ
        +-- moveFluidMesh()
        +-- fluid().evolve(): pimpleFluid, ddt=backward
        +-- updateForce()
              +-- presión Robin pdmsElasticWallPressure
              +-- exporta -t_f y u_f a CSV
              +-- ejecuta dealiiPdmsSolid
              +-- lee d_s, v_s, t_s, a_s
        +-- reportDealIIResiduals() y updateDealIIResidual()
        +-- repetir hasta las dos compuertas o nOuterCorr
        +-- sólo al converger: trial-state.bin -> accepted-state.bin
```

**Participante estructural efectivo — demostrada.** Con `useDealII yes`, `evolve()` llama al camino
deal.II y no ejecuta `solid().evolve()`
([robinRobinCouplingInterface.C:1269](../src/robinRobinCoupling/robinRobinCouplingInterface.C)).
El modelo solids4foam se construye y sirve como participante auxiliar para propiedades, patches y
transferencias, pero no gobierna la evolución constitutiva en este caso. El ejecutable externo y la
promoción del estado están configurados en
[fsiProperties:35](../cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/constant/fsiProperties).

## Ecuaciones realmente implementadas

### Fluido y condición Robin

El fluido usa `default backward`
([fluid/fvSchemes:19](../cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/system/fluid/fvSchemes)),
es decir BDF1 al arranque y esencialmente BDF2 con paso constante después. La condición de presión
construye

\[
  p^{n+1}+\frac{\rho_s h_s}{\rho_f}\,\partial_n p^{n+1}
  =p_\Gamma^{*}-\rho_s h_s\,a_{s,n}^{*}
\]

para presión dinámica, o la forma dividida por `rhoFluid` para presión cinemática. Los coeficientes
y el signo están en
[pdmsElasticWallPressureFvPatchScalarField.C:308](../src/robinRobinCoupling/pdmsElasticWallPressureFvPatchScalarField.C).
La implementación usa `coeff0=1`, `coeff1=rhoSolidHs/rhoFluid` y RHS
`prevPressure-rhoSolidHs*prevDdtUn`.

`rhoSolidHs()` obtiene `rho` y `impK` del sólido auxiliar, calcula
`a_p=sqrt(impK/rho)` y, si no hay espesor constante, toma `h_s=a_p Δt`
([pdmsElasticWallPressureFvPatchScalarField.C:33](../src/robinRobinCoupling/pdmsElasticWallPressureFvPatchScalarField.C)).
El caso fija `constantHs=3.469008593e-6 m`, por lo que la masa superficial es explícita y no cambia
automáticamente con `Δt`. **Inferida:** el valor equivale aproximadamente a media profundidad de onda
durante `1e-7 s`; no es el espesor físico completo.

### Sólido deal.II

La forma residual implementada es

\[
 \int_{\Omega_s} P(F):\nabla\delta d+\rho_s\,a_s^{n+1}\cdot\delta d\,d\Omega
 -\int_{\Gamma}t_s\cdot\delta d\,d\Gamma=0,
\]

con neo-Hooke compresible y

\[
 v_s^{n+1}=\frac{d_s^{n+1}-d_s^n}{\Delta t},\qquad
 a_s^{n+1}=\frac{d_s^{n+1}-d_s^n-\Delta t v_s^n}{\Delta t^2}.
\]

La versión auditada originalmente empleaba esas expresiones BE. Desde la actualización del
2026-08-23, el caso seleccionado declara `time integration = bdf2` y, después del arranque BE, usa

\[
v_s^{n+1}=\frac{3d_s^{n+1}-4d_s^n+d_s^{n-1}}{2\Delta t},\qquad
a_s^{n+1}=\frac{3v_s^{n+1}-4v_s^n+v_s^{n-1}}{2\Delta t}.
\]

Por consiguiente el Jacobiano inercial pasa de `M_s/deltaT^2` durante el arranque a
`9 M_s/(4 deltaT^2)` en BDF2. Esto está demostrado por
[dealiiPdmsSolid.cc:195](../src/dealiiPdmsSolid/dealiiPdmsSolid.cc) y
[dealiiPdmsSolid.cc:240](../src/dealiiPdmsSolid/dealiiPdmsSolid.cc). El Jacobiano de Newton contiene
el coeficiente inercial seleccionado más la tangente material
([dealiiPdmsSolid.cc:322](../src/dealiiPdmsSolid/dealiiPdmsSolid.cc)).

Los estados nuevos almacenan `(d^n,v^n,d^{n-1},v^{n-1},historyDepth)`. Los estados legados de dos
vectores continúan siendo legibles y realizan un paso BE antes de habilitar BDF2.

En la frontera deal.II se impone

\[
t_s=t_{in}+Z(u_{in}-v_s).
\]

El caso usa `Z=0`; por tanto es **Neumann puro**, `t_in=-t_f`
([dealiiPdmsSolid.cc:265](../src/dealiiPdmsSolid/dealiiPdmsSolid.cc),
[solid.prm:13](../cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/dealiiSolid/solid.prm)).
No existe compensación estructural explícita de la masa Robin del fluido. **Inferida:** esa masa actúa
como precondicionador/splitting unilateral; no está duplicada en el sólido, pero tampoco se demostró
la cancelación algebraica exigida por un AMP clásico.

### Significado real de β

**Demostrada:** `semiImplicitBeta_` se lee y sólo se comprueba que valga uno
([robinRobinCouplingInterface.C:45](../src/robinRobinCoupling/robinRobinCouplingInterface.C),
[robinRobinCouplingInterface.C:186](../src/robinRobinCoupling/robinRobinCouplingInterface.C)). No hay
otro uso. No existe predictor β, partición de presión por β ni paso estructural cinemáticamente
acoplado reconocible. `predictSolid no` además desactiva el predictor opcional. El nombre de la
configuración no constituye evidencia algorítmica.

## Consistencia DAE y BDF2–Euler

Las restricciones DAE se hacen cumplir de modo iterativo, no monolítico:

- `D u_f=0`: PIMPLE/presión del fluido;
- `dot(d_s)=v_s`: BDF2 dentro de deal.II después del arranque BE;
- `u_f-v_s=0`: compuerta física de interfaz;
- `t_f+t_s=0`: entrada Neumann, pero no comprobación independiente cuando `Z=0`;
- equilibrio constitutivo/inercial estructural: Newton deal.II.

En la versión inicialmente auditada, la aceleración que volvía a la Robin era BE mientras el fluido
usaba `backward`. Para
paso constante, el defecto es

\[
D_{BDF2}v^{n+1}-D_{BE}v^{n+1}
=\frac{v^{n+1}-2v^n+v^{n-1}}{2\Delta t}=O(\Delta t).
\]

**Demostrada para la versión anterior:** no había término de código que lo cancelara. **Corregida en
el caso seleccionado:** después del arranque, la velocidad y aceleración estructurales usan BDF2 y
son las que se escriben en `robin-out.csv`; el acoplador devuelve esa aceleración a la presión Robin.
El defecto BE–BDF2 anterior desaparece para paso constante. Continúa pendiente derivar una identidad
energética para el splitting completo. En la versión anterior era un defecto de splitting de primer
orden que alteraba la identidad energética. No se
verificó experimentalmente el orden porque no existe aún una terna convergida común
`Δt,Δt/2,Δt/4`. El primer paso combina el arranque `backward` de OpenFOAM con estado inicial deal.II
BE. Paso variable no es compatible: OpenFOAM podría variarlo, pero deal.II lee un `delta t` fijo; el
caso lo evita con `adjustTimeStep no`
([controlDict:20](../cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/system/controlDict)).

## IQN-ILS y convergencia de interfaz

La variable IQN es el desplazamiento de puntos de la interfaz fluida. El residuo acelerado es

\[
r_d^k=d_{s,\Gamma}^k-d_{f,\Gamma}^k.
\]

No se acelera el vector apilado `[u_f-v_s,t_f+t_s]` en la configuración `iqnils`; ese camino sólo
existe como opción diferente `iqnilsCombined`. Las secantes son
`V_i=Δr_d`, `W_i=Δd_s`, y la corrección aproxima la inversa del Jacobiano del punto fijo
([robinRobinCouplingInterface.C:423](../src/robinRobinCoupling/robinRobinCouplingInterface.C)).

La implementación:

- reinicia historia por paso porque `iqnReuseAcrossTimeSteps no`;
- retiene hasta 12 modos;
- filtra dependencia por Gram–Schmidt con umbral relativo `1e-8`;
- resuelve una QR modificada, con corte adicional `1e-12 max(diag(R))`;
- limita la corrección a `0.25 max|r_d|`;
- usa relajación fija 0.25 durante el arranque.

**Defecto demostrado:** no escala componentes/zonas ni registra `ΔX,ΔR`, valores singulares,
condición reducida o una estimación de `rho(I-BJ)`. Las columnas tienen las mismas unidades al ser
sólo desplazamientos, pero el producto Euclídeo no pondera área ni masa. La prueba rápida nueva
confirma que una historia colineal de rango 2/3 se detecta en el análogo reducido; esto prueba el
test, no una SVD runtime del caso.

La condición local necesaria sigue siendo

\[
\rho(I-B_kJ_\Gamma)<1.
\]

**No verificada numéricamente** para la corrida porque el código no persistía las secantes. La
oscilación y posterior estancamiento de la razón física por encima de uno es evidencia empírica de
que el operador corregido no es contractivo en ese paso, no una medición de su radio espectral.

### Residuales y aceptación

El código calcula

\[
r_u=u_f-v_s,\quad r_t=t_f+t_s,\quad r_R=r_t-Zr_u,
\]

y potencia puntual `t_f·u_f+t_s·v_s`
([robinRobinCouplingInterface.C:1020](../src/robinRobinCoupling/robinRobinCouplingInterface.C)). Los
relativos usan máximos de magnitud y pisos configurables. Los pisos `3.0303e-3 m/s` y `100 Pa`
eliminan el cociente cercano a uno en reposo. La prueba nueva obtiene `rel_u=5.72e-12` y
`rel_t=1.73e-14` para perturbaciones de redondeo.

El paso sólo se acepta si `r_d <= outerCorrTolerance` **y** todas las razones físicas son menores o
iguales a uno. Al agotar 100 iteraciones aborta; el estado estructural trial no se promueve
([robinRobinCouplingInterface.C:1298](../src/robinRobinCoupling/robinRobinCouplingInterface.C)).

**Limitación importante:** con `Z=0`, el acoplador reconstruye
`t_s=-t_f` directamente
([robinRobinCouplingInterface.C:385](../src/robinRobinCoupling/robinRobinCouplingInterface.C)). Por
eso `max|r_t|=0` es tautológico y no valida independientemente la tracción calculada por el FEM.

## Rollback e historias

deal.II lee `accepted-state.bin` en cada invocación, escribe `trial-state.bin`, y la interfaz copia
trial a accepted sólo tras ambas compuertas
([dealiiPdmsSolid.cc:170](../src/dealiiPdmsSolid/dealiiPdmsSolid.cc),
[robinRobinCouplingInterface.C:1317](../src/robinRobinCoupling/robinRobinCouplingInterface.C)). Esto
demuestra rollback estructural entre iteraciones externas. La prueba de hash de sólo lectura verificó
que inspeccionar el trial no cambia accepted (`d606aee0d88bec88…`).

**Parcial:** los campos fluidos no vuelven al estado de comienzo del paso en cada corrección; PIMPLE
continúa desde el iterado anterior, que es la práctica de fixed-point habitual. La malla se mueve por
incrementos respecto al iterado previo. No hay checkpoint transaccional completo de fluido+malla si
el proceso aborta.

## Potencia, energía y mapeo

La identidad utilizada es

\[
P_\Gamma=\int_\Gamma[t_f\cdot(u_f-v_s)+(t_f+t_s)\cdot v_s]d\Gamma.
\]

Se añadió al camino deal.II el registro de fuerza residual integrada, potencia integrada, escala de
potencia y defecto relativo integrado, y se incorporó este último a la compuerta. La biblioteca
compila. El camino interno ya tenía una integral equivalente. Este diagnóstico se encuentra en
[robinRobinCouplingInterface.C:1078](../src/robinRobinCoupling/robinRobinCouplingInterface.C).

Un smoke aislado posterior a recompilar cuantificó `integrated traction defect=(0 0 0) N`,
`integrated power defect=8.825653233e-18 W` y defecto relativo integrado `1.060868361e-11` en el
primer paso. Una primera versión normalizada sólo por el trabajo instantáneo produjo una razón
espuria cerca del reposo; se corrigió usando el piso físico `AΓ * tractionScale * velocityScale` y se
repitió el smoke con PASS. Aun así, no están implementadas las energías globales cinética del fluido, cinética y elástica
del sólido, disipación viscosa y trabajo externo. Por tanto no existe cierre del balance energético
total; sólo el defecto de potencia interfacial puntual en logs antiguos.

El caso reducido transfiere tracción a cuadratura deal.II por vecino más cercano y reconstruye
velocidad de cara promediando valores en puntos
([robinRobinCouplingInterface.C:283](../src/robinRobinCoupling/robinRobinCouplingInterface.C)). No usa
un par de operadores primal/dual que satisfaga `T_t = M_s^{-1} T_u^T M_f`. **FAIL demostrado por
estructura del algoritmo:** no se puede garantizar conservación de trabajo discreta. Conservar fuerza
global sería insuficiente. En el caso grande, AMI transfiere campos en ambas direcciones, pero tampoco
se encontró una prueba de adjunción ponderada.

## Condicionamiento y geometría

### Participante reducido deal.II

La malla tiene 3000 celdas y 12948 grados de libertad según los logs. El Newton vigente converge de
`|R|≈1.42e-8` a `3.75e-16` en una actualización, con 24 iteraciones CG, búsqueda lineal `alpha=1` y
`minJ=1`. Hay guardia `J > 1e-8` y búsqueda lineal de hasta diez reducciones
([dealiiPdmsSolid.cc:298](../src/dealiiPdmsSolid/dealiiPdmsSolid.cc)). No se registra autovalor mínimo
ni estimación de condición de `J_s`; la definitud es **inferida**, no demostrada, a partir de que CG
converge rápidamente en este régimen casi sin deformación.

### Caso grande solids4foam

Las cifras conocidas están reproducidas:

- 496200 celdas sólidas;
- 27426 celdas con determinante de calidad menor que 0.001;
- determinante mínimo `4.034028669e-6`;
- aspecto máximo `220.9138358`;
- interfaz fluida 94798 caras, sólida 50250 caras.

Evidencia: [log.checkMesh.solid](../cases/pdmsMicrochannelFSI/singlePhase/youngPDMS/log.checkMesh.solid),
[log.checkMesh.fluid](../cases/pdmsMicrochannelFSI/singlePhase/youngPDMS/log.checkMesh.fluid). Son
**amplificadores de condicionamiento y sensibilidad local**, no prueba por sí solos de inestabilidad.

El fallo histórico pedido también se reproduce documentalmente: a `t=0.00012`, iteración Robin 7,
el desplazamiento acumulado informado es `2.635662528`; después el neo-Hookeano aborta con
`J_min=-0.3152353905`
([log.solids4Foam.restart:448](../cases/pdmsMicrochannelFSI/runHistory/youngPDMS_before_clean_20260820_0854/log.solids4Foam.restart)).
La cadena “corrección grande/desfasada → inversión local → aborto” es **fuertemente inferida** por la
secuencia temporal. No está demostrado que IQN sea la única causa; Courant máximo 2.54, presión
iterativa costosa y la malla mala también amplifican el evento.

`interfaceDeformationLimit 0` no limita el incremento: la comparación de solids4foam es
`maxDelta < limit`, de modo que con cero siempre se elige movimiento de malla completo
([fluidSolidInterface.C:857](/home/ezamora/Workspace/solids4foam/src/solids4FoamModels/fluidSolidInterfaces/fluidSolidInterface/fluidSolidInterface.C)).
La relajación sólida 0.02 y hasta 20000 correctores pertenecen al sólido auxiliar en el caso deal.II;
no controlan su Newton. Sí son relevantes para el caso grande gobernado por solids4foam.

## Resultados reproducibles

### Pruebas rápidas ejecutadas

```bash
python3 cases/pdmsMicrochannelFSI/tests/fsi_audit_quick.py \
  --log cases/pdmsMicrochannelFSI/stabilityStudy/semiImplicitBetaIQNILS/log.run
openfoam2512 -c 'cd src/robinRobinCoupling && wmake libso'
```

Resultados: normalización cerca de reposo PASS; historia deficiente PASS; invariante de estado
accepted PASS; compilación PASS. Un smoke de un paso en `/tmp/fsi-audit-smoke.QKiVmm` terminó con
`End`, razón física 0.003791 y defecto integrado relativo `1.061e-11`; no modificó el caso original.
El analizador encuentra 12 tiempos, último `1.2e-6`, fatal sí,
razón física final 4.138 y residuo final `4.3876e-5`.

La actualización BDF2 fue validada además con `src/dealiiPdmsSolid/tests/run-smoke.sh`: el primer
paso informó `backwardEuler startup`, el segundo `BDF2`, y las pruebas material, restart y límites
Robin pasaron. Un smoke FSI limpio de dos pasos en `/tmp/fsi-bdf2-smoke.YKl81V` terminó con `End`:
razones físicas 0.003791 y 0.020178; en el segundo paso deal.II informó BDF2, Newton convergió en dos
iteraciones, el defecto de potencia integrado fue `1.196e-16 W` y el relativo `1.437e-10`.

### Corridas existentes

| Corrida | Alcance | Resultado |
|---|---:|---|
| `reports/log.semiImplicitBeta.gate10.mass0p5x` | 10 pasos, `dt=1e-7` | PASS hasta `1e-6`; razón final 0.561 |
| `stabilityStudy/semiImplicitBetaIQNILS/log.run` | objetivo `1e-4` | FAIL en `1.2e-6`, iteración 100; razón 4.138 |
| histórico `youngPDMS.../log.solids4Foam.restart` | caso grande | FAIL con `J_min=-0.3152` a `1.2e-4` |

El gate corto costó 177 s de reloj. Un solo punto de 12 pasos llegó a unas 100 correcciones en el
último paso; los barridos de tolerancia, tiempo, densidad, `Z_s`, β y relajación completos requieren
horas. No se lanzaron de forma inadvertida sobre el caso productivo.

### Matriz de pruebas requerida

| Prueba | Estado | Evidencia/razón |
|---|---|---|
| Residuales cerca del reposo | PASS | prueba algebraica con pisos configurados |
| Equilibrio de tracciones | PARTIAL | signo y `Z=0` demostrados; residual runtime no independiente |
| Conservación de trabajo no conforme | FAIL | nearest-neighbour no adjunto; sin identidad matricial |
| Rollback IQN | PARTIAL | estado deal.II PASS; fluido/malla no transaccionales |
| Historia IQN de rango deficiente | PASS | prueba reducida; filtro source presente |
| Tolerancias `1e-4,1e-6,1e-8` | NOT VERIFIED | no hay tres corridas comparables |
| Refinamiento `dt,dt/2,dt/4` | NOT VERIFIED | requiere sincronizar `controlDict`, `solid.prm` y `constantHs` |
| Orden temporal observado | NOT VERIFIED | no hay terna convergida; esperado ≈1 por defecto temporal |
| Barrido `rho_s/rho_f` | NOT VERIFIED | no ejecutado por coste |
| Barrido `Z_s`, β, relajación | PARTIAL | logs de masa virtual existentes; β distinto de 1 es rechazado |
| Energía sin trabajo externo | NOT VERIFIED | smoke de potencia pasa; faltan energías de volumen y caso sin forzamiento |
| Inversión de J | PASS | fallo archivado reproducido y localizado |

Para analizar secantes sin matrices densas, el script acepta un NPZ con matrices delgadas `dX,dR`:

```bash
python3 cases/pdmsMicrochannelFSI/tests/fsi_audit_quick.py --secant histories.npz
```

Calcula SVD sólo en el espacio multisecante y proyecta `Jhat=ΔR(ΔX)^+`; no forma el operador de
interfaz completo. Falta instrumentar la exportación runtime de esas columnas antes de poder estimar
`rho(I-B Jhat)`.

## Defectos y cambios

1. **FAIL:** la etiqueta β no controla las ecuaciones. No se corrigió silenciosamente porque hacerlo
   reemplazaría el esquema, fuera del alcance permitido.
2. **CORREGIDO para el caso seleccionado:** BDF2 consistente en las dos ecuaciones estructurales,
   arranque BE y estado histórico compatible. Falta el estudio de orden y energía global.
3. **FAIL:** tracción residual deal.II con `Z=0` es construido como cero.
4. **FAIL:** mapeo externo nearest-neighbour no es adjunto/conservativo en trabajo.
5. **FAIL:** la corrida vigente pierde convergencia cinemática a `1.2e-6`.
6. **PARTIAL:** normalización cerca del reposo fue corregida previamente mediante pisos y ahora pasa.
7. **PASS:** deal.II tiene guardia `J`, Newton consistente y búsqueda lineal.
8. **Cambio realizado:** diagnóstico de fuerza y potencia integradas en el camino deal.II, incluido en
   la compuerta física; biblioteca recompilada.
9. **Cambio realizado:** prueba rápida no destructiva y analizador reducido de secantes/logs.

## Riesgos y próximos pasos priorizados

1. Exportar las columnas IQN delgadas y registrar singular values, rango, `cond(V)` y una estimación
   reducida de contracción por iteración. Reiniciar historia ante crecimiento persistente del residual.
2. Sustituir nearest-neighbour por un par de transferencias primal/dual adjuntas y probar
   `u_s^T M_s t_s + u_f^T M_f t_f` con campos aleatorios y constantes.
3. Hacer independiente la tracción sólida: integrar el `P N` FEM real en la interfaz en vez de
   reconstruirla desde la entrada.
4. Definir y derivar formalmente una única discretización temporal. Si se conserva BDF2 del fluido,
   añadir el término de consistencia o llevar el sólido/interfaz al mismo operador; después repetir
   el estudio de orden.
5. Añadir energías de volumen y trabajo externo antes de afirmar estabilidad. Ejecutar primero un
   caso sin forzamiento de 10 pasos.
6. Ejecutar la matriz larga en copias aisladas, con un presupuesto por corrida y gate temprano:
   tolerancias; tres `dt`; densidad; `rho_s h_s`; relajación/cap IQN. β no puede barrerse hasta que sea
   una variable matemática real.
7. Para el caso grande, activar límite de deformación, aplicar trust region/búsqueda lineal al
   incremento exterior y reparar/refinar las 27426 celdas deficientes antes de repetir la inversión.

## Tabla final

| Condición | Evidencia | Estado | Riesgo | Acción |
|---|---|---|---|---|
| Participante efectivo identificado | rama `useDealII`, no `solid().evolve()` | PASS | bajo | mantener prueba runtime |
| Esquema β exacto | β sólo validado, no usado | FAIL | alto | derivar/implementar por separado |
| AMP exacto | Robin unilateral, BE/BDF2 y masa parcial | FAIL | alto | clasificar AMP-inspired |
| Cierre cinemático | falla en `1.2e-6`, razón 4.138 | FAIL | alto | diagnosticar operador IQN |
| Cierre dinámico independiente | tracción reconstruida con `Z=0` | FAIL | alto | integrar tracción FEM |
| Newton deal.II, régimen reducido | 2 iteraciones, CG 24, `minJ=1` | PASS | bajo localmente | monitorear condición |
| Convergencia IQN local | gate 10 pasos pasa; paso 12 falla | PARTIAL | alto | secantes/SVD y restart adaptativo |
| Rollback estructural | accepted/trial y promoción tardía | PASS | medio | añadir prueba ejecutable aislada |
| Rollback fluido+malla | sin checkpoint completo | PARTIAL | medio | documentar/restart seguro |
| Normalización en reposo | pisos + prueba rápida | PASS | bajo | conservar escalas físicas |
| Trabajo no conforme | operador nearest no adjunto | FAIL | alto | transferencia dual conservativa |
| Potencia de interfaz integrada | código añadido, compila y smoke PASS | PASS | medio | extender a corrida larga |
| Balance energético total | energías de volumen ausentes | NOT VERIFIED | alto | instrumentar y cerrar balance |
| Consistencia temporal local fluido–sólido | BDF2/BDF2 tras arranque; smoke 2 pasos | PASS | medio | medir orden y energía global |
| Orden temporal ≈1 | inferencia, sin terna | NOT VERIFIED | medio | estudio `dt/2/dt/4` |
| Masa añadida/`Z_s` óptimos | un gate corto con media masa | PARTIAL | alto | barrido controlado |
| Protección geométrica externa | cap IQN 0.25; límite malla 0 | PARTIAL | alto | límite/trust region y gate J |
| Inversión histórica | `J=-0.3152`, paso/iteración exactos | PASS | alto | malla + line search exterior |
| Estabilidad incondicional | hipótesis del teorema no se cumplen | FAIL | crítico | no declarar; demostrar tras rediseño |
