# G2-B.2 — conexión dual viva FVM–FEM

Fecha: 2026-08-26

Commit de entrada: `f88684c877ee91a88f14356a1c7ee939c6652a66`

Estado histórico: **BLOCKED en preflight; bloqueo corregido posteriormente por
G2-B.1-RUNTIME. G2-B.2 continúa pendiente y no fue reintentado.**

## Resultado del preflight obligatorio

La API de G2-B.1 existe como clase C++ reutilizable en
`src/dealiiPdmsSolid/GeneralizedInterfaceLoad.H`, pero todavía no forma parte de
la ruta productiva del ejecutable estructural:

```text
rg GeneralizedLoad src/dealiiPdmsSolid
  tests/test_interface_sparsity.cc: ... GeneralizedLoad ...
  GeneralizedInterfaceLoad.H: class GeneralizedLoad
```

No hay ninguna inclusión ni llamada a `GeneralizedLoad::receive()` o
`add_to_newton()` desde `dealiiPdmsSolid.cc`. En particular:

- `PdmsSolid::assemble_newton()` siempre ensambla la carga Robin legacy con
  `FEFaceValues`, `nearest()`, `JxW` y muestras en Pa;
- el vector `rhs` y la matriz productivos nunca reciben `f_Gamma` en N ni
  `J_Gamma` en N/m;
- la matriz ampliada puede prepararse en `setup()`, pero no existe un canal
  productivo que entregue los mensajes sellados;
- `PdmsSolid::run()` aborta explícitamente si `interface transfer` no es
  `legacyNearestNeighbour`.

Por tanto, los siguientes elementos aprobados en G2-B.1 están demostrados sólo
en el arnés deal.II, no en la API/ruta estructural viva:

| requisito de preflight | estado productivo |
|---|---|
| fuerza generalizada directa en N | **ausente** |
| ruta sin `FEFaceValues/JxW` | **ausente** |
| tangente dispersa en Newton productivo | **ausente** |
| protocolo `timeIndex/outerCorrector/hashes` | clase probada, **sin transporte productivo** |
| matriz replicada por rank | arnés probado, **sin recepción productiva** |
| diferencia finita | arnés probado, **no sobre la llamada productiva** |

Esto satisface exactamente la condición de detención indicada en el apartado
0 del gate: «si alguno existe solamente en el arnés y no en la API C++
estructural real, detenerse y declarar G2-B.2 bloqueado».

## Consecuencias

No se modificó el acoplador, no se intentó el handshake, no se activó
`dualConservative`, no se ejecutaron las etapas A–I y no se creó un caso dual.
Tampoco se ejecutó el baseline fuerte. Hacerlo habría obligado a introducir un
canal temporal o a caer en la cuadratura legacy, ambos prohibidos.

G2-B.2, G2-B y G2 global permanecen **no aprobados**. G2-B.1 requiere una
remediación acotada antes de reabrir este gate.

## Cambio mínimo requerido

Crear un subgate **G2-B.1-RUNTIME** que promueva el contrato ya validado a la
ruta estructural real, sin conectar aún OpenFOAM:

1. añadir a `PdmsSolid` un estado de protocolo propietario del corrector;
2. exponer un transporte versionado/no temporal para recibir el manifiesto,
   `ForceMessage` y `TangentMessage` completos;
3. separar `assemble_newton_legacy()` de la rama generalizada;
4. en la rama generalizada llamar realmente a `receive()` y
   `add_to_newton(rhs,matrix)` después del ensamblaje interno;
5. impedir que esa rama cree `FEFaceValues`, consulte `nearest()` o multiplique
   por `JxW`;
6. implementar el ciclo provisional/aceptado y reconstrucción de restart;
7. verificar sobre el ejecutable productivo, no sólo sobre el arnés, fuerza
   unitaria/nula, tangente, diferencias finitas, mensajes inválidos y réplicas
   MPI;
8. mantener el aborto si falta cualquier mensaje o si los modos/hashes no
   coinciden.

Sólo después G2-B.2 puede implementar el handshake en dos fases y conectar el
fluido real. Este cambio no requiere A_f, masa virtual, IQN-ILS ni ajustes de
física/tolerancias.

## Archivos y evidencia

En este cierre bloqueado sólo se añade documentación. El código y los casos no
cambian. Comandos de auditoría:

```bash
git rev-parse HEAD
git status --short
rg -n 'GeneralizedInterfaceLoad|GeneralizedLoad|add_to_newton' src/dealiiPdmsSolid
rg -n 'dualConservative production transport|assemble_newton' \
  src/dealiiPdmsSolid/dealiiPdmsSolid.cc
```

No se inicia G3.

## Resolución del bloqueo

G2-B.1-RUNTIME conectó `GeneralizedLoad` a `PdmsSolid::assemble_newton()` y
validó el camino real en serial/MPI. Por tanto, el bloqueo descrito aquí ya no
impide reabrir G2-B.2 en un gate posterior. Este trabajo no reintentó handshake
ni modificó el acoplador fluido.
