# G2-B.1-RUNTIME — integración Newton real

Fecha: 2026-08-26

Commit de entrada: `11fc7942615e3e9da1b08727a1e4aa4394baa7c0`

Estado: **PASS; G2-B.2 queda pendiente**

El subgate posterior G2-B.1-EXEC-PROTOCOL expuso esta misma ruta mediante el
binario separado, con socket versionado y participante externo. Véase
`docs/G2-B.1-executable-protocol.md`.

## Ruta y signo

`PdmsSolid::assemble_newton()` comienza cada evaluación con `matrix=0` y
`rhs=0`. Ensambla material neo-Hookean e inercia BE/BDF2 por celda y aplica
`AffineConstraints::distribute_local_to_global()`. El Newton resuelve
`J delta_d=-R`; por ello la implementación histórica construye `rhs=-R`.

Las rutas de interfaz son excluyentes:

- `legacyNearestNeighbour`: conserva `FEFaceValues`, `nearest()` y `JxW`;
- `dualConservative`: omite completamente ese bucle y, después del ensamblaje
  interno, ejecuta una sola vez `GeneralizedLoad::add_to_newton(rhs,matrix)`.

La convención B aprobada queda materializada: `R=R_s-f_Gamma`, por lo que
`f_Gamma [N]` se suma a `rhs=-R`; `J_Gamma=-df_Gamma/dd [N/m]` se suma a J.
La ruta dual no integra, multiplica por área ni busca vecinos.

## Estado runtime

`PdmsSolid` posee manifiesto/mapeo de DoFs, componentes, `Expected`, carga y
tangente provisionales y `GeneralizedLoad`. Los setters separados validan modo,
sellos `(timeIndex,outerCorrector,operatorVersion)`, versiones Z, hashes,
unidades, checksums, IDs/componentes y duplicados. Carga y tangente sólo se
activan cuando ambas existen y coinciden. La validez es explícita; cero es una
carga válida.

Los mensajes permanecen disponibles para todos los Newton del corrector. Como
matrix/rhs se reinician, repetir `assemble_newton()` no acumula contribuciones.
`clear_provisional_generalized_interface_state()` elimina ambos mensajes.
`begin_generalized_interface_corrector()` exige el mismo hashGraph, limpia el
estado anterior y requiere mensajes nuevos; cambiar el soporte después de
`reinit()` aborta.

No se serializan mensajes provisionales. Un restart reconstruye manifiesto,
patrón y hashes y carece de carga/tangente hasta una recepción nueva.

## Orden dual

1. distribuir DoFs;
2. construir constraints;
3. leer/verificar manifiesto;
4. ampliar `DynamicSparsityPattern`;
5. crear `SparsityPattern` y ejecutar `SparseMatrix::reinit()`;
6. crear el estado runtime esperado;
7. recibir fuerza y tangente;
8. permitir `assemble_newton()`.

El ejecutable standalone aún aborta dual porque no tiene transporte vivo. La
API real permite inyección programática versionada; dual sin datos aborta antes
del ensamblaje. Ningún caso FSI fue cambiado.

## Prueba sobre assemble_newton real

`test_runtime_newton.cc` incluye la implementación productiva y llama los
métodos reales de `PdmsSolid`. El participante inyecta una fuerza unitaria y
una tangente diagonal sobre el DoF libre 12.480.

| ejecución | rhs del DoF | error repetición | mejor error FD central | error con signo invertido |
|---|---:|---:|---:|---:|
| serial | 1 | 0 | 1,02721e-11 | 2 |
| MPI-2 | 1 | 0 | 1,02721e-11 | 2 |
| MPI-4 | 1 | 0 | 1,02721e-11 | 2 |

`rhs=+1` confirma el signo y ausencia de integración adicional. El valor no se
multiplica por ranks. Cada rank construye la misma matriz serial replicada una
vez. La diferencia finita usa el residual estructural completo real y la ley
coherente `f(d)=f0-J_Gamma d`; el error es inferior a 1e-8/1e-7.

También se verifican carga nula válida, estado sin tangente, mensaje obsoleto,
hash/versiones mediante la API base, transición de corrector, constraints,
patrón y llamadas Newton repetidas.

## Regresiones

- build del sólido: PASS;
- GeneralizedLoad y patrón G2-B.1a serial/MPI: PASS;
- BE/BDF2, signo y material: PASS;
- `constantHs`: PASS;
- mínimo/restart legacy continuo-partido serial/MPI: PASS;
- equivalencia serial/MPI `rtol=5e-5`: PASS;
- sin NaN/Inf;
- patrón legacy: 322.992 nonzeros, sin cambios;
- SHA legacy idénticos: serial `8bc7ad2d...80e0`/`cceaabe1...cc1`, MPI
  `7bc98c11...a840`/`b6e5fa56...2f4a`.

Artefactos: `/tmp/g2b1-runtime-final-20260826` y
`/tmp/g2b1-runtime-legacy-restart-20260826`. La campaña final del arnés,
incluidas las pruebas negativas de hash, mensaje obsoleto, tangente ausente,
signo invertido y fuerza nula válida, quedó en
`/tmp/g2b1-runtime-final2-20260826`.

## Decisión

G2-B.1-RUNTIME y G2-B.1 quedan **PASS**. La condición previa que bloqueó
G2-B.2 queda corregida, pero G2-B.2 no fue reintentado. Para reabrirlo debe
implementarse el handshake/transporte del fluido hacia estos setters sin
archivos temporales, manteniendo sellos y hashes. No se inició G3.
