# G2-B.1-EXEC-PROTOCOL — protocolo del ejecutable estructural

## Dependencia pendiente para el handshake del operador

El transporte permanece PASS, pero no se añadieron todavía
`FluidPointManifest/HpsRowManifest`. La auditoría determinó que la referencia
G2-A no ofrece una evaluación Q1 deal.II reutilizable: agregar sólo los schemas
ocultaría la ausencia del productor canónico de filas. El bloqueo pertenece al
constructor de operador, no al framing, checksum o difusión MPI de este gate.

## Extensión cinemática requerida por G2-B.2-LIVE

El protocolo expone ahora `StructuralStateMessage` mediante una solicitud
versionada, usando el mismo canal Unix y la difusión colectiva rank-0. El
payload canónico contiene estado aceptado/provisional, BE/BDF2, unidades y
`(global_dof_id, componente, d, v, a)`. Las pruebas externas serial/MPI-2/MPI-4
pasan. Se añadieron mensajes explícitos de aceptación, rechazo y rollback; no
se persiste estado provisional.

Esta extensión no cambia el PASS de EXEC-PROTOCOL. La auditoría LIVE detectó
separadamente que `DofManifest` aún no contiene conectividad/soportes Q1
suficientes para construir `H_ps`; esa limitación bloquea G2-B.2, no la entrega
de fuerza/tangente ya aprobada aquí.

Fecha: 2026-08-26

Commit de entrada: `d77d9e6a493ca74436bdd849a918836bceb20889`

Estado: **PASS; G2-B.2 permanece pendiente**

## Transporte

Se añadió un socket Unix bidireccional local al ejecutable
`dealiiPdmsSolid`. Es coherente con el modelo existente de participantes como
procesos separados y evita compartir objetos PETSc entre las ABI distintas de
OpenFOAM y deal.II. No usa polling ni ficheros de mensajes. `protocol endpoint`
es explícito; rank 0 crea/acepta el socket y aplica un timeout configurable.
Ausencia de participante, lectura parcial o desconexión producen fallo
diagnóstico y eliminación del socket propio.

El participante simulado es un proceso Python independiente. No enlaza
deal.II, no incluye `PdmsSolid` y utiliza exactamente el framing que podrá usar
el futuro adaptador OpenFOAM.

## Sobre canónico

Cada frame consta de una línea ASCII y exactamente `payloadLength` bytes:

```text
G2B1_EXEC schemaVersion messageType sequenceNumber payloadLength
payloadChecksum producerId transferMode timeIndex outerCorrector
operatorVersion messageChecksum\n
payload
```

Los enteros son decimales sin signo y los escalares de payload se serializan
con 17 dígitos. Al ser texto no existe ambigüedad de endianness. Ambos
checksums son FNV-1a de 64 bits: uno cubre el payload y el otro la cabecera
canónica sin `messageChecksum`, el separador y el payload. El límite es 64 MiB.
Se rechazan magic/esquema/tipo, longitud, checksums, truncamiento, trailing
bytes, NaN/Inf y entradas no canónicas antes de mutar el estado.

Tipos implementados: `Hello`, `Capabilities`, `DofManifest`,
`OperatorManifest`, `ForceMessage`, `TangentMessage`,
`ActivateCorrectorState`, `ClearProvisionalState`, `Ready`, `Ack`, `Nack`,
`AssemblyResult` y `Shutdown`.

## Handshake y patrón

La máquina de estados ejecuta:

1. `AwaitingHello`: acuerda dimensión 3, scalar de 8 bytes, unidades N/N/m,
   signo `R=Rs-fGamma;J=Js+JGamma`, BE/BDF2 y matriz replicada;
2. distribuye DoFs y construye `AffineConstraints` y el patrón FEM base;
3. exporta `DofManifest` canónico con ID, componente, coordenada, owner lógico
   0, restricción y `dofManifestHash`;
4. recibe el `OperatorManifest` completo, con H/W, hashes, versión, unidades,
   signo y longitud explícita;
5. amplía el `DynamicSparsityPattern` antes de copiar/reinicializar la matriz;
6. responde `Ready` con los tres hashes;
7. sólo entonces acepta estados de corrector.

El refactor separa `prepare_setup()` de `finish_setup()`. Legacy llama ambas de
forma consecutiva y no requiere manifiesto ni socket.

## Activación, replay y MPI

Fuerza y tangente se deserializan y validan por separado, pero permanecen en
`pendingForce/pendingTangent`. `ActivateCorrectorState` referencia los
checksums exactos de ambos frames y sus sellos/hashes. Sólo después llama
atómicamente `begin_generalized_interface_corrector()` y los setters runtime
reales. Fuerza cero es válida. `ClearProvisionalState` exige la transición
`outer+1` o nuevo tiempo/outer 1; no se reutiliza el estado anterior.

Un `sequenceNumber` repetido con el mismo checksum recibe ACK idempotente; con
contenido diferente recibe NACK. Correctores obsoletos/futuros, versiones
distintas y cambios de `hashGraph` se rechazan.

Sólo rank 0 lee el socket. Difunde una vez los bytes canónicos; cada rank
deserializa los mismos bytes y ensambla una vez en su matriz serial replicada.
Los errores de I/O se difunden antes de salir para evitar ranks bloqueados. El
ACK ocurre después de la validación/activación colectiva. No hay reducción de
fuerza en el sólido.

## Resultados

El proceso externo inyectó `DeltaR` de norma 1 N y `DeltaJ` de norma 100 N/m,
una entrada, y ejecutó dos ensamblajes reales por corrector:

| ranks | `||DeltaR||` | `||DeltaJ||` | error repetición R | error repetición J |
|---:|---:|---:|---:|---:|
| 1 | 1 | 100 | 0 | 0 |
| 2 | 1 | 100 | 0 | 0 |
| 4 | 1 | 100 | 0 | 0 |

También pasó replay idempotente, `ClearProvisionalState`, transición al
corrector siguiente y cierre/reapertura. Por la misma frontera de proceso se
probaron fuerza distribuida (`||DeltaR||=sqrt(5)`), tangente diagonal variable
(`||DeltaJ||=sqrt(50000)`), fuerza/tangente exactamente nulas válidas y un DoF
restringido correctamente eliminado por `AffineConstraints`, en
serial/MPI-2/MPI-4. La regresión G2-B.1-RUNTIME conserva
error FD `1,02721e-11`, error de repetición 0 y error 2 con signo invertido en
serial/MPI-2/MPI-4.

Las 27 pruebas negativas cubren magic, esquema, tipo, longitud, checksum,
truncamiento, timeout, NaN/Inf, IDs/duplicados, sparsity, unidades, los tres
hashes, correctores obsoleto/futuro, replay conflictivo, fuerza/tangente
incompletas, versiones distintas, activación incorrecta y desconexión parcial.
Todas produjeron NACK o cierre controlado sin ensamblaje parcial.

CTest de material, BE/BDF2, signo y framing: 4/4 PASS. `constantHs`: PASS.
Legacy continuo y 5+5 serial/MPI: PASS, `rtol=5e-5`, sin NaN/Inf. SHA-256:

- serial: `8bc7ad2d...80e0`, `cceaabe1...cc1`;
- MPI: `7bc98c11...a840`, `b6e5fa56...2f4a`.

## Archivos y continuidad

Código: `ExecutableProtocol.H`, `dealiiPdmsSolid.cc`,
`InterfaceSparsityExtension.H`, CMake y test C++. Arnés externo:
`g2b1_exec_participant.py`, `g2b1_exec_negative.py` y
`test-g2b1-exec-protocol.sh`. Campañas principales:
`/tmp/g2b1-exec-final-20260826`,
`/tmp/g2b1-exec-runtime-regression-20260826` y
`/tmp/g2b1-exec-legacy-restart2-20260826`.
El binario final instalado en el prefijo aislado tiene SHA-256
`02eb99287b0c03376a467ec0a860b5c9b9d3d00f89fb14c9fec1352f87280ec4`.

`src/robinRobinCoupling` y los casos FSI permanecen intactos y en legacy.
G2-B.1-EXEC-PROTOCOL es PASS. La ruta estructural ya está lista para reabrir
G2-B.2, que deberá implementar el mismo cliente en el acoplador. No se reintentó
G2-B.2 ni se inició G3.
