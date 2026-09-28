# Nomenclatura IUPAC de compuestos orgánicos (JSME + PHP) · v4.2

Aplicación web didáctica de la Escuela de Química (UNACHI). Genera estructuras o recibe las que dibuja el estudiante en **JSME**, calcula su **nombre IUPAC preferido (PIN)** a partir del grafo molecular según las *Recomendaciones IUPAC 2013* y presenta las demás formas aceptadas: prefijos tradicionales y criterios de 1979/1993, prefijos sistemáticos de 1979, nombres de clase funcional y nombres comunes.

## Novedades de la v4.2 (fase 3A: heterociclos y azufre)

- **Dos casillas nuevas**: «Heterociclos» y «Azufre». Combinada con «Dos anillos», la primera genera heterociclos bicíclicos.
- **Lector de estructuras**: acepta S y los heteroátomos aromáticos de JSME (n, o, s, [nH]). Cada sistema aromático se convierte en su estructura de Kekulé; el N–H de pirroles e imidazoles, el O y el S no participan en los dobles enlaces.
- **Heteromonociclos con nombre retenido** (P-22.2.1):
  - insaturados: pirrol, furano, tiofeno, imidazol, pirazol, piridina, piridazina, pirimidina, pirazina, pirano, tiopirano;
  - saturados: pirrolidina, pirazolidina, imidazolidina, piperidina, piperazina, morfolina.
- **Nomenclatura de Hantzsch-Widman** (P-22.2.2) para anillos de 3 a 10 átomos: `oxirano`, `aziridina`, `tiirano`, `oxetano`, `azetidina`, `oxolano`, `tiolano`, `1,3-dioxolano`, `oxano`, `1,4-dioxano`, `1,3-oxazol`, `1,2-oxazol`, `1,3-tiazol`, `1,3,5-triazina`, `1H-1,2,4-triazol`, `azepano`.
- **Numeración de los heteroátomos**: reciben los localizadores más bajos y, si hay elección, primero O, luego S y luego N (P-31.1.4.2.2).
- **Hidrógeno indicado y añadido en heterociclos**: `2H-pirano`, `piridin-2(1H)-ona`, `pirimidina-2,4(1H,3H)-diona` (uracilo), `1,3,7-trimetil-3,7-dihidro-1H-purina-2,6-diona` (cafeína).
- **Heterociclos fusionados con nombre retenido**:
  - indol, isoindol, 1-benzofurano, 2-benzofurano, 1-benzotiofeno, 1H-bencimidazol, 1H-indazol, 1,3-benzoxazol, 1,2-benzoxazol, 1,3-benzotiazol, purina;
  - quinolina, isoquinolina, quinazolina, quinoxalina, cinolina, ftalazina, 1-benzopirano, 2-benzopirano, 1-benzotiopirano, 1,4-benzodioxina.
- **Heterociclos parcialmente hidrogenados**: prefijos hidro sobre el nombre insaturado, como en `2,3-dihidro-1H-indol`, `1,2,3,4-tetrahidroquinolina`, `3,4-dihidro-2H-1-benzopirano` y `3,4-dihidro-2H-pirano`.
- **Prefijos de reemplazo en puentes y espiro** (P-15.4): `1-azabiciclo[2.2.2]octano` (quinuclidina), `7-oxabiciclo[2.2.1]heptano`, `1,4-dioxaespiro[4.5]decano`.
- **Lactonas y lactamas** como cetonas del heterociclo: `oxolan-2-ona`, `pirrolidin-2-ona`, `azepan-2-ona`, `2H-1-benzopiran-2-ona` (cumarina).
- **N-acil heterociclos**, nombrados como cetonas o con sufijos de anillo: `1-(piperidin-1-il)etan-1-ona`, `piperidina-1-carbaldehído`, `pirrolidina-1-carboxilato de etilo`.
- **Antigüedad de anillos con heteroátomos** (P-44.2.1): un heterociclo es preferido a un carbociclo, y uno con N a los demás. Por ejemplo, `2-fenilpiridina`, `3-(1-metilpirrolidin-2-il)piridina` (nicotina) y `2,2'-bipiridina`.
- **Compuestos de azufre**:
  - tioles: sufijo `-tiol`, prefijo `sulfanil` (`etanotiol`, `bencenotiol`);
  - sulfuros: `(metilsulfanil)metano`;
  - sulfóxidos y sulfonas: `(metanosulfinil)metano` (DMSO), `(metanosulfonil)metano`;
  - ácidos sulfónicos, que en el orden de prioridad van justo después de los carboxílicos: `ácido bencenosulfónico`, `ácido 4-metilbenceno-1-sulfónico`.
- **Nombres de clase funcional y comunes**: sulfuro de dimetilo, DMSO, THF, cafeína, uracilo, adenina, nicotinamida, furfural, cumarina, quinuclidina, tropano, tiofenol, ácido p-toluenosulfónico.

## Novedades de la v4.1 (fase 2: dos anillos, fusionados o no)

- **Casilla nueva «Dos anillos»**: el generador produce naftalenos, indanos, tetralinas, decalinas, norbornanos, biciclo[2.2.2]octanos, espiranos, bifenilos, difenilmetanos, difenil éteres y otros esqueletos, combinables con todas las familias de grupos funcionales.
- **Biciclos orto-fusionados con nombre retenido** (P-25.1): pentaleno, indeno, azuleno, naftaleno y heptaleno, con su numeración fija (3a, 4a, 8a…). El programa elige la orientación y aplica:
  - hidrógeno indicado: `1H-indeno`;
  - prefijos hidro: `2,3-dihidro-1H-indeno`, `1,2,3,4-tetrahidronaftaleno`, `decahidronaftaleno`, `2,3,3a,4,5,6,7,7a-octahidro-1H-indeno`;
  - hidrógeno añadido: `3,4-dihidronaftalen-1(2H)-ona`.
- **Biciclos con puente (von Baeyer, P-23.2)**: `biciclo[2.2.1]heptano`, `biciclo[2.2.1]hept-2-eno`, `1,7,7-trimetilbiciclo[2.2.1]heptan-2-ona` (alcanfor), `biciclo[4.1.0]heptano`. También los biciclos fusionados con un anillo de 3 o 4 miembros, como `biciclo[4.2.0]octa-1,3,5-trieno`, con elección de la estructura de Kekulé.
- **Espiro (P-24.2)**: `espiro[4.5]decano`, `espiro[4.5]decan-8-ona`.
- **Ensamblajes de anillos (P-28)**: `1,1'-bifenilo`, `[1,1'-bifenil]-4-ol`, `1,1'-bi(ciclohexano)`, `[1,1'-binaftaleno]-2,2'-diol`.
- **Nomenclatura multiplicativa (P-15.3)** para unidades idénticas: `1,1'-metilendibenceno`, `1,1'-oxidibenceno`, `4,4'-(propano-2,2-diil)difenol` (bisfenol A), `4,4'-metilendianilina`, `ácido 4,4'-carbonildibenzoico`, `1,1'-metilenbis(4-clorobenceno)`.
- **Antigüedad de anillos (P-44.2, P-44.4)**: más anillos, luego más átomos, luego menor hidrogenación (`1-bencilnaftaleno`, `fenilciclooctano`, `ciclohexilbenceno`).
- **Sistemas de anillos como sustituyentes**: `naftalen-2-il` (tradicional: `2-naftil`), `2,3-dihidro-1H-inden-1-il`, `biciclo[2.2.1]heptan-2-il`, `naftaleno-1-carbonil`.
- **Numeración en JSME**: el dibujo muestra los localizadores enteros. Los de fusión (4a, 8a) y los que llevan prima se explican en el texto, porque los mapas de átomo de JSME solo admiten números enteros.

## Novedades de la v4.0 (fase 1)

- **Nuevo motor** `lib/IupacOrganica.php`: nomenclatura sustitutiva general con grupo principal, sufijos y prefijos.
- **Familias**: alquenos, alquinos, benceno y sus derivados, alcoholes, fenoles, éteres, aminas, aldehídos, cetonas, ácidos carboxílicos, ésteres, amidas, nitrilos, nitroderivados y halogenuros.
- **Casillas por tipo de compuesto**, en el mismo estilo de trabajo que la casilla *Incluir halógenos* de la v3.0. El navegador recuerda la selección.
- **Orden de prioridad de clases (P-41)**: el grupo principal va como sufijo y los demás como prefijos. Los niveles intermedio y avanzado combinan grupos para practicarlo.
- **Criterio de 2013 para la cadena principal** (primero la longitud y después la insaturación), con una nota explicativa cuando las reglas de 1979/1993 elegirían otra cadena. Ejemplo: `3-metilidenhexano` (2013) frente a `2-etilpent-1-eno` (1993).
- **Nombres retenidos preferidos**: fenol, anilina, ácido benzoico, benzaldehído, benzamida, benzonitrilo, ácido fórmico, ácido acético, ácido oxálico, formaldehído, acetaldehído, acetamida, acetonitrilo, tolueno, xileno, anisol y acetileno.
- **Verificación ampliada**. Acepta el PIN, las formas tradicional y 1979, los nombres de clase funcional (`alcohol isopropílico`, `etil metil éter`, `etil(metil)amina`) y los nombres comunes (`acetona`, `ácido acético`, `estireno`), en español o en inglés. También detecta el formato anterior a 1993 (`2-butanol` → `butan-2-ol`).
- **Regresión**: en el dominio de la v3.0 (alcanos, cicloalcanos y haloalcanos) el motor nuevo produce los mismos nombres que `IupacAlcanos.php`. La única diferencia es intencional: ya no se escribe guion entre prefijos sin localizador (`bromoclorometil`, no `bromo-clorometil`).

## Funcionamiento

| Botón | Acción |
|---|---|
| **Generar** | `api.php` crea una estructura aleatoria según el esqueleto, el nivel y las familias marcadas, calcula sus nombres y los guarda **solo en la sesión PHP**. Al navegador se envían únicamente el SMILES y la fórmula. |
| **Dibujar** | Limpia el editor y activa el modo de edición de JSME. |
| **Verificar** | Compara la respuesta con todos los nombres válidos, con retroalimentación graduada: correcto (indica qué forma se usó), puntuación, localizadores o incorrecto. |
| **Desplegar** | Muestra el PIN, las otras formas aceptadas, la equivalencia de los sustituyentes, el nombre en inglés, la justificación paso a paso con la regla IUPAC de cada decisión y la numeración sobre la estructura. |

### Niveles

| Nivel | Cadena | Anillo | Grupos funcionales |
|---|---|---|---|
| Básico | 2–6 C, 0–1 ramas | 3–6 C | 1 |
| Intermedio | 3–8 C, hasta 2 ramas | 4–7 C | 1–2 |
| Avanzado | 4–12 C, 1–3 ramas | 5–8 C | 2–3 (combinados) |

Sin casillas marcadas se generan alcanos y cicloalcanos. *Aromáticos* o *Fenoles* sustituyen el anillo por un benceno o añaden un fenilo a la cadena.

## Reglas implementadas (IUPAC 2013)

| Aspecto | Regla | Ejemplo |
|---|---|---|
| Orden de prioridad de clases | P-41 | ácido > éster > amida > nitrilo > aldehído > cetona > alcohol/fenol > amina |
| Máximo número de grupos principales | P-44.1.1 | `4-(hidroximetil)fenol` |
| Anillo preferido a la cadena | P-44.1.2.2 | `decilbenceno`, `fenilmetanol` |
| Cadena más larga; después, más enlaces múltiples y más dobles enlaces | P-44.3, P-44.4.1 | `3-metilidenhexano` |
| Numeración: grupo principal → enlaces múltiples → dobles enlaces → prefijos → primer prefijo citado | P-31.1.4 | `hex-3-en-1-ino`, `pent-3-en-2-ol` |
| Sufijos en anillo con carbono externo | P-65, P-66 | `ácido ciclohexanocarboxílico`, `ciclohexanocarbaldehído` |
| Aldehído o cetona no principal: oxo si el carbono está en la cadena y formil si no lo está | P-64, P-66.6 | `ácido 3-oxopropanoico`, `2-formilbenzamida` |
| Ésteres | P-65.6 | `propanoato de etilo`, `pentanodioato de 1-etilo y 5-metilo` |
| Aminas y amidas con localizador N | P-62, P-66.1 | `N,N-dietiletanamina`, `N-fenilacetamida` |
| Éteres como prefijos alcoxi | P-63.2 | `2-metoxi-2-metilpropano`, `(propan-2-il)oxi` |
| Omisión de localizadores | P-14.3.4 | `etanol`, `propeno`, `ciclohexeno`, `3-metilciclohex-1-eno` |
| Orden alfanumérico y multiplicadores | P-14.5, P-16.3 | `N,N,3-trimetil…`, `bis(2-metilpropil)` |

## Estructura

```
nomenclatura/
├── index.php, server.php   Interfaz (emite el token CSRF y la sesión)
├── api.php                  Servicio JSON: generar / analizar / desplegar / verificar
├── lib/IupacOrganica.php    Motor v4.0: lector SMILES, grupos, estructura principal, nombres, explicación, generador
├── lib/IupacAlcanos.php     Motor v3.0 (solo para la prueba de regresión)
├── assets/app.js            Cliente e integración con JSME
├── assets/estilos.css       Estilos adaptables
├── jsme/                    JSME 2024.04.29 (BSD-3), local, sin CDN
└── tests/                   Validación (acceso web bloqueado)
```

## Instalación

Requisitos: PHP ≥ 7.4 (probado en 8.3 y 8.4) con `mbstring` y sesiones. No usa base de datos. El `Dockerfile` no cambia (Dokploy en el VPS de Hostinger): basta con hacer *push* y volver a desplegar. Con Apache, los `.htaccess` bloquean `lib/` y `tests/`. Con Nginx:

```nginx
location ~ ^/(lib|tests)/ { deny all; return 404; }
```

## Validación realizada

Con **OPSIN** (nombre → estructura) y **RDKit** (SMILES canónico), sobre los nombres en inglés en los estilos PIN, tradicional, 1979 y de clase funcional:

| Prueba | Estructuras | Nombres | Discrepancias |
|---|---|---|---|
| `validar_org.py gen`: estructuras del generador, 3 niveles y combinaciones aleatorias de familias | 2995 | 9148 | **0** |
| `aleatorias_org.py`: estructuras arbitrarias como las que dibujaría un estudiante (cadenas, un anillo o benceno, insaturaciones y 0–3 grupos), 4 semillas | 9000 | 27 207 | **0** |
| Invariancia ante el orden de los átomos (SMILES aleatorios y Kekulé) | 7200 | — | **0** |
| `esperados_es.tsv`: 44 compuestos de referencia con su nombre en español (aspirina, ibuprofeno, paracetamol, ácido cinámico…) | 44 | 44 | **0** |
| **Fase 2** · `aleatorias_bic.py`: estructuras de dos anillos con grupos funcionales (2 semillas) | 3 756 | 11 330 | **0** * |
| **Fase 2** · invariancia ante el orden de los átomos y la forma de Kekulé | 11 268 | — | **0** |
| **Fase 2** · `esperados_fase2.tsv`: 34 compuestos de referencia en español (tetralona, alcanfor, bisfenol A, BINOL…) | 34 | 34 | **0** |
| **Fase 3A** · `aleatorias_het.py`: heterociclos y compuestos de azufre con grupos funcionales (3 semillas) | 5 714 | 17 176 | **0** |
| **Fase 3A** · invariancia ante el orden de los átomos y la forma de Kekulé | 17 142 | — | **0** |
| **Fase 3A** · `esperados_fase3a.tsv`: 43 compuestos de referencia en español (cafeína, nicotina, uracilo, adenina, DMSO, cumarina, quinuclidina…) | 43 | 43 | **0** |

\* En 157 nombres de pentaleno y heptaleno, OPSIN devuelve otra estructura de Kekulé del mismo compuesto. Estos sistemas no son aromáticos para RDKit y el nombre IUPAC no distingue entre sus estructuras de Kekulé, así que la prueba los acepta como equivalentes.
| `regresion_alcanos.php`: igualdad con el motor v3.0 (6 variantes por estructura) | 3000 | 18 000 | **0** |

```bash
pip install rdkit pyopsin
export OPSIN_JAR=$(python3 -c "import pyopsin,os;print(os.path.join(os.path.dirname(pyopsin.__file__),'opsin_cli.jar'))")
cd tests
python3 validar_org.py gen 2000 1
python3 aleatorias_org.py 2000 1
python3 aleatorias_bic.py 1500 1        # fase 2
python3 aleatorias_het.py 1500 1        # fase 3A
php regresion_alcanos.php 1000
```

OPSIN comprueba que cada nombre corresponde exactamente a la estructura. La elección del nombre *preferido* la garantizan los criterios de la búsqueda exhaustiva del motor, que se comprobaron con los 44 casos de referencia.

## Alcance y limitaciones (v4.2)

- **Esqueleto**: cadenas de hasta 20 C; anillos simples; sistemas de dos anillos fusionados, con puente o espiro; ensamblajes y unidades multiplicadas de dos sistemas idénticos. Puede haber cualquier número de sistemas de anillos como sustituyentes.
- **Pendiente (fases 3B y 3C)**:
  - estereodescriptores E/Z y R/S, que se ignoran si se dibujan;
  - haluros de acilo, anhídridos, imidas, iminas, oximas, carbamatos y ureas;
  - sistemas de tres o más anillos;
  - fusionados sin nombre retenido (p. ej., `benzo[7]anuleno`, furo[3,2-b]piridina, naftiridinas);
  - heterociclos con S oxidado en el anillo (sulfolano), disulfuros, sulfonamidas y ésteres sulfónicos.
- **Sulfóxidos y sulfonas**: los prefijos `alcanosulfinil` y `alcanosulfonil` se construyen para grupos alquilo, cicloalquilo, arilo y heteroarilo. Para grupos R ramificados poco comunes (p. ej., isobutilo) en el nombre preferido, el programa lo informa.
- **Tropano**: el programa escribe `8-metil-8-azabiciclo[3.2.1]octano`, con el localizador del anillo en el sustituyente del N.
- **Nomenclatura multiplicativa**: se aplica cuando las dos unidades son idénticas, con los mismos prefijos. Si llevan prefijos distintos, el programa usa nomenclatura sustitutiva, que es un nombre válido. Conviene revisar esos casos con las reglas P-15.3.
- **Ensamblajes con componentes parcialmente hidrogenados** o con hidrógeno indicado (p. ej., biindeno): se nombran de forma sustitutiva.
- Con más de dos grupos –COOH, –CHO o –CN en una cadena acíclica se usan prefijos carboxi, formil o ciano en lugar de la excepción P-65.1.2.2.
- `acetileno` y `anisol` se usan solo sin sustituir.

## Referencias

- Favre, H. A.; Powell, W. H. *Nomenclature of Organic Chemistry. IUPAC Recommendations and Preferred Names 2013*. Cambridge: RSC, 2014. doi:10.1039/9781849733069 (versión en línea: iupac.qmul.ac.uk/BlueBook).
- IUPAC. *Nomenclature of Organic Chemistry, Sections A–H* (1979). Oxford: Pergamon, 1979.
- IUPAC. *A Guide to IUPAC Nomenclature of Organic Compounds (Recommendations 1993)*. Oxford: Blackwell, 1993.
- Bienfait, B.; Ertl, P. JSME: a free molecule editor in JavaScript. *J. Cheminform.* 2013, 5, 24. doi:10.1186/1758-2946-5-24
- Lowe, D. M. et al. Chemical name to structure: OPSIN. *J. Chem. Inf. Model.* 2011, 51, 739–753. doi:10.1021/ci100384d
