# Nomenclatura IUPAC de compuestos orgánicos (JSME + PHP) · v4.0

Aplicación web didáctica de la Escuela de Química (UNACHI). Genera estructuras o recibe las que dibuja el estudiante en **JSME**, calcula su **nombre IUPAC preferido (PIN)** a partir del grafo molecular según las *Recomendaciones IUPAC 2013* y presenta las demás formas aceptadas: prefijos tradicionales y criterios de 1979/1993, prefijos sistemáticos de 1979, nombres de clase funcional y nombres comunes.

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
| `regresion_alcanos.php`: igualdad con el motor v3.0 (6 variantes por estructura) | 3000 | 18 000 | **0** |

```bash
pip install rdkit pyopsin
export OPSIN_JAR=$(python3 -c "import pyopsin,os;print(os.path.join(os.path.dirname(pyopsin.__file__),'opsin_cli.jar'))")
cd tests
python3 validar_org.py gen 2000 1
python3 aleatorias_org.py 2000 1
php regresion_alcanos.php 1000
```

OPSIN comprueba que cada nombre corresponde exactamente a la estructura. La elección del nombre *preferido* la garantizan los criterios de la búsqueda exhaustiva del motor, que se comprobaron con los 44 casos de referencia.

## Alcance y limitaciones de la fase 1

- **Esqueleto**: cadenas de hasta 20 C, un anillo carbocíclico (3–20 C, saturado, insaturado o benceno) como estructura principal y cualquier número de anillos aislados como sustituyentes cuando la estructura principal es una cadena (`difenilmetanona`).
- **Pendiente para la fase 2**: anillos fusionados (naftaleno, indano, tetralina, decalina), con puente (norbornano), espiro, ensamblajes (bifenilo) y nomenclatura multiplicativa (`1,1'-metilendibenceno`).
- **Pendiente para la fase 3**: heterociclos (piridina, furano, lactonas), haluros de acilo, anhídridos, imidas, iminas, tioles y estereodescriptores E/Z y R/S (si se dibujan, se ignoran).
- Con más de dos grupos del tipo –COOH, –CHO o –CN sobre una cadena acíclica, el programa usa prefijos carboxi, formil o ciano en lugar de la excepción P-65.1.2.2 (`propano-1,2,3-tricarboxílico`).
- `acetileno` y `anisol` se usan solo sin sustituir. Sus derivados se nombran de forma sistemática (`1-metoxi-4-nitrobenceno`).

## Referencias

- Favre, H. A.; Powell, W. H. *Nomenclature of Organic Chemistry. IUPAC Recommendations and Preferred Names 2013*. Cambridge: RSC, 2014. doi:10.1039/9781849733069 (versión en línea: iupac.qmul.ac.uk/BlueBook).
- IUPAC. *Nomenclature of Organic Chemistry, Sections A–H* (1979). Oxford: Pergamon, 1979.
- IUPAC. *A Guide to IUPAC Nomenclature of Organic Compounds (Recommendations 1993)*. Oxford: Blackwell, 1993.
- Bienfait, B.; Ertl, P. JSME: a free molecule editor in JavaScript. *J. Cheminform.* 2013, 5, 24. doi:10.1186/1758-2946-5-24
- Lowe, D. M. et al. Chemical name to structure: OPSIN. *J. Chem. Inf. Model.* 2011, 51, 739–753. doi:10.1021/ci100384d
