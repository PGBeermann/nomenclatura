# Nomenclatura IUPAC de alcanos, cicloalcanos y haloalcanos (JSME + PHP) · v3.0

Aplicación web didáctica que genera estructuras de alcanos lineales, ramificados y cicloalcanos (con o sin ramificaciones), o recibe las que dibuja el estudiante en **JSME**. Calcula su **nombre IUPAC** a partir del grafo molecular, según las *Recomendaciones IUPAC 2013*, y presenta además las formas alternativas de nombrar los sustituyentes que aceptan las reglas de 1979/1993.

## Novedades de la v3.0

- **Haloalcanos (F, Cl, Br, I)**: la casilla *Incluir halógenos* agrega de 1 a 4 halógenos según el nivel. En el modo *Dibujar*, el estudiante también puede usar los botones F, Cl, Br e I de JSME.
- **Nombre radicofuncional** para monohaloalcanos con grupo R sencillo: `bromuro de isopropilo`, `cloruro de terc-butilo`, `cloruro de ciclohexilo`, junto con su versión inglesa (`isopropyl bromide`).
- **Prefijos sec- y neo-** en la verificación: se aceptan las variantes de escritura frecuentes (`s-butil`, `neo-pentil`, `iso-propil`, `iodo`/`yodo`), que se normalizan antes de comparar.

## Funcionamiento

| Botón | Acción |
|---|---|
| **Generar** | El servidor (`api.php`) crea una estructura aleatoria según el tipo y el nivel elegidos, calcula sus nombres y los guarda **solo en la sesión PHP**. Al navegador se envía únicamente el SMILES y la fórmula. |
| **Dibujar** | Limpia el editor y activa el modo de edición de JSME. El estudiante dibuja su propio alcano o cicloalcano, escribe el nombre y pulsa **Verificar** o **Desplegar**. |
| **Verificar** | Compara la respuesta del estudiante con todos los nombres válidos (preferido 2013, retenido y sistemático 1979, en español o en inglés). La retroalimentación es graduada: correcto (indica qué forma usó), error de puntuación o error de localizadores. |
| **Desplegar** | Muestra el nombre IUPAC preferido y las otras formas aceptadas, por ejemplo `(propan-2-il)` / `isopropil` / `(1-metiletil)`. Incluye una tabla de equivalencias de los sustituyentes, el nombre en inglés, la justificación paso a paso con la referencia a cada regla y la numeración dibujada sobre la estructura. |

**Estructuras dibujadas o modificadas.** Cada vez que se pulsa *Verificar* o *Desplegar*, el programa compara lo que hay en el editor con la última estructura cargada. Si el estudiante la cambió, la envía al servidor (acción `analizar`) y la nombra. Si no es un alcano, el estudiante recibe un mensaje explicativo: enlaces dobles o triples, heteroátomos, anillos aromáticos, más de un anillo, más de una molécula, cadena principal de más de 20 C o carbonos con más de cuatro enlaces.

### Niveles

| Nivel | Cadena principal | Ramificaciones | Sustituyentes |
|---|---|---|---|
| Básico | 4–7 C | 1–2 | metil, etil |
| Intermedio | 5–10 C | 1–3 | + propil, isopropil, butil, sec‑butil, isobutil, terc‑butil |
| Avanzado | **8–20 C** | 2–5 | + pentil, isopentil, neopentil, terc‑pentil, pentan‑2‑il, pentan‑3‑il, 2‑metilbutil, 3‑metilbutan‑2‑il, 3,3‑dimetilbutan‑2‑il, 2,3‑dimetilbutil, 2‑metilpentan‑3‑il, 3,3‑dimetilbutil, 2,3‑dimetilbutan‑2‑il |

En el nivel avanzado, los cicloalcanos van de 5 a 12 C. Como la cadena principal se determina después de colocar las ramas, pueden aparecer sustituyentes complejos anidados, por ejemplo `[3,3-dimetil-1-(1-metiletil)butil]`.

## Haloalcanos: reglas aplicadas

| Aspecto | Regla | Ejemplo |
|---|---|---|
| Los halógenos solo se expresan como prefijos (fluoro, cloro, bromo, yodo), nunca como sufijo | P‑61.3.1 | `2-bromopropano` |
| La cadena principal se busca solo entre carbonos (la más larga); los halógenos cuentan en el criterio de mayor número de sustituyentes | P‑44.3, P‑45.2.1 | `1-bromo-2,2-dimetilpropano` |
| Localizadores más bajos para todos los prefijos en conjunto (halógenos y alquilos) y, si hay empate, para el citado primero en orden alfabético | P‑31.1.4, P‑14.5 | `1-bromo-4-cloro-2-metilpentano` |
| Orden alfabético por idioma (en español *yodo* va al final; en inglés *iodo* va entre *ethyl* y *methyl*) | P‑14.5 | `2-cloro-3-yodobutano` |
| Halógenos dentro de un sustituyente: prefijo complejo entre paréntesis | P‑16.5 | `3-(clorometil)-4-etilhexano`, `(2-bromoetil)` |
| Omisión de localizadores: metano, etano monosustituido y sustitución total por un mismo prefijo | P‑14.3.4 | `triclorometano`, `cloroetano`, `hexafluoroetano` |
| Los prefijos retenidos (terc-butil, isopropil…) solo se usan sin sustituir | P‑29.6; 1979 A‑2.25 | `(1-cloro-1-metiletil)`, no «cloroisopropil» |
| Nombre radicofuncional R–X: aceptado en nomenclatura general, no es nombre preferido | P‑61.3.2 | `cloruro de sec-butilo` |

**Sobre sec- y neo-.** Solo existen como parte de prefijos retenidos **no sustituidos**: *sec*-butil (butan-2-il) y neopentil (2,2-dimetilpropil), igual que isopropil, isobutil, isopentil, *terc*-butil y *terc*-pentil (IUPAC 1979, A‑2.25). El programa no usa «sec-pentil» ni «neohexil», que la IUPAC nunca aceptó porque son ambiguos. En el orden alfabético, *sec*- y *terc*- (en cursiva) se ignoran, mientras que iso- y neo- sí cuentan. Por eso *sec*-butil se ordena en la «b» y neopentil en la «n».

## Formas de nomenclatura presentadas

| Estilo | Sustituyente de ejemplo | Fundamento |
|---|---|---|
| **Preferido 2013 (PIN)** | `propan-2-il`, `butan-2-il`, `2-metilpropil`, `3-metilbutan-2-il`, `terc-butil` | P‑29.2 y P‑29.6: cadena más larga que contiene el átomo de unión, con la valencia libre en el localizador más bajo |
| **Retenido / tradicional** | `isopropil`, `sec-butil`, `isobutil`, `terc-butil`, `isopentil`, `neopentil`, `terc-pentil` | IUPAC 1979, A‑2.25 (prefijos retenidos solo sin sustituir); los demás se nombran como en 1979 |
| **Sistemático 1979** | `1-metiletil`, `1-metilpropil`, `1,1-dimetiletil`, `1,2-dimetilpropil` | IUPAC 1979, A‑2.6: cadena más larga que empieza en el átomo de unión (C1) |

El orden alfanumérico y, por lo tanto, la numeración cuando hay empate se recalculan para cada estilo. Por eso el mismo compuesto puede llamarse `1-metil-4-(propan-2-il)ciclohexano` (PIN) y `1-isopropil-4-metilciclohexano` (tradicional).

## Estructura

```
alcanos_iupac/
├── index.php            Interfaz (emite el token CSRF y la sesión)
├── api.php              Servicio JSON: generar / analizar / desplegar / verificar
├── lib/IupacAlcanos.php Motor: lector SMILES, generador y algoritmo de nomenclatura (sin dependencias)
├── assets/app.js        Lógica del cliente e integración con JSME
├── assets/estilos.css   Estilos adaptables (móvil y escritorio)
├── jsme/                JSME 2024.04.29 (BSD‑3), servido localmente, sin CDN
└── tests/               Scripts de validación (acceso web bloqueado)
```

## Instalación en el VPS (Hostinger)

Requisitos: PHP ≥ 7.4 (probado en PHP 8.4) con `mbstring` y sesiones habilitadas. No usa base de datos.

1. Suba la carpeta `alcanos_iupac/` al directorio público, por ejemplo `/var/www/html/alcanos_iupac/`. Si actualiza desde la v1.0, reemplace todos los archivos.
2. **Apache**: los `.htaccess` incluidos bloquean `lib/` y `tests/` (requiere `AllowOverride All`).
   **Nginx**: agregue al `server {}`:
   ```nginx
   location ~ ^/alcanos_iupac/(lib|tests)/ { deny all; return 404; }
   ```
3. Use HTTPS. La cookie de sesión se marca como `Secure`, `HttpOnly` y `SameSite=Strict`.
4. Abra `https://su-dominio/alcanos_iupac/`.

Prueba local: `php -S 127.0.0.1:8080` dentro de la carpeta.

## Reglas implementadas (IUPAC 2013)

- **P‑44.3**: la cadena principal es la cadena continua más larga (hasta 20 C, icosano).
- **P‑44.1.2.2**: en los cicloalcanos, el anillo es preferido a la cadena como hidruro progenitor. Si el estudiante dibuja una cadena más larga que el anillo, la justificación advierte que las reglas de 1979/1993 tomarían la cadena como progenitor.
- **P‑45.2.1 → P‑45.2.3 → P‑45.5**: si hay varias cadenas de igual longitud, se prefiere la que tiene más sustituyentes; luego la de localizadores más bajos; luego la que da el localizador más bajo al prefijo citado primero; por último, la que produce el nombre primero en orden alfanumérico.
- **P‑14.4 / P‑31.1.4**: localizadores más bajos, comparados en el primer punto de diferencia. En los sustituyentes, la valencia libre recibe el localizador más bajo.
- **P‑14.5**: orden alfanumérico. No se consideran di/tri ni sec‑/terc‑; sí se considera iso‑ y los multiplicadores dentro de un prefijo complejo.
- **P‑16.3**: di…icosa para prefijos simples, con paréntesis cuando llevan localizador, por ejemplo `di(propan-2-il)`; bis/tris/tetrakis… para prefijos complejos, por ejemplo `bis(2-metilpropil)`; `di-terc-butil`.
- **P‑16.5.4**: signos de inclusión anidados en el orden ( ) → [ ] → { }.
- **P‑29 / P‑29.6**: prefijos sustituyentes preferidos; los complejos se nombran de forma recursiva.
- **P‑14.3.4**: en un cicloalcano monosustituido se omite el localizador.

## Validación realizada

Con **OPSIN** (convierte nombre → estructura) y **RDKit** (compara el SMILES canónico), sobre los tres estilos en inglés:

| Prueba | Estructuras | Nombres | Discrepancias |
|---|---|---|---|
| `validar.py`: estructuras del generador, los 3 niveles (la mitad con halógenos) | 6000 | 18 000 | **0** |
| `validar.py`: invariancia ante renumeración aleatoria de átomos | 3238 | 12 952 | **0** |
| `arboles_aleatorios.py`: estructuras arbitrarias como las que puede dibujar un estudiante (árboles de 1–26 C y monociclos de 3–10 C, el 60 % con 1–5 halógenos) | 8222 | 24 666 | **0** |
| Nombres radicofuncionales en inglés (PIN, retenido y sistemático) | 800 | 2400 | **0** |

```bash
pip install rdkit pyopsin
export OPSIN_JAR=$(python3 -c "import pyopsin,os;print(os.path.join(os.path.dirname(pyopsin.__file__),'opsin_cli.jar'))")
cd tests
python3 validar.py 2000
python3 arboles_aleatorios.py 2000
```

OPSIN comprueba que cada nombre corresponde exactamente a la estructura. La elección de la cadena y de la numeración preferidas la garantiza la búsqueda exhaustiva del motor, que evalúa todas las cadenas candidatas y todas las numeraciones con los criterios anteriores.

## Alcance y limitaciones

- Solo alcanos y haloalcanos saturados acíclicos (≤ 20 C en la cadena principal, ≤ 60 C en total) y cicloalcanos **monocíclicos** con sustituyentes alquilo o halógeno. No incluye estereoquímica (cis/trans, R/S): si se dibuja, se ignora. Tampoco incluye anillos múltiples ni heteroátomos distintos de los halógenos. No se usan nombres triviales como cloroformo o tetracloruro de carbono.
- El estilo tradicional usa los criterios de 2013 para elegir la cadena principal. En casos muy poco frecuentes, los criterios adicionales de 1979 (C‑13.11: ramas más pequeñas con más carbonos, ramas menos ramificadas) podrían elegir otra cadena de igual longitud.
- Los nombres en español siguen la traducción usual de los prefijos (metil, etil, terc‑butil, propan‑2‑il). El orden alfabético se calcula en cada idioma.

## Referencias

- Favre, H. A.; Powell, W. H. *Nomenclature of Organic Chemistry. IUPAC Recommendations and Preferred Names 2013*. Cambridge: RSC, 2014. doi:10.1039/9781849733069
- IUPAC. *Nomenclature of Organic Chemistry, Sections A, B, C, D, E, F and H* (1979 ed.). Oxford: Pergamon, 1979 (reglas A‑2.6, A‑2.25).
- IUPAC. *A Guide to IUPAC Nomenclature of Organic Compounds (Recommendations 1993)*. Oxford: Blackwell, 1993.
- Bienfait, B.; Ertl, P. JSME: a free molecule editor in JavaScript. *J. Cheminform.* 2013, 5, 24. doi:10.1186/1758-2946-5-24
- Lowe, D. M. et al. Chemical name to structure: OPSIN. *J. Chem. Inf. Model.* 2011, 51, 739–753. doi:10.1021/ci100384d
