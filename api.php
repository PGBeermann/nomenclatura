<?php
/**
 * api.php – Servicio JSON del ejercitador de nomenclatura orgánica (v4.0).
 *
 *  POST accion=generar    {tipo, nivel, familias[]} → {id, smiles, formula, carbonos, longitud}
 *  POST accion=analizar   {smiles}                  → igual que generar, para una estructura dibujada en JSME
 *  POST accion=desplegar  {id}                      → nombres IUPAC + alternativas + explicación + SMILES numerado
 *  POST accion=verificar  {id, respuesta}           → {correcto, mensaje}
 *
 * El nombre correcto se guarda SOLO en la sesión del servidor y no viaja al
 * navegador hasta que el estudiante pulsa <Desplegar>.
 */
declare(strict_types=1);

require __DIR__ . '/lib/IupacOrganica.php';

session_set_cookie_params([
    'lifetime' => 0,
    'path'     => '/',
    'secure'   => !empty($_SERVER['HTTPS']) && $_SERVER['HTTPS'] !== 'off',
    'httponly' => true,
    'samesite' => 'Strict',
]);
session_start();

header('Content-Type: application/json; charset=utf-8');
header('Cache-Control: no-store');
header('X-Content-Type-Options: nosniff');

function responder(array $data, int $codigo = 200): void
{
    http_response_code($codigo);
    echo json_encode($data, JSON_UNESCAPED_UNICODE | JSON_UNESCAPED_SLASHES);
    exit;
}

if ($_SERVER['REQUEST_METHOD'] !== 'POST') {
    responder(['error' => 'Método no permitido'], 405);
}

$entrada = json_decode(file_get_contents('php://input') ?: '[]', true);
if (!is_array($entrada)) { $entrada = []; }

// Protección CSRF: el token se emite en index.php
$token = $_SERVER['HTTP_X_CSRF_TOKEN'] ?? '';
if (empty($_SESSION['csrf']) || !hash_equals($_SESSION['csrf'], $token)) {
    responder(['error' => 'Sesión caducada. Recargue la página.'], 403);
}

$accion = (string)($entrada['accion'] ?? '');

/** Normaliza un nombre para comparar la respuesta del estudiante. */
function normalizar(string $s): string
{
    $s = mb_strtolower(trim($s), 'UTF-8');
    $s = strtr($s, ['á' => 'a', 'é' => 'e', 'í' => 'i', 'ó' => 'o', 'ú' => 'u', 'ü' => 'u', '’' => "'", '′' => "'"]);
    $s = str_replace(['‐', '‑', '‒', '–', '—', '−'], '-', $s);
    $s = preg_replace('/\s+/u', ' ', $s);
    $s = str_replace(['ter-', 'tert-', 't-butil'], ['terc-', 'terc-', 'terc-butil'], $s);
    // variantes frecuentes de los prefijos sec-, neo-, iso- y del yodo
    $s = preg_replace('/(^|[\-\d,(\[{ ])s-but/u', '$1sec-but', $s);   // s-butil → sec-butil
    $s = str_replace(['sec.-', 'secbut', 'neo-', 'iso-'], ['sec-', 'sec-but', 'neo', 'iso'], $s);
    $s = str_replace(['iodo', 'ioduro'], ['yodo', 'yoduro'], $s);
    return str_replace(' ', '', $s);
}

const ETIQUETAS = [
    'pin'   => 'Nombre IUPAC preferido (2013)',
    'trad'  => 'Con prefijos tradicionales retenidos (isopropil, sec-butil, vinil, alil…) y criterios de 1979/1993',
    'sis79' => 'Con prefijos sistemáticos numerados desde el punto de unión (IUPAC 1979)',
    'clase' => 'Nombre de clase funcional o radicofuncional (aceptado en nomenclatura general, no preferido)',
    'comun' => 'Nombre común o retenido aceptado en nomenclatura general (no preferido)',
];

/**
 * Calcula todos los nombres de una molécula, guarda el ejercicio en la sesión
 * y devuelve los datos públicos (sin el nombre).
 */
function crearEjercicio(MolOrg $mol, bool $propio): array
{
    $nom = []; $res = [];
    foreach (['es', 'en'] as $l) {
        foreach (['pin', 'trad', 'sis79'] as $e) {
            $res[$l][$e] = (new IupacOrganica($mol, $l, $e))->nombrar();
            $nom[$l][$e] = $res[$l][$e]['nombre'];
        }
    }
    $pin = $res['es']['pin'];

    // nombres de clase funcional (con prefijos tradicionales y, si difiere, con los preferidos)
    $clase = [];
    foreach (['es', 'en'] as $l) {
        $c = array_values(array_unique(array_filter([$res[$l]['trad']['radicofuncional'], $res[$l]['pin']['radicofuncional']])));
        $clase[$l] = $c;
    }
    // nombres comunes, indexados por el PIN en inglés
    $com = IupacOrganica::COMUNES[$nom['en']['pin']] ?? [[], []];
    $comunes = ['es' => $com[0], 'en' => $com[1]];

    // Equivalencias de los sustituyentes: preferido 2013 / retenido / sistemático 1979
    $trad = new IupacOrganica($mol, 'es', 'trad');
    $sis  = new IupacOrganica($mol, 'es', 'sis79');
    $eq = [];
    foreach ($pin['subs'] as $s) {
        if (($s['atomo'] ?? -1) < 0 || !isset($s['desde'])) { continue; }
        $p = $s['info']['nombre'];
        if (isset($eq[$p])) { continue; }
        $t = $trad->sustituyente($s['atomo'], $s['desde'], $s['orden'] ?? 1)['nombre'];
        $y = $sis->sustituyente($s['atomo'], $s['desde'], $s['orden'] ?? 1)['nombre'];
        if ($p !== $t || $p !== $y) { $eq[$p] = ['pin' => $p, 'trad' => $t, 'sis' => $y]; }
    }

    $explicacion = $pin['explicacion'];
    // Diferencia 2013 frente a 1979/1993 en la elección de la cadena principal
    $cp = $pin['cadena']; $ct = $res['es']['trad']['cadena']; sort($cp); sort($ct);
    if (!$pin['anillo'] && $cp !== $ct) {
        $explicacion[] = 'Nota: con los criterios de 1979/1993, que anteponían el mayor número de enlaces múltiples a la longitud de la cadena, '
            . 'se elegiría otra cadena principal: «' . $nom['es']['trad'] . '». Las recomendaciones de 2013 eligen primero la cadena más larga (P-44.3) '
            . 'y expresan el doble enlace externo con un prefijo «-ilideno».';
    }

    $mapas = [];
    foreach ($pin['cadena'] as $i => $a) { $mapas[$a] = $i + 1; }
    $inicio = $pin['cadena'][0];

    $id = bin2hex(random_bytes(8));
    $_SESSION['ej'][$id] = [
        'nombres'       => $nom,
        'clase'         => $clase,
        'comunes'       => $comunes,
        'equivalencias' => array_values($eq),
        'explicacion'   => $explicacion,
        'smiles_num'    => $mol->smiles($mapas, $inicio),
        'propio'        => $propio,
    ];
    if (count($_SESSION['ej']) > 30) {
        $_SESSION['ej'] = array_slice($_SESSION['ej'], -30, null, true);
    }
    return [
        'id'       => $id,
        'smiles'   => $mol->smiles([], $inicio),
        'formula'  => $mol->formula(),
        'carbonos' => $mol->numCarbonos(),
        'longitud' => max(count($pin['cadena']), 4),
        'ciclico'  => $pin['anillo'],
        'nombre_pin' => $nom['es']['pin'], // uso interno (no se envía)
    ];
}

function obtenerEjercicio(array $entrada): array
{
    $id = (string)($entrada['id'] ?? '');
    $ej = $_SESSION['ej'][$id] ?? null;
    if (!$ej) { responder(['error' => 'Ejercicio no encontrado. Genere o dibuje una nueva estructura.'], 404); }
    return $ej;
}

/** Lista [etiqueta, nombre] de formas alternativas distintas del PIN, en un idioma. */
function alternativas(array $ej, string $l): array
{
    $vistos = [$ej['nombres'][$l]['pin'] => true];
    $out = [];
    foreach (['trad', 'sis79'] as $e) {
        $n = $ej['nombres'][$l][$e];
        if (!isset($vistos[$n])) { $vistos[$n] = true; $out[] = ['etiqueta' => ETIQUETAS[$e], 'nombre' => $n]; }
    }
    foreach ($ej['clase'][$l] as $n) {
        if (!isset($vistos[$n])) { $vistos[$n] = true; $out[] = ['etiqueta' => ETIQUETAS['clase'], 'nombre' => $n]; }
    }
    foreach ($ej['comunes'][$l] as $n) {
        if (!isset($vistos[$n])) { $vistos[$n] = true; $out[] = ['etiqueta' => ETIQUETAS['comun'], 'nombre' => $n]; }
    }
    return $out;
}

try {
    switch ($accion) {

        case 'generar': {
            $tipos = ['aleatorio', 'lineal', 'ramificado', 'ciclico', 'ciclico_ramificado'];
            $tipo  = in_array($entrada['tipo'] ?? '', $tipos, true) ? $entrada['tipo'] : 'aleatorio';
            $nivel = max(1, min(3, (int)($entrada['nivel'] ?? 2)));
            $familias = is_array($entrada['familias'] ?? null) ? array_map('strval', $entrada['familias']) : [];
            if (!empty($entrada['halogenos'])) { $familias[] = 'halogenos'; }        // compatibilidad v3.0
            $familias = array_values(array_unique(array_intersect($familias, GeneradorOrganico::FAMILIAS)));
            $previo = $_SESSION['ultimo_nombre'] ?? '';
            $datos = null;
            for ($intentos = 0; $intentos < 40; $intentos++) {
                try {
                    $mol = GeneradorOrganico::generar($tipo, $nivel, $familias);
                    $d = crearEjercicio($mol, false);
                } catch (NomenclaturaException $e) {
                    continue; // estructura fuera de dominio: se genera otra
                }
                $datos = $d;
                if ($d['nombre_pin'] !== $previo) { break; }
            }
            if (!$datos) { throw new RuntimeException('No se pudo generar la estructura'); }
            $_SESSION['ultimo_nombre'] = $datos['nombre_pin'];
            unset($datos['nombre_pin']);
            responder($datos);
        }

        case 'analizar': {
            // Estructura dibujada por el estudiante en JSME (SMILES)
            $smi = (string)($entrada['smiles'] ?? '');
            try {
                $mol = MolOrg::desdeSmiles($smi);
                $datos = crearEjercicio($mol, true);
            } catch (NomenclaturaException $e) {
                responder(['error' => $e->getMessage(), 'dominio' => true], 422);
            }
            unset($datos['nombre_pin']);
            responder($datos);
        }

        case 'desplegar': {
            $ej = obtenerEjercicio($entrada);
            responder([
                'nombre'          => $ej['nombres']['es']['pin'],
                'alternativas'    => alternativas($ej, 'es'),
                'nombre_en'       => $ej['nombres']['en']['pin'],
                'alternativas_en' => array_map(fn($a) => $a['nombre'], alternativas($ej, 'en')),
                'equivalencias'   => $ej['equivalencias'],
                'explicacion'     => $ej['explicacion'],
                'smiles_num'      => $ej['smiles_num'],
            ]);
        }

        case 'verificar': {
            $ej = obtenerEjercicio($entrada);
            $resp = normalizar(mb_substr((string)($entrada['respuesta'] ?? ''), 0, 300));
            if ($resp === '') { responder(['correcto' => false, 'mensaje' => 'Escriba un nombre antes de verificar.']); }

            $mensajes = [
                'pin'   => '¡Correcto! Coincide con el nombre IUPAC preferido (2013).',
                'trad'  => '¡Correcto! Nombre aceptado con prefijos tradicionales (isopropil, sec-butil, vinil, alil…) o con los criterios de 1979/1993.',
                'sis79' => '¡Correcto! Nombre aceptado con prefijos sistemáticos numerados desde el punto de unión (IUPAC 1979).',
                'clase' => '¡Correcto! Es un nombre de clase funcional (radicofuncional), aceptado en nomenclatura general; el nombre IUPAC preferido es el sustitutivo.',
                'comun' => '¡Correcto! Es un nombre común o retenido aceptado en nomenclatura general, pero no es el nombre IUPAC preferido. Pulse «Desplegar» para verlo.',
            ];
            $validos = [];                                   // normalizado => clave de mensaje
            foreach (['pin', 'trad', 'sis79'] as $e) {
                foreach (['es', 'en'] as $l) { $validos[normalizar($ej['nombres'][$l][$e])] ??= $e; }
            }
            foreach (['es', 'en'] as $l) {
                foreach ($ej['clase'][$l] as $n) {
                    $validos[normalizar($n)] ??= 'clase';
                    $validos[normalizar(str_replace(['(', ')'], '', $n))] ??= 'clase';   // etilmetilamina
                }
                foreach ($ej['comunes'][$l] as $n) { $validos[normalizar($n)] ??= 'comun'; }
            }
            if (isset($validos[$resp])) {
                responder(['correcto' => true, 'mensaje' => $mensajes[$validos[$resp]]]);
            }
            $sinSignos = fn(string $s) => preg_replace('/[\s,\-\(\)\[\]\{\}\']/u', '', $s);
            foreach ($validos as $v => $k) {
                if ($sinSignos($v) === $sinSignos($resp)) {
                    responder(['correcto' => false, 'nivel' => 'casi', 'mensaje' => 'Casi: los componentes son correctos, pero revise comas, guiones y paréntesis.']);
                }
            }
            $sinLocs = fn(string $s) => preg_replace('/(?<=^|[\-,])n(?=[\d,\-\'])|[\d,\-\']/u', '', $s);
            $digitos = function (string $s): string { preg_match_all('/\d+/', $s, $m); $d = $m[0]; sort($d); return implode(',', $d); };
            foreach ($validos as $v => $k) {
                if ($sinLocs($v) === $sinLocs($resp)) {
                    if ($digitos($v) === $digitos($resp)) {
                        responder(['correcto' => false, 'nivel' => 'casi', 'mensaje' => 'Los localizadores son correctos, pero están mal ubicados. Desde 1993 el localizador del sufijo o del enlace múltiple se escribe justo antes de él: «butan-2-ol», no «2-butanol»; «but-1-eno», no «1-buteno».']);
                    }
                    responder(['correcto' => false, 'nivel' => 'casi', 'mensaje' => 'Sustituyentes, grupo principal y cadena correctos; revise los localizadores (numeración).']);
                }
            }
            responder(['correcto' => false, 'nivel' => 'mal', 'mensaje' => 'Incorrecto. Revise el grupo principal (sufijo), la cadena o anillo principal, los prefijos y el orden alfabético.']);
        }

        default:
            responder(['error' => 'Acción no válida'], 400);
    }
} catch (Throwable $e) {
    error_log('[nomenclatura] ' . $e->getMessage());
    responder(['error' => 'Error interno al procesar la estructura.'], 500);
}
