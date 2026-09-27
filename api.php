<?php
/**
 * api.php – Servicio JSON del generador de alcanos.
 *
 *  POST accion=generar    {tipo, nivel}        → {id, smiles, formula, carbonos, longitud}
 *  POST accion=analizar   {smiles}             → igual que generar, para una estructura dibujada en JSME
 *  POST accion=desplegar  {id}                 → nombres IUPAC + explicación + SMILES numerado
 *  POST accion=verificar  {id, respuesta}      → {correcto, mensaje}
 *
 * El nombre correcto se guarda SOLO en la sesión del servidor y no viaja al
 * navegador hasta que el estudiante pulsa <Desplegar>.
 */
declare(strict_types=1);

require __DIR__ . '/lib/IupacAlcanos.php';

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
    $s = str_replace(['‐', '‑', '‒', '–', '—', '−'], '-', $s);
    $s = preg_replace('/\s+/u', '', $s);
    $s = str_replace(['ter-', 'tert-', 't-butil'], ['terc-', 'terc-', 'terc-butil'], $s);
    // variantes frecuentes de los prefijos sec-, neo-, iso- y del yodo
    $s = preg_replace('/(^|[\-\d,(\[{])s-but/u', '$1sec-but', $s);   // s-butil → sec-butil
    $s = str_replace(['sec.-', 'secbut', 'neo-', 'iso-'], ['sec-', 'sec-but', 'neo', 'iso'], $s);
    $s = str_replace('iodo', 'yodo', $s);                              // iodo/yodo (inglés/español)
    $s = str_replace('ioduro', 'yoduro', $s);
    return $s;
}

/**
 * Calcula todos los nombres de una molécula, guarda el ejercicio en la sesión
 * y devuelve los datos públicos (sin el nombre).
 */
function crearEjercicio(Molecula $mol, bool $propio): array
{
    $nom = [];
    $radico = [];
    foreach (['es', 'en'] as $l) {
        foreach (['pin', 'trad', 'sis79'] as $e) {
            $nom[$l][$e] = (new IupacAlcanos($mol, $l, $e))->nombrar();
            if ($nom[$l][$e]['radicofuncional']) { $radico[$l][$e] = $nom[$l][$e]['radicofuncional']; }
        }
    }
    $pin = $nom['es']['pin'];

    // Equivalencias de los sustituyentes: preferido 2013 / retenido / sistemático 1979
    $trad = new IupacAlcanos($mol, 'es', 'trad');
    $sis  = new IupacAlcanos($mol, 'es', 'sis79');
    $eq = [];
    foreach ($pin['subs'] as $s) {
        $p = $s['info']['nombre'];
        if (isset($eq[$p]) || !empty($s['info']['hal'])) { continue; }
        $t = $trad->sustituyente($s['atomo'], $s['desde'])['nombre'];
        $y = $sis->sustituyente($s['atomo'], $s['desde'])['nombre'];
        if ($p !== $t || $p !== $y) { $eq[$p] = ['pin' => $p, 'trad' => $t, 'sis' => $y]; }
    }

    $mapas = [];
    foreach ($pin['cadena'] as $i => $a) { $mapas[$a] = $i + 1; }
    $inicio = $pin['cadena'][0];

    $id = bin2hex(random_bytes(8));
    $_SESSION['ej'][$id] = [
        'nombres'      => [
            'es' => array_map(fn($r) => $r['nombre'], $nom['es']),
            'en' => array_map(fn($r) => $r['nombre'], $nom['en']),
        ],
        'equivalencias' => array_values($eq),
        'radicofuncional' => $radico,
        'explicacion'  => $pin['explicacion'],
        'smiles_num'   => IupacAlcanos::smiles($mol, $mapas, $inicio),
        'propio'       => $propio,
    ];
    if (count($_SESSION['ej']) > 30) {
        $_SESSION['ej'] = array_slice($_SESSION['ej'], -30, null, true);
    }
    return [
        'id'       => $id,
        'smiles'   => IupacAlcanos::smiles($mol, [], $inicio),
        'formula'  => IupacAlcanos::formula($mol),
        'carbonos' => $mol->numCarbonos(),
        'longitud' => count($pin['cadena']),
        'ciclico'  => $pin['anillo'],
        'nombre_pin' => $pin['nombre'], // uso interno (no se envía)
    ];
}

function obtenerEjercicio(array $entrada): array
{
    $id = (string)($entrada['id'] ?? '');
    $ej = $_SESSION['ej'][$id] ?? null;
    if (!$ej) { responder(['error' => 'Ejercicio no encontrado. Genere o dibuje una nueva estructura.'], 404); }
    return $ej;
}

/** Devuelve $b solo si difiere de todos los nombres ya listados. */
function distinto(?string $b, array $previos): ?string
{
    return ($b === null || in_array($b, $previos, true)) ? null : $b;
}

try {
    switch ($accion) {

        case 'generar': {
            $tipos = ['aleatorio', 'lineal', 'ramificado', 'ciclico', 'ciclico_ramificado'];
            $tipo  = in_array($entrada['tipo'] ?? '', $tipos, true) ? $entrada['tipo'] : 'aleatorio';
            $nivel = max(1, min(3, (int)($entrada['nivel'] ?? 2)));
            $conHal = !empty($entrada['halogenos']);
            $previo = $_SESSION['ultimo_nombre'] ?? '';
            $datos = null;
            for ($intentos = 0; $intentos < 15; $intentos++) {
                try {
                    $mol = GeneradorAlcanos::generar($tipo, $nivel, $conHal);
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
                $mol = Molecula::desdeSmiles($smi);
                $datos = crearEjercicio($mol, true);
            } catch (NomenclaturaException $e) {
                responder(['error' => $e->getMessage(), 'dominio' => true], 422);
            }
            unset($datos['nombre_pin']);
            responder($datos);
        }

        case 'desplegar': {
            $ej = obtenerEjercicio($entrada);
            $es = $ej['nombres']['es'];
            $en = $ej['nombres']['en'];
            responder([
                'nombre'        => $es['pin'],
                'nombre_trad'   => distinto($es['trad'], [$es['pin']]),
                'nombre_sis'    => distinto($es['sis79'], [$es['pin'], $es['trad']]),
                'nombre_en'     => $en['pin'],
                'nombre_en_trad'=> distinto($en['trad'], [$en['pin']]),
                'nombre_en_sis' => distinto($en['sis79'], [$en['pin'], $en['trad']]),
                'equivalencias' => $ej['equivalencias'],
                'radicofuncional' => $ej['radicofuncional']['es']['trad'] ?? null,
                'radicofuncional_pin' => distinto($ej['radicofuncional']['es']['pin'] ?? null, [$ej['radicofuncional']['es']['trad'] ?? '']),
                'radicofuncional_en' => $ej['radicofuncional']['en']['trad'] ?? null,
                'explicacion'   => $ej['explicacion'],
                'smiles_num'    => $ej['smiles_num'],
            ]);
        }

        case 'verificar': {
            $ej = obtenerEjercicio($entrada);
            $resp = normalizar(mb_substr((string)($entrada['respuesta'] ?? ''), 0, 250));
            if ($resp === '') { responder(['correcto' => false, 'mensaje' => 'Escriba un nombre antes de verificar.']); }

            $etiquetas = [
                'pin'   => 'Coincide con el nombre IUPAC preferido (2013).',
                'trad'  => 'Nombre aceptado con prefijos tradicionales retenidos (isopropil, sec-butil, isobutil, terc-butil…).',
                'sis79' => 'Nombre aceptado con prefijos sistemáticos numerados desde el punto de unión (IUPAC 1979).',
            ];
            foreach (['pin', 'trad', 'sis79'] as $e) {
                foreach (['es', 'en'] as $l) {
                    if ($resp === normalizar($ej['nombres'][$l][$e])) {
                        responder(['correcto' => true, 'mensaje' => '¡Correcto! ' . $etiquetas[$e]]);
                    }
                }
            }
            foreach ($ej['radicofuncional'] as $porEstilo) {
                foreach ($porEstilo as $n) {
                    if ($resp === normalizar($n)) {
                        responder(['correcto' => true, 'mensaje' => '¡Correcto! Nombre radicofuncional (de clase funcional), aceptado en nomenclatura general; el nombre IUPAC preferido es el sustitutivo con el prefijo halo (P-61.3).']);
                    }
                }
            }
            $validos = [];
            foreach ($ej['nombres'] as $porEstilo) { foreach ($porEstilo as $n) { $validos[] = normalizar($n); } }
            foreach ($ej['radicofuncional'] as $porEstilo) { foreach ($porEstilo as $n) { $validos[] = normalizar($n); } }
            $sinSignos = fn(string $s) => preg_replace('/[\s,\-\(\)\[\]]/u', '', $s);
            foreach ($validos as $v) {
                if ($sinSignos($v) === $sinSignos($resp)) {
                    responder(['correcto' => false, 'nivel' => 'casi', 'mensaje' => 'Casi: los componentes son correctos, pero revise comas, guiones y paréntesis.']);
                }
            }
            $sinLocs = fn(string $s) => preg_replace('/[\d,\-]/', '', $s);
            foreach ($validos as $v) {
                if ($sinLocs($v) === $sinLocs($resp)) {
                    responder(['correcto' => false, 'nivel' => 'casi', 'mensaje' => 'Sustituyentes y cadena principal correctos; revise los localizadores (numeración).']);
                }
            }
            responder(['correcto' => false, 'nivel' => 'mal', 'mensaje' => 'Incorrecto. Revise la cadena principal, los sustituyentes y el orden alfabético.']);
        }

        default:
            responder(['error' => 'Acción no válida'], 400);
    }
} catch (Throwable $e) {
    error_log('[alcanos_iupac] ' . $e->getMessage());
    responder(['error' => 'Error interno al procesar la estructura.'], 500);
}
