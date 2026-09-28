<?php
// Regresión: en el dominio de la v3.0 (alcanos, cicloalcanos, haloalcanos) el motor nuevo
// debe producir exactamente los mismos nombres que IupacAlcanos (6 variantes + radicofuncional).
require __DIR__ . '/../lib/IupacAlcanos.php';
require __DIR__ . '/../lib/IupacOrganica.php';
$N = (int)($argv[1] ?? 3000);
$tipos = ['lineal', 'ramificado', 'ciclico', 'ciclico_ramificado'];
$dif = 0; $tot = 0;
for ($i = 0; $i < $N; $i++) {
    $m = GeneradorAlcanos::generar($tipos[$i % 4], 1 + intdiv($i, 4) % 3, $i % 2 === 1);
    $smi = IupacAlcanos::smiles($m);
    $mo = MolOrg::desdeSmiles($smi);
    foreach (['es', 'en'] as $l) foreach (['pin', 'trad', 'sis79'] as $e) {
        $v = (new IupacAlcanos($m, $l, $e))->nombrar();
        $n = (new IupacOrganica($mo, $l, $e))->nombrar();
        $tot++;
        // corrección intencional de la v4.0: sin localizadores no se escribe guion entre prefijos
        // (v3.0: «bromo-clorometil»; v4.0: «bromoclorometil»)
        $sinGuion = fn($x) => preg_replace('/(?<=[a-z])-(?=[a-z])/', '', $x);
        if ($v['nombre'] !== $n['nombre'] && $sinGuion($v['nombre']) === $n['nombre']) { $orto = ($orto ?? 0) + 1; continue; }
        if ($v['nombre'] !== $n['nombre']) { $dif++; if ($dif < 25) echo "DIF $smi [$l/$e] v3: {$v['nombre']}  v4: {$n['nombre']}\n"; }
        if ($e === 'trad' && ($v['radicofuncional'] ?? null) !== ($n['radicofuncional'] ?? null)) {
            $dif++; if ($dif < 25) echo "DIF-RADICO $smi [$l] v3: " . ($v['radicofuncional'] ?? '-') . "  v4: " . ($n['radicofuncional'] ?? '-') . "\n";
        }
    }
}
echo "comparaciones $tot  diferencias $dif  correcciones ortográficas (guion entre prefijos sin localizador) " . ($orto ?? 0) . "\n";
