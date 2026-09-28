<?php
// Genera N ejercicios con familias aleatorias y vuelca: SMILES \t en-pin \t en-trad \t en-sis79 \t en-clase-funcional(trad) \t es-pin
require __DIR__ . '/../lib/IupacOrganica.php';
$N = (int)($argv[1] ?? 500); mt_srand((int)($argv[2] ?? 1));
$tipos = ['lineal', 'ramificado', 'ciclico', 'ciclico_ramificado', 'aleatorio'];
$fam = GeneradorOrganico::FAMILIAS; $err = 0;
for ($i = 0; $i < $N; $i++) {
    $sel = []; foreach ($fam as $f) { if (random_int(0, 3) === 0) { $sel[] = $f; } }
    $m = GeneradorOrganico::generar($tipos[$i % 5], 1 + intdiv($i, 5) % 3, $sel);
    try {
        $o = [$m->smiles()];
        foreach (['pin', 'trad', 'sis79'] as $e) { $o[] = (new IupacOrganica($m, 'en', $e))->nombrar()['nombre']; }
        $o[] = (new IupacOrganica($m, 'en', 'trad'))->nombreClaseFuncional() ?? '';
        $o[] = (new IupacOrganica($m, 'es', 'pin'))->nombrar()['nombre'];
        echo implode("\t", $o), "\n";
    } catch (NomenclaturaException $e) { $err++; fwrite(STDERR, "ERR\t" . $m->smiles() . "\t" . $e->getMessage() . "\n"); }
}
fwrite(STDERR, "fuera de dominio: $err\n");
