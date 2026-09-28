<?php
// Lee SMILES (uno por línea) y escribe: en-pin, es-pin, en-trad, en-sis79, [en-clase funcional]
require __DIR__ . '/../lib/IupacOrganica.php';
while (($l = fgets(STDIN)) !== false) {
    $l = trim($l); if ($l === '') { continue; }
    try {
        $m = MolOrg::desdeSmiles($l); $o = [];
        foreach ([['en', 'pin'], ['es', 'pin'], ['en', 'trad'], ['en', 'sis79']] as [$a, $b]) { $o[] = (new IupacOrganica($m, $a, $b))->nombrar()['nombre']; }
        $o[] = (new IupacOrganica($m, 'en', 'trad'))->nombreClaseFuncional() ?? '';
        echo implode("\t", $o), "\n";
    } catch (NomenclaturaException $e) { echo "ERR\t", $e->getMessage(), "\n"; }
}
