<?php
// Genera N estructuras y vuelca: SMILES, nombres en inglés (pin, trad, sis79) y en español (pin, trad, sis79)
require __DIR__.'/../lib/IupacAlcanos.php';
$N = (int)($argv[1] ?? 2000);
$tipos=['lineal','ramificado','ciclico','ciclico_ramificado','ramificado','ciclico_ramificado'];
for($i=0;$i<$N;$i++){
  $t=$tipos[$i%6]; $nv=1+intdiv($i,6)%3;
  $m=GeneradorAlcanos::generar($t,$nv,$i%2===1);
  $out=[IupacAlcanos::smiles($m)];
  foreach(['en','es'] as $l) foreach(['pin','trad','sis79'] as $e){ $out[]=(new IupacAlcanos($m,$l,$e))->nombrar()['nombre']; }
  echo implode("\t",$out),"\n";
}
