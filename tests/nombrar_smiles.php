<?php
// Lee SMILES (uno por línea) y escribe: en-pin, es-pin, en-trad, en-sis79
require __DIR__.'/../lib/IupacAlcanos.php';
while(($l=fgets(STDIN))!==false){ $l=trim($l); if($l==='')continue;
  try { $m=Molecula::desdeSmiles($l); $o=[];
    foreach([['en','pin'],['es','pin'],['en','trad'],['en','sis79']] as [$a,$b]) $o[]=(new IupacAlcanos($m,$a,$b))->nombrar()['nombre'];
    echo implode("\t",$o),"\n";
  } catch (NomenclaturaException $e) { echo "ERR\t",$e->getMessage(),"\n"; }
}
