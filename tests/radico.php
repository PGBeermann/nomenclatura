<?php
require __DIR__.'/../lib/IupacAlcanos.php';
while(($l=fgets(STDIN))!==false){ $l=trim($l); if($l==='')continue; $m=Molecula::desdeSmiles($l);
 $o=[$l]; foreach(['pin','trad','sis79'] as $e){ $o[]=(new IupacAlcanos($m,'en',$e))->nombrar()['radicofuncional'] ?? ''; } echo implode("\t",$o),"\n"; }
