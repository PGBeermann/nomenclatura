<?php
/**
 * IupacOrganica.php  ·  v4.1 (fases 1 y 2)
 * Motor de nomenclatura sustitutiva IUPAC para compuestos orgánicos con un
 * esqueleto carbonado acíclico o monocíclico (cicloalcanos, cicloalquenos y
 * benceno) y los grupos característicos más frecuentes en la enseñanza:
 *
 *   ácidos carboxílicos > ésteres > amidas > nitrilos > aldehídos > cetonas
 *   > alcoholes y fenoles > aminas      (orden de clases, P-41)
 *   y, solo como prefijos: éteres (alcoxi), halógenos (halo) y nitro.
 *
 * Fundamento normativo:
 *   Favre, H. A.; Powell, W. H. Nomenclature of Organic Chemistry. IUPAC
 *   Recommendations and Preferred Names 2013. RSC, 2014. doi:10.1039/9781849733069
 *     P-14.3.4 omisión de localizadores · P-14.5 orden alfanumérico
 *     P-16.3 multiplicadores · P-16.5 signos de inclusión · P-29 prefijos sustituyentes
 *     P-31.1.4 numeración · P-41 orden de antigüedad de clases
 *     P-44.1.1 máximo de grupos principales · P-44.1.2.2 anillo antes que cadena
 *     P-44.3 cadena más larga · P-44.4.1 insaturación · P-45.2 criterios adicionales
 *     P-61 (halógenos, nitro) · P-62 aminas · P-63 alcoholes, fenoles, éteres
 *     P-64 cetonas · P-65 ácidos y ésteres · P-66 amidas, nitrilos, aldehídos
 *
 * Fase 2: sistemas de dos anillos — fusionados con nombre retenido (P-25: pentaleno,
 * indeno, azuleno, naftaleno, heptaleno) con hidrógeno indicado, añadido y prefijos hidro
 * (P-14.7, P-31.1.4.2.4); von Baeyer (P-23.2); espiro (P-24.2); ensamblajes (P-28);
 * nomenclatura multiplicativa (P-15.3); antigüedad de anillos (P-44.2, P-44.4).
 *
 * El motor calcula los nombres a partir del grafo molecular (no los memoriza).
 * Los sistemas de tres o más anillos y los heterociclos se rechazan con un mensaje didáctico.
 *
 * Requiere PHP >= 7.4.  UNACHI · Facultad de Ciencias Naturales y Exactas · Escuela de Química.
 */
declare(strict_types=1);

if (!class_exists('NomenclaturaException', false)) {
    /** Error de dominio con mensaje apto para mostrar al estudiante. */
    final class NomenclaturaException extends RuntimeException {}
}

/* =====================================================================
 *  Grafo molecular (átomos pesados; hidrógenos implícitos)
 * ===================================================================== */
final class MolOrg
{
    public const AROM = 4;               // código de enlace aromático (orden 1,5)
    public const HALOGENOS = ['F', 'Cl', 'Br', 'I'];
    private const VALENCIA = ['C' => 4, 'N' => 3, 'O' => 2, 'F' => 1, 'Cl' => 1, 'Br' => 1, 'I' => 1];

    /** @var string[] */ public array $el = [];
    /** @var bool[]   */ public array $arom = [];
    /** @var int[]    */ public array $carga = [];
    /** @var array<int,array<int,int>> a => [b => orden] */ public array $adj = [];
    public int $n = 0;

    public function nuevo(string $el, bool $arom = false, int $carga = 0): int
    {
        $this->el[$this->n] = $el;
        $this->arom[$this->n] = $arom;
        $this->carga[$this->n] = $carga;
        $this->adj[$this->n] = [];
        $this->cacheSis = null;
        return $this->n++;
    }

    public function enlazar(int $a, int $b, int $o = 1): void
    {
        $this->adj[$a][$b] = $o;
        $this->adj[$b][$a] = $o;
        $this->cacheSis = null;
    }

    public function orden(int $a, int $b): int { return $this->adj[$a][$b] ?? 0; }

    /** @return int[] */
    public function vecinos(int $a): array { return array_keys($this->adj[$a]); }

    public function grado(int $a): int { return count($this->adj[$a]); }

    public function valenciaUsada(int $a): float
    {
        $s = 0.0;
        foreach ($this->adj[$a] as $o) { $s += $o === self::AROM ? 1.5 : $o; }
        return $s;
    }

    public function hidrogenos(int $a): int
    {
        $v = self::VALENCIA[$this->el[$a]] + ($this->el[$a] === 'N' ? $this->carga[$a] : 0)
            + ($this->el[$a] === 'O' ? $this->carga[$a] : 0);
        return (int)round($v - $this->valenciaUsada($a));
    }

    public function numCarbonos(): int
    {
        return count(array_filter($this->el, fn($e) => $e === 'C'));
    }

    public function copia(): MolOrg
    {
        $c = new MolOrg();
        $c->el = $this->el; $c->arom = $this->arom; $c->carga = $this->carga; $c->adj = $this->adj; $c->n = $this->n;
        return $c;
    }

    /* ------------------------------------------------------------------
     *  Lector SMILES (subconjunto orgánico: C N O F Cl Br I y c aromático).
     *  Acepta la salida de JSME ([CH2:3], [NH2], [C@@H]…); los mapas de
     *  átomo y la estereoquímica se ignoran.
     * ------------------------------------------------------------------ */
    public static function desdeSmiles(string $smi): MolOrg
    {
        $smi = trim($smi);
        if ($smi === '') { throw new NomenclaturaException('No hay ninguna estructura dibujada.'); }
        if (strlen($smi) > 500) { throw new NomenclaturaException('La estructura es demasiado grande.'); }
        $m = new MolOrg();
        $pila = []; $prev = -1; $bond = 0; $anillos = []; $len = strlen($smi);
        $agregar = function (string $el, bool $ar, int $q) use ($m, &$prev, &$bond): void {
            $a = $m->nuevo($el, $ar, $q);
            if ($prev >= 0) {
                $o = $bond ?: (($ar && $m->arom[$prev]) ? self::AROM : 1);
                $m->enlazar($prev, $a, $o);
            }
            $prev = $a; $bond = 0;
        };
        for ($i = 0; $i < $len; $i++) {
            $ch = $smi[$i]; $sig = $smi[$i + 1] ?? '';
            if ($ch === 'C' && $sig === 'l') { $agregar('Cl', false, 0); $i++; }
            elseif ($ch === 'B' && $sig === 'r') { $agregar('Br', false, 0); $i++; }
            elseif (in_array($ch, ['C', 'N', 'O', 'F', 'I'], true)) { $agregar($ch, false, 0); }
            elseif ($ch === 'c') { $agregar('C', true, 0); }
            elseif (in_array($ch, ['n', 'o', 's', 'p'], true)) {
                throw new NomenclaturaException('La estructura contiene un anillo heteroaromático (piridina, furano, tiofeno…). Los heterociclos se incorporarán en la fase 3.');
            }
            elseif ($ch === 'S' || $ch === 'P' || $ch === 'B') {
                throw new NomenclaturaException('La estructura contiene azufre, fósforo o boro. Por ahora el programa trabaja con C, H, N, O y halógenos.');
            }
            elseif ($ch === '[') {
                $j = strpos($smi, ']', $i);
                if ($j === false) { throw new NomenclaturaException('SMILES mal formado.'); }
                $at = substr($smi, $i + 1, $j - $i - 1);
                if (!preg_match('/^(\d*)(Cl|Br|C|N|O|F|I|c|n|o|s|S|P|B|[A-Z][a-z]?)(@{0,2})(H\d?)?([+-]\d?)?(:\d+)?$/', $at, $mm)) {
                    throw new NomenclaturaException('Átomo no admitido: [' . $at . '].');
                }
                if ($mm[1] !== '') { throw new NomenclaturaException('No se admiten isótopos.'); }
                $sym = $mm[2];
                $q = 0;
                if (!empty($mm[5])) { $q = ($mm[5][0] === '+' ? 1 : -1) * (strlen($mm[5]) > 1 ? (int)substr($mm[5], 1) : 1); }
                if (in_array($sym, ['n', 'o', 's'], true)) {
                    throw new NomenclaturaException('La estructura contiene un anillo heteroaromático. Los heterociclos se incorporarán en la fase 3.');
                }
                if (!in_array($sym, ['C', 'c', 'N', 'O', 'F', 'Cl', 'Br', 'I'], true)) {
                    throw new NomenclaturaException('La estructura contiene átomos que no son C, H, N, O ni halógenos ([' . $at . ']).');
                }
                $agregar($sym === 'c' ? 'C' : $sym, $sym === 'c', $q);
                $i = $j;
            }
            elseif ($ch === '(') { $pila[] = $prev; }
            elseif ($ch === ')') {
                if (!$pila) { throw new NomenclaturaException('SMILES mal formado.'); }
                $prev = array_pop($pila);
            }
            elseif ($ch === '=') { $bond = 2; }
            elseif ($ch === '#') { $bond = 3; }
            elseif ($ch === ':') { $bond = self::AROM; }
            elseif ($ch === '-' || $ch === '/' || $ch === '\\') { $bond = $bond ?: 1; }
            elseif ($ch === '$') { throw new NomenclaturaException('Enlace cuádruple no admitido.'); }
            elseif (ctype_digit($ch) || $ch === '%') {
                $d = $ch;
                if ($ch === '%') { $d = substr($smi, $i + 1, 2); $i += 2; }
                if (isset($anillos[$d])) {
                    [$a0, $b0] = $anillos[$d];
                    if ($a0 === $prev || isset($m->adj[$a0][$prev])) { throw new NomenclaturaException('SMILES mal formado.'); }
                    $o = $bond ?: ($b0 ?: (($m->arom[$a0] && $m->arom[$prev]) ? self::AROM : 1));
                    $m->enlazar($a0, $prev, $o);
                    unset($anillos[$d]);
                } else { $anillos[$d] = [$prev, $bond]; }
                $bond = 0;
            }
            elseif ($ch === '.') { throw new NomenclaturaException('Hay más de una molécula en el editor. Dibuje una sola estructura.'); }
            elseif ($ch === '*') { throw new NomenclaturaException('La estructura contiene átomos genéricos (R, X o *).'); }
            else { throw new NomenclaturaException('Carácter no reconocido en el SMILES: ' . $ch); }
        }
        if ($anillos || $pila) { throw new NomenclaturaException('SMILES mal formado.'); }
        if ($m->n === 0) { throw new NomenclaturaException('No hay ninguna estructura dibujada.'); }
        $m->perceibirAromaticidad();
        $m->validar();
        return $m;
    }

    /** Comprueba valencias y cargas (solo se admiten las del grupo nitro). */
    private function validar(): void
    {
        // nitro pentavalente N(=O)=O → forma con separación de cargas [N+](=O)[O-]
        foreach ($this->el as $a => $e) {
            if ($e === 'N' && $this->carga[$a] === 0 && $this->esNitroN($a)) {
                $dobles = array_keys(array_filter($this->adj[$a], fn($o) => $o === 2));
                if (count($dobles) === 2) {
                    $this->carga[$a] = 1;
                    $this->enlazar($a, $dobles[1], 1);
                    $this->carga[$dobles[1]] = -1;
                }
            }
        }
        foreach ($this->el as $a => $e) {
            if ($this->carga[$a] !== 0) {
                $esNitro = $e === 'N' && $this->carga[$a] === 1 && $this->esNitroN($a);
                $esONitro = $e === 'O' && $this->carga[$a] === -1 && count($this->adj[$a]) === 1
                    && $this->esNitroN($this->vecinos($a)[0]);
                if (!$esNitro && !$esONitro) {
                    throw new NomenclaturaException('La estructura tiene cargas eléctricas (iones). Solo se admite la representación con cargas del grupo nitro.');
                }
            }
            if ($this->hidrogenos($a) < 0) {
                $nom = ['C' => 'carbono', 'N' => 'nitrógeno', 'O' => 'oxígeno'][$e] ?? 'halógeno';
                throw new NomenclaturaException("Un átomo de $nom tiene más enlaces de los permitidos por su valencia.");
            }
        }
    }

    /** N de un grupo nitro: N unido a un C y a dos O terminales (N(=O)=O o [N+](=O)[O-]). */
    public function esNitroN(int $n): bool
    {
        if ($this->el[$n] !== 'N') { return false; }
        $o = 0; $c = 0;
        foreach ($this->adj[$n] as $b => $ord) {
            if ($this->el[$b] === 'O' && count($this->adj[$b]) === 1) { $o++; }
            elseif ($this->el[$b] === 'C') { $c++; }
            else { return false; }
        }
        return $o === 2 && $c === 1;
    }

    /* ---------------- sistemas de anillos ---------------- */

    private ?array $cacheSis = null;

    /**
     * Sistemas de anillos (fase 2): monociclos, biciclos (fusionados o con puente) y espiro de dos anillos.
     * Cada sistema: ['tipo' => mono|bic|espiro, 'atomos' => int[], 'aristas' => [[a,b]…], y además
     *   mono:   'ciclo'   => átomos en orden;
     *   bic:    'cabezas' => [h1, h2], 'puentes' => tres listas de átomos interiores ordenadas de h1 a h2;
     *   espiro: 'espiro'  => átomo espiro, 'ciclos' => dos ciclos ordenados que empiezan en el átomo espiro.
     * Lanza excepción didáctica para heterociclos y sistemas de tres o más anillos.
     */
    public function sistemas(): array
    {
        if ($this->cacheSis !== null) { return $this->cacheSis; }
        $disc = []; $low = []; $t = 0; $pilaE = []; $comps = [];
        $dfs = function (int $u, int $p) use (&$dfs, &$disc, &$low, &$t, &$pilaE, &$comps): void {
            $disc[$u] = $low[$u] = ++$t;
            foreach ($this->adj[$u] as $v => $o) {
                if ($v === $p) { continue; }
                if (!isset($disc[$v])) {
                    $pilaE[] = [$u, $v];
                    $dfs($v, $u);
                    $low[$u] = min($low[$u], $low[$v]);
                    if ($low[$v] >= $disc[$u]) {
                        $c = [];
                        do { $e = array_pop($pilaE); $c[] = $e; } while ($e !== [$u, $v]);
                        $comps[] = $c;
                    }
                } elseif ($disc[$v] < $disc[$u]) {
                    $pilaE[] = [$u, $v];
                    $low[$u] = min($low[$u], $disc[$v]);
                }
            }
        };
        for ($a = 0; $a < $this->n; $a++) { if (!isset($disc[$a])) { $dfs($a, -1); } }
        $cic = [];
        foreach ($comps as $c) {
            $V = [];
            foreach ($c as [$x, $y]) { $V[$x] = true; $V[$y] = true; }
            $nE = count($c); $nV = count($V);
            if ($nE < 3) { continue; }                       // enlace acíclico
            foreach (array_keys($V) as $a) {
                if ($this->el[$a] !== 'C') {
                    throw new NomenclaturaException('El anillo contiene un heteroátomo (' . $this->el[$a] . '): es un heterociclo (lactona, lactama, éter cíclico, amina cíclica…). Los heterociclos se incorporarán en la fase 3.');
                }
            }
            if ($nE > $nV + 1) {
                throw new NomenclaturaException('La estructura tiene un sistema de tres o más anillos (antraceno, fenantreno, adamantano, esteroides…). Por ahora el programa admite hasta dos anillos por sistema.');
            }
            $cic[] = ['V' => $V, 'E' => $c, 'nE' => $nE, 'nV' => $nV];
        }
        // agrupar componentes que comparten un átomo (espiro)
        $grupo = range(0, max(0, count($cic) - 1));
        $raizG = function (int $i) use (&$grupo, &$raizG): int { return $grupo[$i] === $i ? $i : ($grupo[$i] = $raizG($grupo[$i])); };
        $duenio = [];
        foreach ($cic as $i => $c) {
            foreach (array_keys($c['V']) as $a) {
                if (isset($duenio[$a])) { $grupo[$raizG($i)] = $raizG($duenio[$a]); } else { $duenio[$a] = $i; }
            }
        }
        $grupos = [];
        foreach ($cic as $i => $c) { $grupos[$raizG($i)][] = $i; }
        $sis = [];
        foreach ($grupos as $miembros) {
            if (count($miembros) === 1) {
                $c = $cic[$miembros[0]];
                $sis[] = $c['nE'] === $c['nV'] ? $this->sistemaMono($c) : $this->sistemaBic($c);
            } elseif (count($miembros) === 2 && $cic[$miembros[0]]['nE'] === $cic[$miembros[0]]['nV'] && $cic[$miembros[1]]['nE'] === $cic[$miembros[1]]['nV']) {
                $A = $this->sistemaMono($cic[$miembros[0]]); $B = $this->sistemaMono($cic[$miembros[1]]);
                $com = array_values(array_intersect($A['atomos'], $B['atomos']));
                if (count($com) !== 1) { throw new NomenclaturaException('Sistema de anillos no admitido.'); }
                $s0 = $com[0];
                $rot = function (array $ciclo) use ($s0): array { $i = array_search($s0, $ciclo, true); return array_merge(array_slice($ciclo, $i), array_slice($ciclo, 0, $i)); };
                $sis[] = ['tipo' => 'espiro', 'atomos' => array_values(array_unique(array_merge($A['atomos'], $B['atomos']))),
                    'aristas' => array_merge($A['aristas'], $B['aristas']), 'espiro' => $s0, 'ciclos' => [$rot($A['ciclo']), $rot($B['ciclo'])]];
            } else {
                throw new NomenclaturaException('La estructura tiene un sistema de más de dos anillos unidos en espiro o combinados con puentes. Por ahora el programa admite hasta dos anillos por sistema.');
            }
        }
        return $this->cacheSis = $sis;
    }

    private function sistemaMono(array $c): array
    {
        $V = $c['V'];
        $ini = array_key_first($V); $orden = [$ini]; $prev = -1; $act = $ini;
        while (true) {
            $sig = -1;
            foreach ($this->adj[$act] as $b => $o) { if (isset($V[$b]) && $b !== $prev) { $sig = $b; break; } }
            if ($sig === $ini || $sig === -1) { break; }
            $orden[] = $sig; $prev = $act; $act = $sig;
        }
        $ar = [];
        $N = count($orden);
        for ($i = 0; $i < $N; $i++) { $ar[] = [$orden[$i], $orden[($i + 1) % $N]]; }
        return ['tipo' => 'mono', 'atomos' => $orden, 'aristas' => $ar, 'ciclo' => $orden];
    }

    private function sistemaBic(array $c): array
    {
        $V = $c['V'];
        $vec = [];
        foreach ($c['E'] as [$x, $y]) { $vec[$x][] = $y; $vec[$y][] = $x; }
        $cab = array_values(array_filter(array_keys($V), fn($a) => count($vec[$a]) === 3));
        if (count($cab) !== 2) { throw new NomenclaturaException('Sistema bicíclico no reconocido.'); }
        [$h1, $h2] = $cab;
        $puentes = [];
        foreach ($vec[$h1] as $n) {
            $ruta = []; $prev = $h1; $act = $n;
            while ($act !== $h2) {
                $ruta[] = $act;
                $sig = null;
                foreach ($vec[$act] as $w) { if ($w !== $prev) { $sig = $w; } }
                $prev = $act; $act = $sig;
            }
            $puentes[] = $ruta;
        }
        return ['tipo' => 'bic', 'atomos' => array_keys($V), 'aristas' => $c['E'], 'cabezas' => [$h1, $h2], 'puentes' => $puentes];
    }

    /** Ciclos elementales (lista de átomos ordenados) de todos los sistemas. */
    public function anillos(): array
    {
        $r = [];
        foreach ($this->sistemas() as $s) {
            if ($s['tipo'] === 'mono') { $r[] = $s['ciclo']; }
            elseif ($s['tipo'] === 'espiro') { $r[] = $s['ciclos'][0]; $r[] = $s['ciclos'][1]; }
            else {
                [$h1, $h2] = $s['cabezas']; $P = $s['puentes'];
                for ($i = 0; $i < 3; $i++) {
                    for ($j = $i + 1; $j < 3; $j++) { $r[] = array_merge([$h1], $P[$i], [$h2], array_reverse($P[$j])); }
                }
            }
        }
        return $r;
    }

    /**
     * Emparejamientos perfectos (estructuras de Kekulé) de un conjunto de átomos usando las aristas dadas.
     * @return array<int,array<int,array{0:int,1:int}>> hasta $limite emparejamientos
     */
    public static function emparejamientos(array $atomos, array $aristas, int $limite = 64): array
    {
        $vec = [];
        foreach ($atomos as $a) { $vec[$a] = []; }
        foreach ($aristas as [$x, $y]) { if (isset($vec[$x], $vec[$y])) { $vec[$x][] = $y; $vec[$y][] = $x; } }
        $res = [];
        $rec = function (array $libres, array $acc) use (&$rec, &$res, $vec, $limite): void {
            if (count($res) >= $limite) { return; }
            if (!$libres) { $res[] = $acc; return; }
            $a = array_key_first($libres);
            foreach ($vec[$a] as $b) {
                if (!isset($libres[$b])) { continue; }
                $l2 = $libres; unset($l2[$a], $l2[$b]);
                $rec($l2, array_merge($acc, [[$a, $b]]));
            }
        };
        $rec(array_fill_keys($atomos, true), []);
        return $res;
    }

    /** Percepción de aromaticidad: benceno de Kekulé → aromático; validación del subgrafo aromático. */
    private function perceibirAromaticidad(): void
    {
        $this->cacheSis = null;
        $sis = $this->sistemas();
        foreach ($sis as $sm) {
            if ($sm['tipo'] !== 'mono') { continue; }           // solo el benceno aislado se trata como aromático
            $r = $sm['ciclo'];
            $N = count($r);
            if ($N !== 6) { continue; }
            $dobles = 0; $alterna = true; $todosArom = true;
            for ($i = 0; $i < $N; $i++) {
                $o = $this->orden($r[$i], $r[($i + 1) % $N]);
                if ($o !== self::AROM) { $todosArom = false; }
                if ($o === 2) { $dobles++; }
                $o2 = $this->orden($r[($i + 1) % $N], $r[($i + 2) % $N]);
                if ($o === 2 && $o2 === 2) { $alterna = false; }
            }
            if (!$todosArom && $dobles === 3 && $alterna) {
                $ok = true;
                foreach ($r as $a) {
                    foreach ($this->adj[$a] as $b => $o) { if (!in_array($b, $r, true) && $o >= 2) { $ok = false; } }
                }
                if ($ok) {
                    for ($i = 0; $i < $N; $i++) { $this->enlazar($r[$i], $r[($i + 1) % $N], self::AROM); }
                    foreach ($r as $a) { $this->arom[$a] = true; }
                }
            }
        }
        // cada átomo aromático debe pertenecer a un sistema de anillos y el conjunto debe ser kekulizable
        $sisDe = [];
        foreach ($sis as $k => $s) { foreach ($s['atomos'] as $a) { $sisDe[$a] = $k; } }
        foreach ($this->arom as $a => $ar) {
            foreach ($this->adj[$a] as $b => $o) {
                if ($o === self::AROM && (!isset($sisDe[$a], $sisDe[$b]) || $sisDe[$a] !== $sisDe[$b])) {
                    throw new NomenclaturaException('Enlace aromático fuera de un anillo.');
                }
            }
            if ($ar && !isset($sisDe[$a])) { throw new NomenclaturaException('Hay átomos aromáticos fuera de un anillo.'); }
        }
        foreach ($sis as $s) {
            $ar = array_values(array_filter($s['atomos'], fn($a) => $this->arom[$a]));
            if (!$ar) { continue; }
            $ars = array_values(array_filter($s['aristas'], fn($e) => $this->orden($e[0], $e[1]) === self::AROM));
            $k = self::emparejamientos($ar, $ars, 1);
            if (!$k) { throw new NomenclaturaException('El anillo aromático dibujado no es válido (no admite una estructura de Kekulé).'); }
            // salvo el benceno aislado, los sistemas aromáticos (naftaleno, azuleno, benzociclobuteno,
            // ciclooctatetraeno…) se representan con su estructura de Kekulé explícita
            if (!($s['tipo'] === 'mono' && count($s['atomos']) === 6)) {
                foreach ($ars as [$x, $y]) { $this->enlazar($x, $y, 1); }
                foreach ($k[0] as [$x, $y]) { $this->enlazar($x, $y, 2); }
                foreach ($ar as $a) { $this->arom[$a] = false; }
            }
        }
        $this->cacheSis = null;
    }

    /* ---------------- salida ---------------- */

    /** Fórmula molecular (orden de Hill). */
    public function formula(): string
    {
        $cnt = []; $h = 0;
        foreach ($this->el as $a => $e) { $cnt[$e] = ($cnt[$e] ?? 0) + 1; $h += max(0, $this->hidrogenos($a)); }
        $f = '';
        if (!empty($cnt['C'])) { $f .= 'C' . ($cnt['C'] > 1 ? $cnt['C'] : ''); }
        if ($h > 0) { $f .= 'H' . ($h > 1 ? $h : ''); }
        $otros = array_diff(array_keys($cnt), ['C']);
        sort($otros, SORT_STRING);
        foreach ($otros as $e) { $f .= $e . ($cnt[$e] > 1 ? $cnt[$e] : ''); }
        return $f;
    }

    /** SMILES con mapas de átomo opcionales ($mapas[átomo] = número mostrado). */
    public function smiles(array $mapas = [], int $inicio = 0): string
    {
        $vis = []; $cierres = []; $padre = [$inicio => -1];
        $dfs = function (int $a) use (&$dfs, &$vis, &$cierres, &$padre): void {
            $vis[$a] = true;
            foreach ($this->adj[$a] as $b => $o) {
                if ($b === $padre[$a]) { continue; }
                if (isset($vis[$b])) {
                    $k = min($a, $b) . '-' . max($a, $b);
                    if (!isset($cierres[$k])) { $cierres[$k] = [$b, $a]; }
                    continue;
                }
                $padre[$b] = $a;
                $dfs($b);
            }
        };
        $dfs($inicio);
        $dig = []; $d = 1;
        foreach ($cierres as $k => [$x, $y]) { $dig[$x][] = [$d, $y]; $dig[$y][] = [$d, $x]; $d++; }
        $simEnlace = function (int $a, int $b): string {
            $o = $this->orden($a, $b);
            if ($o === 2) { return '='; }
            if ($o === 3) { return '#'; }
            if ($o === 1 && $this->arom[$a] && $this->arom[$b]) { return '-'; }
            return '';
        };
        $atomo = function (int $a) use ($mapas): string {
            $e = $this->el[$a]; $q = $this->carga[$a];
            $sym = $this->arom[$a] ? strtolower($e) : $e;
            if (isset($mapas[$a]) || $q !== 0) {
                $h = max(0, $this->hidrogenos($a));
                $s = '[' . $sym . ($h > 0 ? 'H' . ($h > 1 ? $h : '') : '');
                if ($q !== 0) { $s .= $q > 0 ? '+' : '-'; }
                if (isset($mapas[$a])) { $s .= ':' . $mapas[$a]; }
                return $s . ']';
            }
            return $sym;
        };
        $emit = function (int $a, int $p) use (&$emit, $dig, $cierres, $simEnlace, $atomo): string {
            $s = $atomo($a);
            foreach ($dig[$a] ?? [] as [$dg, $otro]) {
                $s .= $simEnlace($a, $otro) . ($dg > 9 ? '%' . $dg : (string)$dg);
            }
            $hijos = [];
            foreach ($this->adj[$a] as $b => $o) {
                if ($b === $p) { continue; }
                if (isset($cierres[min($a, $b) . '-' . max($a, $b)])) { continue; }
                $hijos[] = $b;
            }
            $nh = count($hijos);
            foreach ($hijos as $i => $b) {
                $sub = $simEnlace($a, $b) . $emit($b, $a);
                $s .= ($i < $nh - 1) ? '(' . $sub . ')' : $sub;
            }
            return $s;
        };
        // cierres: el dígito solo se escribe una vez con su símbolo de enlace en la apertura
        return $emit($inicio, -1);
    }
}

/* =====================================================================
 *  Motor de nomenclatura
 * ===================================================================== */
final class IupacOrganica
{
    private const RAIZ = [
        1 => 'met', 2 => 'et', 3 => 'prop', 4 => 'but', 5 => 'pent', 6 => 'hex', 7 => 'hept', 8 => 'oct',
        9 => 'non', 10 => 'dec', 11 => 'undec', 12 => 'dodec', 13 => 'tridec', 14 => 'tetradec',
        15 => 'pentadec', 16 => 'hexadec', 17 => 'heptadec', 18 => 'octadec', 19 => 'nonadec', 20 => 'icos',
    ];
    private const MULT_SIMPLE = [2 => 'di', 3 => 'tri', 4 => 'tetra', 5 => 'penta', 6 => 'hexa', 7 => 'hepta', 8 => 'octa',
        9 => 'nona', 10 => 'deca', 11 => 'undeca', 12 => 'dodeca', 13 => 'trideca', 14 => 'tetradeca', 15 => 'pentadeca',
        16 => 'hexadeca', 17 => 'heptadeca', 18 => 'octadeca', 19 => 'nonadeca', 20 => 'icosa'];
    private const MULT_COMPLEJO = [2 => 'bis', 3 => 'tris', 4 => 'tetrakis', 5 => 'pentakis', 6 => 'hexakis', 7 => 'heptakis',
        8 => 'octakis', 9 => 'nonakis', 10 => 'decakis'];
    private const HALO = [
        'es' => ['F' => 'fluoro', 'Cl' => 'cloro', 'Br' => 'bromo', 'I' => 'yodo'],
        'en' => ['F' => 'fluoro', 'Cl' => 'chloro', 'Br' => 'bromo', 'I' => 'iodo'],
    ];
    private const HALURO = [
        'es' => ['F' => 'fluoruro', 'Cl' => 'cloruro', 'Br' => 'bromuro', 'I' => 'yoduro'],
        'en' => ['F' => 'fluoride', 'Cl' => 'chloride', 'Br' => 'bromide', 'I' => 'iodide'],
    ];

    /** Orden de antigüedad de las clases que pueden expresarse como sufijo (P-41). */
    public const CLASES = ['acido', 'ester', 'amida', 'nitrilo', 'aldehido', 'cetona', 'alcohol', 'amina'];
    public const NOMBRE_CLASE = [
        'acido' => 'ácido carboxílico', 'ester' => 'éster', 'amida' => 'amida', 'nitrilo' => 'nitrilo',
        'aldehido' => 'aldehído', 'cetona' => 'cetona', 'alcohol' => 'alcohol o fenol', 'amina' => 'amina',
        'eter' => 'éter', 'halogeno' => 'halógeno', 'nitro' => 'nitro', 'alqueno' => 'doble enlace C=C', 'alquino' => 'triple enlace C≡C',
    ];
    /** [sufijo en cadena, sufijo en anillo, prefijo, regla] */
    private const INFO_CLASE = [
        'acido'    => ['«ácido …-oico»', '«ácido …carboxílico»', 'carboxi', 'P-65.1'],
        'ester'    => ['«…-oato de alquilo»', '«…carboxilato de alquilo»', 'alcoxicarbonil / aciloxi', 'P-65.6'],
        'amida'    => ['«-amida»', '«-carboxamida»', 'carbamoil / acilamino', 'P-66.1'],
        'nitrilo'  => ['«-nitrilo»', '«-carbonitrilo»', 'ciano', 'P-66.5'],
        'aldehido' => ['«-al»', '«-carbaldehído»', 'oxo / formil', 'P-66.6'],
        'cetona'   => ['«-ona»', '«-ona»', 'oxo', 'P-64'],
        'alcohol'  => ['«-ol»', '«-ol»', 'hidroxi', 'P-63.1'],
        'amina'    => ['«-amina»', '«-amina»', 'amino', 'P-62'],
    ];

    private MolOrg $mol;
    private string $lang;
    private string $estilo;              // pin | trad | sis79
    /** @var array<int,array> carbonos funcionales (C de COOH, COOR, CONR2, CN) */
    private array $fc = [];
    /** @var array<int,int> carbonos de aldehído => carbono de anclaje (-1 si no hay) */
    private array $ald = [];
    /** @var array<int,int> átomo => índice del sistema de anillos */
    private array $anilloDe = [];
    /** @var array<int,array> sistemas de anillos (MolOrg::sistemas) */
    private array $sistemas = [];
    /** @var array<int,array> caché de numeraciones y datos de nomenclatura por sistema */
    private array $infoSis = [];
    private ?string $principal = null;
    /** @var array<string,int> clases presentes (incluye éter, halógeno, nitro, insaturaciones) */
    public array $clases = [];
    private array $memo = [];
    /** Prefijos N (N-metil…) de sufijos amina/amida: se numeran N, N1, N2… */
    private int $nSufN = 0;

    public function __construct(MolOrg $mol, string $lang = 'es', string $estilo = 'pin')
    {
        $this->mol = $mol;
        $this->lang = $lang;
        $this->estilo = $estilo;
        $this->sistemas = $mol->sistemas();
        foreach ($this->sistemas as $i => $r) { foreach ($r['atomos'] as $a) { $this->anilloDe[$a] = $i; } }
        $this->clasificar();
    }

    public function principal(): ?string { return $this->principal; }

    /** Nombres comunes o retenidos aceptados en nomenclatura general (no PIN), indexados por el PIN en inglés. */
    public const COMUNES = [
        'propan-2-one' => [['acetona'], ['acetone']],
        'acetic acid' => [['ácido etanoico'], ['ethanoic acid']],
        'formic acid' => [['ácido metanoico'], ['methanoic acid']],
        'formaldehyde' => [['metanal'], ['methanal']],
        'acetaldehyde' => [['etanal'], ['ethanal']],
        'formamide' => [['metanamida'], ['methanamide']],
        'acetamide' => [['etanamida'], ['ethanamide']],
        'acetonitrile' => [['etanonitrilo'], ['ethanenitrile']],
        'methyl acetate' => [['etanoato de metilo'], ['methyl ethanoate']],
        'ethyl acetate' => [['etanoato de etilo'], ['ethyl ethanoate']],
        'propanoic acid' => [['ácido propiónico'], ['propionic acid']],
        'butanoic acid' => [['ácido butírico'], ['butyric acid']],
        'propanedioic acid' => [['ácido malónico'], ['malonic acid']],
        'butanedioic acid' => [['ácido succínico'], ['succinic acid']],
        'hexanedioic acid' => [['ácido adípico'], ['adipic acid']],
        'oxalic acid' => [['ácido etanodioico'], ['ethanedioic acid']],
        'prop-2-enoic acid' => [['ácido acrílico'], ['acrylic acid']],
        '2-hydroxypropanoic acid' => [['ácido láctico'], ['lactic acid']],
        'benzoic acid' => [['ácido bencenocarboxílico'], ['benzenecarboxylic acid']],
        '2-hydroxybenzoic acid' => [['ácido salicílico'], ['salicylic acid']],
        'benzene-1,2-dicarboxylic acid' => [['ácido ftálico'], ['phthalic acid']],
        'benzene-1,4-dicarboxylic acid' => [['ácido tereftálico'], ['terephthalic acid']],
        'benzaldehyde' => [['bencenocarbaldehído'], ['benzenecarbaldehyde']],
        'benzamide' => [['bencenocarboxamida'], ['benzenecarboxamide']],
        'benzonitrile' => [['bencenocarbonitrilo'], ['benzenecarbonitrile']],
        'N-phenylacetamide' => [['acetanilida'], ['acetanilide']],
        'aniline' => [['bencenamina', 'fenilamina'], ['benzenamine', 'phenylamine']],
        'phenol' => [['hidroxibenceno'], ['hydroxybenzene']],
        '2-methylphenol' => [['o-cresol'], ['o-cresol']],
        '3-methylphenol' => [['m-cresol'], ['m-cresol']],
        '4-methylphenol' => [['p-cresol'], ['p-cresol']],
        'benzene-1,2-diol' => [['catecol', 'pirocatecol'], ['catechol', 'pyrocatechol']],
        'benzene-1,3-diol' => [['resorcinol'], ['resorcinol']],
        'benzene-1,4-diol' => [['hidroquinona'], ['hydroquinone']],
        'toluene' => [['metilbenceno'], ['methylbenzene']],
        'anisole' => [['metoxibenceno'], ['methoxybenzene']],
        '1,2-xylene' => [['o-xileno', '1,2-dimetilbenceno'], ['o-xylene', '1,2-dimethylbenzene']],
        '1,3-xylene' => [['m-xileno', '1,3-dimetilbenceno'], ['m-xylene', '1,3-dimethylbenzene']],
        '1,4-xylene' => [['p-xileno', '1,4-dimetilbenceno'], ['p-xylene', '1,4-dimethylbenzene']],
        'ethenylbenzene' => [['estireno'], ['styrene']],
        '1-phenylethan-1-one' => [['acetofenona'], ['acetophenone']],
        'diphenylmethanone' => [['benzofenona'], ['benzophenone']],
        'cyclohexa-2,5-diene-1,4-dione' => [['p-benzoquinona', '1,4-benzoquinona'], ['p-benzoquinone', '1,4-benzoquinone']],
        'ethane-1,2-diol' => [['etilenglicol'], ['ethylene glycol']],
        'propane-1,2,3-triol' => [['glicerol', 'glicerina'], ['glycerol']],
        'propan-2-ol' => [['isopropanol'], ['isopropanol']],
        'ethene' => [['etileno'], ['ethylene']],
        'propene' => [['propileno', 'prop-1-eno'], ['propylene', 'prop-1-ene']],
        'acetylene' => [['etino'], ['ethyne']],
        'propyne' => [['prop-1-ino'], ['prop-1-yne']],
        'buta-1,3-diene' => [['1,3-butadieno'], ['1,3-butadiene']],
        '2-methylbuta-1,3-diene' => [['isopreno'], ['isoprene']],
        'ethoxyethane' => [['éter etílico', 'éter dietílico'], ['ether']],
        'trichloromethane' => [['cloroformo'], ['chloroform']],
        'tetrachloromethane' => [['tetracloruro de carbono'], ['carbon tetrachloride']],
        'dichloromethane' => [['cloruro de metileno'], ['methylene chloride']],
        'prop-2-enal' => [['acroleína'], ['acrolein']],
        'prop-2-enenitrile' => [['acrilonitrilo'], ['acrylonitrile']],
        // fase 2: sistemas de dos anillos
        'naphthalen-1-ol' => [['1-naftol', 'α-naftol'], ['1-naphthol']],
        'naphthalen-2-ol' => [['2-naftol', 'β-naftol'], ['2-naphthol']],
        'naphthalene-1-carboxylic acid' => [['ácido 1-naftoico'], ['1-naphthoic acid']],
        'naphthalene-2-carboxylic acid' => [['ácido 2-naftoico'], ['2-naphthoic acid']],
        'naphthalene-1,4-dione' => [['1,4-naftoquinona'], ['1,4-naphthoquinone']],
        'decahydronaphthalene' => [['decalina'], ['decalin']],
        '1,2,3,4-tetrahydronaphthalene' => [['tetralina'], ['tetralin']],
        '3,4-dihydronaphthalen-1(2H)-one' => [['1-tetralona', 'α-tetralona'], ['1-tetralone']],
        '2,3-dihydro-1H-indene' => [['indano'], ['indane']],
        '2,3-dihydro-1H-inden-1-one' => [['1-indanona'], ['1-indanone']],
        'bicyclo[2.2.1]heptane' => [['norbornano'], ['norbornane']],
        'bicyclo[2.2.1]hept-2-ene' => [['norborneno'], ['norbornene']],
        '1,7,7-trimethylbicyclo[2.2.1]heptan-2-one' => [['alcanfor'], ['camphor']],
        'bicyclo[4.2.0]octa-1,3,5-triene' => [['benzociclobuteno'], ['benzocyclobutene']],
        "1,1'-biphenyl" => [['bifenilo'], ['biphenyl']],
        "1,1'-bi(cyclohexane)" => [['biciclohexilo'], ['bicyclohexyl']],
        "1,1'-methylenedibenzene" => [['difenilmetano'], ['diphenylmethane']],
        "1,1'-oxydibenzene" => [['difenil éter', 'éter difenílico'], ['diphenyl ether']],
        "1,1'-(ethane-1,2-diyl)dibenzene" => [['bibencilo', '1,2-difeniletano'], ['bibenzyl', '1,2-diphenylethane']],
        "4,4'-(propane-2,2-diyl)diphenol" => [['bisfenol A'], ['bisphenol A']],
        "4,4'-methylenedianiline" => [["4,4'-diaminodifenilmetano"], ["4,4'-diaminodiphenylmethane"]],
        "[1,1'-binaphthalene]-2,2'-diol" => [['BINOL'], ['BINOL']],
        'N-phenylaniline' => [['difenilamina'], ['diphenylamine']],
    ];

    /* ======================= clasificación de grupos ======================= */

    private function clasificar(): void
    {
        $m = $this->mol;
        $cl = [];
        foreach ($m->el as $a => $e) {
            if ($e === 'O') {
                $vs = $m->adj[$a];
                if (count($vs) === 0) { throw new NomenclaturaException('Hay un átomo de oxígeno aislado.'); }
                foreach ($vs as $b => $o) {
                    if ($m->el[$b] === 'O') { throw new NomenclaturaException('La estructura contiene un enlace O–O (peróxido). No se incluye en esta fase.'); }
                    if ($m->el[$b] === 'N' && !$m->esNitroN($b)) { throw new NomenclaturaException('La estructura contiene un enlace N–O (hidroxilamina, óxido de amina…). No se incluye en esta fase.'); }
                    if (MolOrg::HALOGENOS !== [] && in_array($m->el[$b], MolOrg::HALOGENOS, true)) { throw new NomenclaturaException('Enlace oxígeno–halógeno no admitido.'); }
                }
            } elseif ($e === 'N') {
                if ($m->esNitroN($a)) { $cl['nitro'] = ($cl['nitro'] ?? 0) + 1; continue; }
                foreach ($m->adj[$a] as $b => $o) {
                    if ($m->el[$b] !== 'C') { throw new NomenclaturaException('El nitrógeno solo puede estar unido a carbonos (se excluyen hidracinas, azo, N-halo…).'); }
                    if ($o === 2) { throw new NomenclaturaException('La estructura contiene un doble enlace C=N (imina u oxima). Las iminas se incorporarán en una fase posterior.'); }
                }
            } elseif (in_array($e, MolOrg::HALOGENOS, true)) {
                $vs = $m->vecinos($a);
                if (count($vs) !== 1 || $m->el[$vs[0]] !== 'C') { throw new NomenclaturaException('Cada halógeno debe estar unido a un solo átomo de carbono.'); }
                $cl['halogeno'] = ($cl['halogeno'] ?? 0) + 1;
            }
        }
        foreach ($m->el as $c => $e) {
            if ($e !== 'C') { continue; }
            $oxo = 0; $oh = 0; $or = []; $nN = []; $nTri = 0; $hal = 0; $carb = []; $nitro = 0;
            foreach ($m->adj[$c] as $b => $o) {
                $eb = $m->el[$b];
                if ($eb === 'C') { $carb[] = $b; }
                elseif ($eb === 'O') {
                    if ($o === 2) { $oxo++; }
                    elseif (count($m->adj[$b]) === 1) { $oh++; }
                    else { $or[] = $b; }
                } elseif ($eb === 'N') {
                    if ($m->esNitroN($b)) { $nitro++; continue; }
                    if ($o === 3) { $nTri++; } else { $nN[] = $b; }
                } else { $hal++; }
            }
            if ($oxo > 1) { throw new NomenclaturaException('Un carbono tiene dos dobles enlaces a oxígeno.'); }
            if ($nTri > 0) {
                if ($oxo || $oh || $or || $nN || $hal) { throw new NomenclaturaException('Grupo funcional no admitido en esta fase (cianato, cianamida…).'); }
                $n = array_keys(array_filter($m->adj[$c], fn($o) => $o === 3))[0];
                if (count($m->adj[$n]) !== 1) { throw new NomenclaturaException('Grupo C≡N no terminal no admitido.'); }
                $this->fc[$c] = ['t' => 'nitrilo', 'anc' => $carb[0] ?? -1, 'n' => $n];
                continue;
            }
            if ($oxo === 1) {
                $het = $oh + count($or) + count($nN) + $hal + $nitro;
                if ($het === 0) {
                    if ($m->hidrogenos($c) >= 1) { $this->ald[$c] = $carb[0] ?? -1; }
                    continue;                                        // cetona: carbono del esqueleto
                }
                if ($het > 1 || $hal || $nitro) {
                    throw new NomenclaturaException('El carbono carbonílico lleva dos heteroátomos o un halógeno (haluro de acilo, carbonato, carbamato, urea, anhídrido…). Estos grupos se incorporarán en una fase posterior.');
                }
                if ($oh) { $this->fc[$c] = ['t' => 'acido', 'anc' => $carb[0] ?? -1]; continue; }
                if ($or) {
                    $o = $or[0];
                    $r = null;
                    foreach ($m->adj[$o] as $b => $x) { if ($b !== $c) { $r = $b; } }
                    if ($m->el[$r] !== 'C') { throw new NomenclaturaException('Éster no admitido.'); }
                    foreach ($m->adj[$r] as $b => $x) {
                        if ($m->el[$b] === 'O' && $x === 2) { throw new NomenclaturaException('La estructura contiene un anhídrido de ácido. Se incorporará en una fase posterior.'); }
                    }
                    $this->fc[$c] = ['t' => 'ester', 'anc' => $carb[0] ?? -1, 'o' => $o, 'r' => $r];
                    continue;
                }
                $n = $nN[0];
                $ns = [];
                foreach ($m->adj[$n] as $b => $x) { if ($b !== $c) { $ns[] = $b; } }
                $this->fc[$c] = ['t' => 'amida', 'anc' => $carb[0] ?? -1, 'n' => $n, 'ns' => $ns];
            }
        }
        // amidas: el N no puede estar unido a dos carbonos acilo (imida)
        foreach ($this->fc as $c => $d) {
            if ($d['t'] !== 'amida') { continue; }
            foreach ($d['ns'] as $b) {
                if (isset($this->fc[$b])) { throw new NomenclaturaException('La estructura contiene una imida (N unido a dos grupos acilo). Se incorporará en una fase posterior.'); }
            }
        }
        foreach ($this->fc as $c => $d) {
            $cl[$d['t']] = ($cl[$d['t']] ?? 0) + 1;
        }
        foreach ($this->ald as $c => $x) { $cl['aldehido'] = ($cl['aldehido'] ?? 0) + 1; }
        foreach ($m->el as $a => $e) {
            if ($e === 'O' && count($m->adj[$a]) === 1) {
                $c = $m->vecinos($a)[0];
                if ($m->el[$c] !== 'C' || isset($this->fc[$c])) { continue; }
                if ($m->orden($a, $c) === 2) {
                    if (!isset($this->ald[$c])) { $cl['cetona'] = ($cl['cetona'] ?? 0) + 1; }
                } else { $cl['alcohol'] = ($cl['alcohol'] ?? 0) + 1; }
            } elseif ($e === 'O' && count($m->adj[$a]) === 2) {
                $es = false;
                foreach ($m->adj[$a] as $b => $x) { if (isset($this->fc[$b]) && $this->fc[$b]['t'] === 'ester' && $this->fc[$b]['o'] === $a) { $es = true; } }
                if (!$es) { $cl['eter'] = ($cl['eter'] ?? 0) + 1; }
            } elseif ($e === 'N' && !$m->esNitroN($a) && !$this->esNAmida($a) && !$this->esNNitrilo($a)) {
                $cl['amina'] = ($cl['amina'] ?? 0) + 1;
            }
        }
        foreach ($m->adj as $a => $vs) {
            foreach ($vs as $b => $o) {
                if ($a < $b && $m->el[$a] === 'C' && $m->el[$b] === 'C') {
                    if ($o === 2) { $cl['alqueno'] = ($cl['alqueno'] ?? 0) + 1; }
                    if ($o === 3) { $cl['alquino'] = ($cl['alquino'] ?? 0) + 1; }
                }
            }
        }
        $this->clases = $cl;
        foreach (self::CLASES as $k) { if (!empty($cl[$k])) { $this->principal = $k; break; } }
    }

    private function esNAmida(int $n): bool
    {
        foreach ($this->mol->adj[$n] as $b => $o) { if (isset($this->fc[$b]) && $this->fc[$b]['t'] === 'amida' && $this->fc[$b]['n'] === $n) { return true; } }
        return false;
    }

    private function esNNitrilo(int $n): bool
    {
        foreach ($this->mol->adj[$n] as $b => $o) { if ($o === 3) { return true; } }
        return false;
    }

    /** Tipo de carbono funcional en el contexto indicado (el aldehído principal actúa como tal solo en la estructura principal). */
    private function tipoFC(int $x, bool $padre): ?string
    {
        if (isset($this->fc[$x])) { return $this->fc[$x]['t']; }
        if ($padre && $this->principal === 'aldehido' && isset($this->ald[$x])) { return 'aldehido'; }
        return null;
    }

    private function anclaFC(int $x): int
    {
        return isset($this->fc[$x]) ? $this->fc[$x]['anc'] : ($this->ald[$x] ?? -1);
    }

    /** Carbono del esqueleto (candidato a formar parte de cadenas). */
    private function esEsqueleto(int $x, bool $padre): bool
    {
        return $this->mol->el[$x] === 'C' && $this->tipoFC($x, $padre) === null;
    }

    /* ======================= utilidades léxicas ======================= */

    private function es(): bool { return $this->lang === 'es'; }

    private function raiz(int $n): string
    {
        if ($n > 20) { throw new NomenclaturaException("La cadena principal tiene $n carbonos; el programa admite hasta 20 (icosano)."); }
        $r = self::RAIZ[$n];
        if (!$this->es() && $n <= 2) { $r .= 'h'; }       // meth, eth
        return $r;
    }

    private function yl(): string { return $this->es() ? 'il' : 'yl'; }

    /** Clave alfanumérica (P-14.5): se ignoran sec-, terc-/tert-, localizadores y signos. */
    private static function claveAlfa(string $nombre): string
    {
        $s = preg_replace('/\b(tert|terc|sec)-/', '', $nombre);
        $s = strtr($s, ['á' => 'a', 'é' => 'e', 'í' => 'i', 'ó' => 'o', 'ú' => 'u']);
        return preg_replace('/[^a-z]/', '', strtolower($s));
    }

    private function info(string $nombre, string $tipo, int $longitud = 0, array $extra = []): array
    {
        return ['nombre' => $nombre, 'tipo' => $tipo, 'clave' => self::claveAlfa($nombre), 'longitud' => $longitud] + $extra;
    }

    /** Encierra un prefijo en ( ), [ ] o { } según lo que ya contenga (P-16.5.4). */
    private static function encerrar(string $nm): string
    {
        if (strpos($nm, '[') !== false) { return '{' . $nm . '}'; }
        if (strpos($nm, '(') !== false) { return '[' . $nm . ']'; }
        return '(' . $nm . ')';
    }

    /**
     * Ensambla raíz + insaturación + sufijo aplicando la elisión de la vocal final
     * («e» en inglés, «o» en español) ante vocal (P-16.7). $sufijo ya incluye sus localizadores.
     */
    private function ensamblar(string $raiz, array $dob, array $trip, bool $locs, string $sufijo = ''): string
    {
        $p = $raiz;
        if (!$dob && !$trip) {
            $p .= 'an§';
        } else {
            $primero = $dob ?: $trip;
            if (count($primero) > 1) { $p .= 'a'; }
            if ($dob) {
                $p .= ($locs ? '-' . implode(',', $dob) . '-' : '') . (count($dob) > 1 ? self::MULT_SIMPLE[count($dob)] : '') . 'en§';
            }
            if ($trip) {
                $p .= ($locs ? '-' . implode(',', $trip) . '-' : '') . (count($trip) > 1 ? self::MULT_SIMPLE[count($trip)] : '') . ($this->es() ? 'in' : 'yn') . '§';
            }
        }
        $p .= $sufijo;
        return $this->elidir($p);
    }

    private function elidir(string $p): string
    {
        $vocal = $this->es() ? 'aeiouáéíóú' : 'aeiouy';
        $fin = $this->es() ? 'o' : 'e';
        return preg_replace_callback("/§((?:[\\-\\d,']|(?<=\\d)[a-z]|\\(\\d+[a-z]?'*H(?:,\\d+[a-z]?'*H)*\\))*)(.?)/u", function ($m) use ($vocal, $fin) {
            $sig = $m[2];
            if ($sig === '' || mb_strpos($vocal, mb_strtolower($sig)) === false) { return $fin . $m[1] . $sig; }
            return $m[1] . $sig;
        }, $p);
    }

    private static function cmpLocs(array $x, array $y): int
    {
        $n = min(count($x), count($y));
        for ($i = 0; $i < $n; $i++) {
            $a = self::v($x[$i]); $b = self::v($y[$i]);
            if ($a !== $b) { return $a <=> $b; }
        }
        return count($x) <=> count($y);
    }

    /**
     * Valor numérico de un localizador para compararlo (P-31.1.4): 4 < 4a < 5; 1 < 1' < 2;
     * localizador compuesto 1(6) > 1. Los localizadores N se excluyen (−1).
     */
    private static function v($l): float
    {
        if (is_int($l)) { return (float)$l; }
        if (preg_match("/^(\\d+)([a-z]?)('*)(?:\\((\\d+)\\))?$/", (string)$l, $mm)) {
            return (int)$mm[1] + ($mm[2] !== '' ? 0.5 : 0) + 0.01 * strlen($mm[3]) + (isset($mm[4]) ? 0.0001 * (int)$mm[4] : 0);
        }
        return -1.0;
    }

    /** Localizador de nitrógeno (N, N1, N'…). */
    private static function esN($l): bool { return is_string($l) && $l !== '' && $l[0] === 'N'; }

    private static function ordenarLocs(array $l): array
    {
        usort($l, fn($a, $b) => self::v($a) <=> self::v($b));
        return $l;
    }

    /* ======================= evaluación de una cadena o anillo ======================= */

    /**
     * Evalúa una secuencia ordenada de átomos (índice 0 → localizador 1).
     * Opciones: padre (bool), anillo (bool), excluir (set), fv (localizador de la valencia libre),
     *           acilo (átomo carbonílico cuyos heteroátomos pertenecen al grupo acilo).
     */
    private function evaluar(array $pos, array $op): array
    {
        $m = $this->mol;
        $padre = !empty($op['padre']);
        $anillo = !empty($op['anillo']);
        $excl = $op['excluir'] ?? [];
        $acilo = $op['acilo'] ?? -1;
        $en = array_flip($pos);
        $N = count($pos);
        $subs = []; $princ = []; $mult = []; $dob = []; $trip = []; $invalido = false;
        $locsOp = $op['locs'] ?? null;
        $aristasOp = $op['aristas'] ?? null;
        $ordenes = $op['ordenes'] ?? [];
        for ($i = 0; $i < $N; $i++) {
            $a = $pos[$i]; $L = $locsOp !== null ? $locsOp[$i] : $i + 1;
            if ($aristasOp === null && ($i < $N - 1 || ($anillo && $N > 2))) {
                $b = $pos[($i + 1) % $N];
                $o = $m->orden($a, $b);
                $locB = ($i === $N - 1) ? $N : $L;
                if ($o === 2) { $mult[] = $locB; $dob[] = $locB; }
                if ($o === 3) { $mult[] = $locB; $trip[] = $locB; }
            }
            $tA = $padre ? $this->tipoFC($a, true) : null;
            if ($tA !== null) {                         // carbono funcional incorporado a la cadena principal
                $princ[] = ['loc' => $L, 't' => $tA, 'fc' => $a];
            }
            foreach ($m->adj[$a] as $x => $o) {
                if (isset($en[$x]) || isset($excl[$x])) { continue; }
                if (($a === $acilo || $tA !== null) && $m->el[$x] !== 'C') { continue; }
                $ex = $m->el[$x];
                if ($ex === 'C') {
                    $tx = $this->tipoFC($x, $padre);
                    if ($tx !== null) {
                        if ($padre && $anillo && $tx === $this->principal) {
                            $princ[] = ['loc' => $L, 't' => $tx, 'fc' => $x];
                        } else {
                            $subs[] = ['loc' => $L, 'info' => $this->prefijoFC($x, $a), 'atomo' => $x, 'desde' => $a, 'orden' => 1];
                        }
                        continue;
                    }
                    if ($o === 3) { $invalido = true; continue; }
                    $ord = $o === 2 ? 2 : 1;
                    $subs[] = ['loc' => $L, 'info' => $this->sustituyente($x, $a, $ord), 'atomo' => $x, 'desde' => $a, 'orden' => $ord];
                } elseif (in_array($ex, MolOrg::HALOGENOS, true)) {
                    $subs[] = ['loc' => $L, 'info' => $this->info(self::HALO[$this->lang][$ex], 'simple', 0, ['hal' => true]), 'atomo' => -1];
                } elseif ($ex === 'N') {
                    if ($m->esNitroN($x)) {
                        $subs[] = ['loc' => $L, 'info' => $this->info('nitro', 'simple', 0, ['hal' => true]), 'atomo' => -1];
                    } elseif ($this->esNAmida($x)) {
                        $subs[] = ['loc' => $L, 'info' => $this->acilamino($x, $a), 'atomo' => -1];
                    } elseif ($padre && $this->principal === 'amina') {
                        $princ[] = ['loc' => $L, 't' => 'amina', 'n' => $x, 'desde' => $a];
                    } else {
                        $subs[] = ['loc' => $L, 'info' => $this->amino($x, $a), 'atomo' => -1];
                    }
                } elseif ($ex === 'O') {
                    if ($o === 2) {
                        if ($padre && $this->principal === 'cetona' && !isset($this->ald[$a])) {
                            $princ[] = ['loc' => $L, 't' => 'cetona'];
                        } else {
                            $subs[] = ['loc' => $L, 'info' => $this->info('oxo', 'simple', 0, ['hal' => true]), 'atomo' => -1];
                        }
                        continue;
                    }
                    $otros = array_values(array_filter($m->vecinos($x), fn($v) => $v !== $a));
                    if (!$otros) {
                        if ($padre && $this->principal === 'alcohol') { $princ[] = ['loc' => $L, 't' => 'alcohol']; }
                        else { $subs[] = ['loc' => $L, 'info' => $this->info($this->es() ? 'hidroxi' : 'hydroxy', 'simple', 0, ['hal' => true]), 'atomo' => -1]; }
                    } else {
                        $y = $otros[0];
                        if (isset($this->fc[$y]) && $this->fc[$y]['t'] === 'ester' && $this->fc[$y]['o'] === $x) {
                            $subs[] = ['loc' => $L, 'info' => $this->aciloxi($y, $x), 'atomo' => -1];
                        } else {
                            $subs[] = ['loc' => $L, 'info' => $this->alcoxi($y, $x), 'atomo' => -1];
                        }
                    }
                }
            }
        }
        // prefijos N de los sufijos amina/amida (N, o N1, N2… si hay varios)
        if ($padre) {
            $conN = array_values(array_filter($princ, fn($q) => $q['t'] === 'amina' || $q['t'] === 'amida'));
            // N sin prima para el nitrógeno con más sustituyentes (localizadores más bajos), luego por orden alfabético
            $claveN = function (array $q) use ($m): array {
                $n = $q['t'] === 'amina' ? $q['n'] : $this->fc[$q['fc']]['n'];
                $ex = $q['t'] === 'amina' ? $q['desde'] : $q['fc'];
                $nm = [];
                foreach ($m->adj[$n] as $c => $o) { if ($c !== $ex) { $nm[] = $this->sustituyente($c, $n, 1)['clave']; } }
                sort($nm);
                return [$q['loc'], -count($nm), implode(',', $nm)];
            };
            usort($conN, fn($x, $y) => $claveN($x) <=> $claveN($y));
            foreach ($conN as $q) {
                if ($q['t'] === 'amina') { $n = $q['n']; $ex = [$q['desde']]; }
                else { $n = $this->fc[$q['fc']]['n']; $ex = [$q['fc']]; }
                $loc = 'N';
                if (count($conN) > 1) {
                    $locsN = array_map(fn($z) => $z['loc'], $conN);
                    if (count(array_unique($locsN)) === count($locsN)) { $loc = 'N' . $q['loc']; }
                    else { $loc = 'N' . str_repeat("'", (int)array_search($q, $conN, true)); }
                }
                foreach ($m->adj[$n] as $c => $o) {
                    if (in_array($c, $ex, true)) { continue; }
                    $subs[] = ['loc' => $loc, 'info' => $this->sustituyente($c, $n, 1), 'atomo' => $c, 'desde' => $n, 'orden' => 1, 'nsub' => true];
                }
            }
        }
        if ($aristasOp !== null) {                  // enlaces múltiples de un sistema de anillos
            foreach ($aristasOp as [$i, $j]) {
                $x = $pos[$i]; $y = $pos[$j];
                $o = $ordenes[min($x, $y) . '-' . max($x, $y)] ?? $m->orden($x, $y);
                if ($o !== 2 && $o !== 3) { continue; }
                $li = $locsOp[$i]; $lj = $locsOp[$j];
                if (self::v($li) > self::v($lj)) { [$li, $lj] = [$lj, $li]; }
                $lb = (is_int($li) && is_int($lj) && $lj - $li === 1) ? $li : $li . '(' . $lj . ')';
                $mult[] = $lb;
                if ($o === 2) { $dob[] = $lb; } else { $trip[] = $lb; }
            }
        }
        $mult = self::ordenarLocs($mult); $dob = self::ordenarLocs($dob); $trip = self::ordenarLocs($trip);
        $pl = self::ordenarLocs(array_map(fn($p) => $p['loc'], $princ));
        return ['pos' => $pos, 'subs' => $subs, 'princ' => $princ, 'plocs' => $pl, 'mlocs' => $mult, 'dlocs' => $dob,
            'tlocs' => $trip, 'mn' => count($mult), 'dn' => count($dob), 'anillo' => $anillo, 'len' => $N, 'invalido' => $invalido];
    }

    /** Agrupa y ordena los prefijos (P-14.5, P-16.3); calcula los conjuntos de localizadores. */
    private function armar(array $subs, bool $omitirLoc): array
    {
        $grupos = [];
        foreach ($subs as $s) {
            $nm = $s['info']['nombre'];
            if (!isset($grupos[$nm])) { $grupos[$nm] = ['info' => $s['info'], 'locs' => []]; }
            $grupos[$nm]['locs'][] = $s['loc'];
        }
        uasort($grupos, function ($a, $b) {
            $c = strcmp($a['info']['clave'], $b['info']['clave']);
            return $c !== 0 ? $c : strcmp($a['info']['nombre'], $b['info']['nombre']);
        });
        $partes = []; $citados = []; $locs = [];
        foreach ($grupos as $nm => $g) {
            usort($g['locs'], function ($x, $y) {
                $xs = self::esN($x); $ys = self::esN($y);
                if ($xs !== $ys) { return $xs ? -1 : 1; }
                return $xs ? strnatcmp($x, $y) : self::v($x) <=> self::v($y);
            });
            $cnt = count($g['locs']);
            $info = $g['info'];
            foreach ($g['locs'] as $l) { if (!self::esN($l)) { $citados[] = $l; $locs[] = $l; } }
            $paren = $info['tipo'] !== 'simple';
            if ($cnt === 1) { $mult = ''; }
            elseif ($info['tipo'] === 'complejo') { $mult = self::MULT_COMPLEJO[$cnt] ?? ('(' . $cnt . ')'); }
            else { $mult = self::MULT_SIMPLE[$cnt]; }
            if ($mult !== '' && !$paren && preg_match('/^(sec|terc|tert)-/', $nm)) { $mult .= '-'; }
            $txt = $mult . ($paren ? self::encerrar($nm) : $nm);
            $locTxt = array_filter($g['locs'], fn($l) => self::esN($l) || !$omitirLoc);
            $partes[] = [$locTxt ? implode(',', $locTxt) . '-' : '', $txt];
        }
        $locs = self::ordenarLocs($locs);
        $pref = '';
        $atomicos = '/^(di|tri|tetra|penta|hexa)?(fluoro|chloro|cloro|bromo|iodo|yodo|nitro|oxo|hydroxy|hidroxi|cyano|ciano|carboxy|carboxi|formyl|formil|amino)$/';
        foreach ($partes as [$lt, $txt]) {
            // sin localizadores los prefijos se escriben seguidos (bromoclorometano); un prefijo que admite
            // sustitución se encierra entre paréntesis para evitar ambigüedad: bromo(etoxi)metil, cloro(fenil)acético
            if ($pref !== '' && $lt === '' && !preg_match($atomicos, $txt) && $txt[0] !== '(' && $txt[0] !== '[' && $txt[0] !== '{') {
                $txt = self::encerrar($txt);
            }
            $pref .= ($pref !== '' && $lt !== '' ? '-' : '') . $lt . $txt;
        }
        return ['prefijos' => $pref, 'locs' => $locs, 'citados' => $citados, 'grupos' => $grupos,
            'cuenta' => count(array_filter($subs, fn($s) => !self::esN($s['loc'])))];
    }

    /** Compara dos candidatos. Devuelve [signo, criterio]. */
    private function comparar(array $a, array $b, ?array $seq = null): array
    {
        $seq = $seq ?? $this->secuencia();
        foreach ($seq as $k) {
            switch ($k) {
                case 'pn': $c = count($b['princ']) <=> count($a['princ']); break;
                case 'anillo': $c = (int)$b['anillo'] <=> (int)$a['anillo']; break;
                case 'len': case 'mn': case 'dn': case 'cuenta': case 'nr': case 'insat': $c = ($b[$k] ?? 0) <=> ($a[$k] ?? 0); break;
                case 'plocs': case 'mlocs': case 'dlocs': case 'locs': case 'citados':
                case 'ilocs': case 'alocs': case 'hlocs': case 'union':
                    $c = self::cmpLocs($a[$k] ?? [], $b[$k] ?? []); break;
                case 'nombre': $c = strcmp($a['nombreCmp'], $b['nombreCmp']); $c = $c <=> 0; break;
                default: $c = 0;
            }
            if ($c !== 0) { return [$c, $k]; }
        }
        return [0, ''];
    }

    /**
     * Orden de los criterios: selección de la estructura principal (P-44.1, P-44.2, P-44.3, P-44.4)
     * y numeración (P-31.1.4): hidrógeno indicado → grupo principal/valencia libre → punto de unión
     * (multiplicativa) → hidrógeno añadido → hidro → enlaces múltiples → prefijos.
     */
    private function secuencia(): array
    {
        $sel = $this->estilo === 'pin' ? ['pn', 'anillo', 'nr', 'len', 'insat', 'mn', 'dn'] : ['pn', 'anillo', 'nr', 'mn', 'dn', 'len', 'insat'];
        return array_merge($sel, ['ilocs', 'plocs', 'union', 'alocs', 'hlocs', 'mlocs', 'dlocs', 'cuenta', 'locs', 'citados', 'nombre']);
    }

    /* ======================= sustituyentes ======================= */

    /**
     * Prefijo del sustituyente que nace en el carbono $r unido a $p (enlace sencillo: -il;
     * doble: -ilideno). Devuelve info: nombre, tipo (simple|compuesto|complejo), clave.
     */
    public function sustituyente(int $r, int $p, int $ord = 1): array
    {
        $k = "$r|$p|$ord";
        if (isset($this->memo[$k])) { return $this->memo[$k]; }
        if (isset($this->fc[$r])) { return $this->memo[$k] = $this->prefijoFC($r, $p); }
        if (isset($this->anilloDe[$r])) { return $this->memo[$k] = $this->sustSistema($r, $p, $ord); }
        $m = $this->mol;
        if ($ord === 1) {
            foreach ($m->adj[$r] as $x => $o) {
                if ($m->el[$x] === 'O' && $o === 2) { return $this->memo[$k] = $this->acilo($r, $p); }
            }
        }
        return $this->memo[$k] = $this->sustCadena($r, $p, $ord);
    }

    /** Átomos acíclicos del esqueleto alcanzables desde $r sin pasar por $p. */
    private function subarbolAciclico(int $r, int $p, bool $padre = false): array
    {
        $vis = [$r => true, $p => true];
        $pila = [$r]; $res = [];
        while ($pila) {
            $a = array_pop($pila); $res[] = $a;
            foreach ($this->mol->adj[$a] as $b => $o) {
                if (isset($vis[$b])) { continue; }
                if (!$this->esEsqueleto($b, $padre) || isset($this->anilloDe[$b])) { continue; }
                $vis[$b] = true; $pila[] = $b;
            }
        }
        return $res;
    }

    /** Camino único entre u y v dentro de un conjunto acíclico permitido. */
    private function camino(int $u, int $v, array $perm): array
    {
        $padre = [$u => -1]; $cola = [$u];
        while ($cola) {
            $a = array_shift($cola);
            if ($a === $v) { break; }
            foreach ($this->mol->adj[$a] as $b => $o) {
                if (isset($perm[$b]) && !array_key_exists($b, $padre)) { $padre[$b] = $a; $cola[] = $b; }
            }
        }
        if (!array_key_exists($v, $padre)) { return []; }
        $ruta = [];
        for ($x = $v; $x !== -1; $x = $padre[$x]) { $ruta[] = $x; }
        return array_reverse($ruta);
    }

    private function sustCadena(int $r, int $p, int $ord): array
    {
        $atomos = $this->subarbolAciclico($r, $p);
        $perm = array_fill_keys($atomos, true);
        $ext = [$r];
        foreach ($atomos as $a) {
            if ($a === $r) { continue; }
            $g = 0;
            foreach ($this->mol->adj[$a] as $b => $o) { if (isset($perm[$b])) { $g++; } }
            if ($g <= 1) { $ext[] = $a; }
        }
        $desdeUnion = $this->estilo !== 'pin';
        $cands = [];
        foreach ($ext as $u) {
            foreach ($ext as $v) {
                if ($desdeUnion && $u !== $r) { continue; }
                if ($u === $v && $u !== $r) { continue; }
                $ruta = $u === $v ? [$u] : $this->camino($u, $v, $perm);
                $fv = array_search($r, $ruta, true);
                if ($fv === false) { continue; }
                $d = $this->evaluar($ruta, ['excluir' => [$p => true]]);
                if ($d['invalido']) { continue; }
                $d['plocs'] = [$fv + 1];
                $d['princ'] = [];
                $d = $this->completar($d, false);
                $d['nombreCmp'] = $d['prefijos'];
                $cands[] = $d;
            }
        }
        if (!$cands) { throw new NomenclaturaException('Sustituyente unido por un triple enlace: no admitido.'); }
        usort($cands, fn($a, $b) => $this->comparar($a, $b)[0]);
        $best = $cands[0];
        $n = $best['len'];
        $fv = $best['plocs'][0];
        if ($n === 1) { $best = $this->completar($best, true); }
        $yl = $this->yl();
        $suf = $ord === 2 ? ($this->es() ? 'iliden' : 'ylidene') : $yl;
        $insat = $best['dlocs'] || $best['tlocs'];
        $raiz = $this->raiz($n);
        if (!$insat && $fv === 1) {
            $base = $raiz . $suf;
        } elseif ($n <= 2 && !($n === 2 && $fv === 2)) {
            $base = $this->ensamblar($raiz, $best['dlocs'], $best['tlocs'], false, $suf);
        } elseif ($fv === 1 && $desdeUnion) {
            $base = $this->ensamblar($raiz, $best['dlocs'], $best['tlocs'], true, $suf);
        } else {
            $base = $this->ensamblar($raiz, $best['dlocs'], $best['tlocs'], true, '-' . $fv . '-' . $suf);
        }
        $nombre = $best['prefijos'] . $base;

        // prefijos retenidos
        $en = !$this->es();
        if ($this->estilo !== 'sis79') {
            $tb = $en ? ['2-methylpropan-2-yl' => 'tert-butyl', '1,1-dimethylethyl' => 'tert-butyl']
                      : ['2-metilpropan-2-il' => 'terc-butil', '1,1-dimetiletil' => 'terc-butil'];
            if (isset($tb[$nombre])) { return $this->info($tb[$nombre], 'simple', 3); }
            // bencilo (P-29.6): C6H5-CH2- sin sustituir
            if ($n === 1 && count($best['subs']) === 1 && $best['subs'][0]['info']['nombre'] === ($en ? 'phenyl' : 'fenil')) {
                return $this->info(($en ? 'benzyl' : 'bencil') . ($ord === 2 ? ($en ? 'idene' : 'iden') : ''), 'simple', 1);
            }
        }
        if ($this->estilo === 'trad') {
            $ret = $en
                ? ['1-methylethyl' => 'isopropyl', '1-methylpropyl' => 'sec-butyl', '2-methylpropyl' => 'isobutyl',
                   '3-methylbutyl' => 'isopentyl', '2,2-dimethylpropyl' => 'neopentyl', '1,1-dimethylpropyl' => 'tert-pentyl',
                   'ethenyl' => 'vinyl', 'prop-2-enyl' => 'allyl', '1-methylethenyl' => 'isopropenyl',
                   '1-methylethylidene' => 'isopropylidene', 'methylidene' => 'methylene']
                : ['1-metiletil' => 'isopropil', '1-metilpropil' => 'sec-butil', '2-metilpropil' => 'isobutil',
                   '3-metilbutil' => 'isopentil', '2,2-dimetilpropil' => 'neopentil', '1,1-dimetilpropil' => 'terc-pentil',
                   'etenil' => 'vinil', 'prop-2-enil' => 'alil', '1-metiletenil' => 'isopropenil',
                   '1-metiletiliden' => 'isopropiliden', 'metiliden' => 'metilen'];
            if (isset($ret[$nombre])) { return $this->info($ret[$nombre], 'simple', $n); }
        }
        if ($this->estilo === 'sis79' && $nombre === ($en ? 'methylidene' : 'metiliden')) { $nombre = $en ? 'methylene' : 'metilen'; }
        if ($best['cuenta'] > 0) { return $this->info($nombre, 'complejo', $n); }
        return $this->info($nombre, preg_match('/\d/', $base) ? 'compuesto' : 'simple', $n);
    }

    /* ======================= sistemas de anillos (fase 2) ======================= */

    /** Hidruros fusionados con nombre retenido (P-25.1.1): tamaños → [inglés, español, perímetro desde C1, átomos de fusión]. */
    private const FUSIONADOS = [
        '5,5' => ['pentalen', 'pentalen', [1, 2, 3, '3a', 4, 5, 6, '6a'], ['3a', '6a']],
        '5,6' => ['inden', 'inden', [1, 2, 3, '3a', 4, 5, 6, 7, '7a'], ['3a', '7a']],
        '5,7' => ['azulen', 'azulen', [1, 2, 3, '3a', 4, 5, 6, 7, 8, '8a'], ['3a', '8a']],
        '6,6' => ['naphthalen', 'naftalen', [1, 2, 3, 4, '4a', 5, 6, 7, 8, '8a'], ['4a', '8a']],
        '7,7' => ['heptalen', 'heptalen', [1, 2, 3, 4, 5, '5a', 6, 7, 8, 9, 10, '10a'], ['5a', '10a']],
    ];

    /**
     * Datos de nomenclatura de un sistema: clase (mono|fus|vb|espiro), numeraciones válidas
     * (pos, locs, aristas como pares de índices), plantilla o descriptor entre corchetes.
     */
    private function datosSistema(int $k): array
    {
        if (isset($this->infoSis[$k])) { return $this->infoSis[$k]; }
        $S = $this->sistemas[$k];
        $m = $this->mol;
        $nums = [];
        $info = ['n' => count($S['atomos']), 'atomos' => $S['atomos']];
        $idxAristas = function (array $pos) use ($S): array {
            $ix = array_flip($pos); $r = [];
            foreach ($S['aristas'] as [$x, $y]) { $r[] = [$ix[$x], $ix[$y]]; }
            return $r;
        };
        if ($S['tipo'] === 'mono') {
            $info['clase'] = 'mono';
            $c = $S['ciclo']; $N = count($c);
            $info['arom'] = $N === 6 && $m->arom[$c[0]] && count(array_filter($c, fn($a) => $m->arom[$a])) === 6;
            for ($s0 = 0; $s0 < $N; $s0++) {
                foreach ([1, -1] as $dir) {
                    $pos = [];
                    for ($i = 0; $i < $N; $i++) { $pos[] = $c[(($s0 + $dir * $i) % $N + $N) % $N]; }
                    $nums[] = ['pos' => $pos, 'locs' => range(1, $N), 'aristas' => $idxAristas($pos)];
                }
            }
        } elseif ($S['tipo'] === 'espiro') {
            $info['clase'] = 'espiro';
            [$c1, $c2] = $S['ciclos'];
            $x = count($c1) - 1; $y = count($c2) - 1;
            $info['corchete'] = '[' . min($x, $y) . '.' . max($x, $y) . ']';
            $ordenes = $x === $y ? [[$c1, $c2], [$c2, $c1]] : ($x < $y ? [[$c1, $c2]] : [[$c2, $c1]]);
            foreach ($ordenes as [$pe, $gr]) {
                $restoP = array_slice($pe, 1); $restoG = array_slice($gr, 1);
                foreach ([$restoP, array_reverse($restoP)] as $rp) {
                    foreach ([$restoG, array_reverse($restoG)] as $rg) {
                        $pos = array_merge($rp, [$pe[0]], $rg);
                        $nums[] = ['pos' => $pos, 'locs' => range(1, count($pos)), 'aristas' => $idxAristas($pos)];
                    }
                }
            }
        } else {
            [$h1, $h2] = $S['cabezas'];
            $P = $S['puentes'];
            $tam = array_map('count', $P);
            $cero = array_search(0, $tam, true);
            $anillosT = [];
            if ($cero !== false) {
                $otros = array_values(array_diff_key($P, [$cero => 1]));
                $anillosT = [count($otros[0]) + 2, count($otros[1]) + 2];
                sort($anillosT);
            }
            $clave = implode(',', $anillosT);
            if ($cero !== false && $anillosT[0] >= 5 && isset(self::FUSIONADOS[$clave])) {
                // ---- sistema orto-fusionado con nombre retenido (P-25.1) ----
                $info['clase'] = 'fus';
                $info['plantilla'] = $clave;
                [$en, $es, $T, $F] = self::FUSIONADOS[$clave];
                $otros = array_values(array_diff_key($P, [$cero => 1]));
                $per = array_merge([$h1], $otros[0], [$h2], array_reverse($otros[1]));
                $n = count($per);
                $iF0 = array_search($F[0], $T, true); $iF1 = array_search($F[1], $T, true);
                for ($s0 = 0; $s0 < $n; $s0++) {
                    foreach ([1, -1] as $dir) {
                        $pos = [];
                        for ($i = 0; $i < $n; $i++) { $pos[] = $per[(($s0 + $dir * $i) % $n + $n) % $n]; }
                        $f = [$pos[$iF0], $pos[$iF1]]; sort($f); $h = [$h1, $h2]; sort($h);
                        if ($f !== $h) { continue; }
                        // la plantilla también fija el tamaño de cada anillo
                        $ok = true;
                        $ix = array_flip($pos);
                        foreach ($S['aristas'] as [$x, $y]) {
                            $d = abs($ix[$x] - $ix[$y]);
                            if (!($d === 1 || $d === $n - 1 || ([$ix[$x], $ix[$y]] == [$iF0, $iF1]) || ([$ix[$y], $ix[$x]] == [$iF0, $iF1]))) { $ok = false; }
                        }
                        if ($ok) { $nums[] = ['pos' => $pos, 'locs' => $T, 'aristas' => $idxAristas($pos)]; }
                    }
                }
                // posiciones válidas de hidrógeno indicado del hidruro progenitor (número impar de átomos)
                $info['ih'] = [];
                if ($n % 2 === 1) {
                    foreach ($S['atomos'] as $a) {
                        $resto = array_values(array_diff($S['atomos'], [$a]));
                        if (MolOrg::emparejamientos($resto, $S['aristas'], 1)) { $info['ih'][] = $a; }
                    }
                }
            } elseif ($cero !== false && $anillosT[0] >= 5) {
                throw new NomenclaturaException('Sistema fusionado de anillos de ' . $anillosT[0] . ' y ' . $anillosT[1] . ' miembros: su nombre de fusión (p. ej., benzo[7]anuleno) no es un nombre retenido y se incorporará en una fase posterior.');
            } else {
                // ---- von Baeyer (P-23.2) ----
                $info['clase'] = 'vb';
                $ord = $tam; rsort($ord);
                $info['corchete'] = '[' . implode('.', $ord) . ']';
                foreach ([[$h1, $h2, $P], [$h2, $h1, array_map('array_reverse', $P)]] as [$c1, $c2, $PP]) {
                    foreach ([[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] as [$i, $j, $l]) {
                        if (!(count($PP[$i]) >= count($PP[$j]) && count($PP[$j]) >= count($PP[$l]))) { continue; }
                        $pos = array_merge([$c1], $PP[$i], [$c2], array_reverse($PP[$j]), $PP[$l]);
                        $nums[] = ['pos' => $pos, 'locs' => range(1, count($pos)), 'aristas' => $idxAristas($pos)];
                    }
                }
            }
        }
        // átomos con insaturación endocíclica (para la antigüedad de anillos P-44.4.1 y el patrón «hidro»)
        $ins = [];
        foreach ($S['aristas'] as [$x, $y]) {
            $o = $m->orden($x, $y);
            if ($o === 2 || $o === 3 || $o === MolOrg::AROM) { $ins[$x] = true; $ins[$y] = true; }
        }
        $info['insat'] = $ins;
        $info['nums'] = $nums;
        return $this->infoSis[$k] = $info;
    }

    /** Descriptor que identifica sistemas idénticos (ensamblajes y nomenclatura multiplicativa). */
    private function firmaSistema(int $k): ?string
    {
        $d = $this->datosSistema($k);
        $nIns = count($d['insat']);
        if ($d['clase'] === 'mono') {
            if (!empty($d['arom'])) { return 'benceno'; }
            return $nIns === 0 ? 'ciclo' . $d['n'] : null;
        }
        if ($d['clase'] === 'fus') {
            if ($nIns === $d['n'] - ($d['ih'] ? 1 : 0)) { return 'fus' . $d['plantilla']; }
            return $nIns === 0 ? 'fusat' . $d['plantilla'] : null;
        }
        return $nIns === 0 ? $d['clase'] . ($d['corchete'] ?? '') : null;
    }

    /**
     * Candidatos (numeración × hidrógeno indicado/añadido/hidro × Kekulé) de un sistema de anillos.
     * $op: padre, excluir, fv (átomo con valencia libre), ordFv (1|2), union (átomo unido al conector multiplicativo).
     */
    private function candidatosSistema(int $k, array $op): array
    {
        $m = $this->mol;
        $D = $this->datosSistema($k);
        $S = $this->sistemas[$k];
        $padre = !empty($op['padre']);
        $fv = $op['fv'] ?? null;
        $ordFv = $op['ordFv'] ?? 1;
        $union = $op['union'] ?? null;
        $cands = [];
        // estructuras de Kekulé alternativas de los anillos «bencénicos» de sistemas von Baeyer o espiro
        // (benzociclobuteno: octa-1,3,5-trieno frente a octa-1(6),2,4-trieno)
        $variantes = [[]];
        if (in_array($D['clase'], ['vb', 'espiro'], true)) {
            $U = []; $EU = [];
            foreach ($m->anillos() as $r) {
                if (count($r) !== 6 || array_diff($r, $S['atomos'])) { continue; }
                $dob = 0;
                for ($i = 0; $i < 6; $i++) { if ($m->orden($r[$i], $r[($i + 1) % 6]) === 2) { $dob++; } }
                if ($dob !== 3) { continue; }
                foreach ($r as $a) { $U[$a] = true; }
                for ($i = 0; $i < 6; $i++) { $EU[] = [$r[$i], $r[($i + 1) % 6]]; }
            }
            if ($U) {
                $variantes = [];
                foreach (MolOrg::emparejamientos(array_keys($U), $EU, 16) as $mt) {
                    $ov = [];
                    foreach ($EU as [$x, $y]) { $ov[min($x, $y) . '-' . max($x, $y)] = 1; }
                    foreach ($mt as [$x, $y]) { $ov[min($x, $y) . '-' . max($x, $y)] = 2; }
                    $variantes[] = $ov;
                }
            }
        }
        $sat = array_values(array_diff($S['atomos'], array_keys($D['insat'])));
        // átomos que deben ser sp3 en el progenitor por el sufijo (-ona, -ilideno): hidrógeno indicado o añadido (P-14.7)
        $K = [];
        if ($D['clase'] === 'fus') {
            foreach ($S['atomos'] as $a) {
                if ($padre && $this->principal === 'cetona') {
                    foreach ($m->adj[$a] as $x => $o) { if ($m->el[$x] === 'O' && $o === 2) { $K[] = $a; } }
                }
            }
            if ($fv !== null && $ordFv === 2) { $K[] = $fv; }
        }
        foreach ($D['nums'] as $num) {
            foreach ($variantes as $ov) {
                $d = $this->evaluar($num['pos'], ['padre' => $padre, 'anillo' => true, 'excluir' => $op['excluir'] ?? [],
                    'locs' => $num['locs'], 'aristas' => $num['aristas'], 'ordenes' => $ov]);
                if ($d['invalido']) { continue; }
                $lab = array_combine($num['pos'], $num['locs']);
                $d['sis'] = $k; $d['clase'] = $D['clase']; $d['lab'] = $lab;
                $d['nr'] = $D['clase'] === 'mono' ? 1 : 2;
                $d['len'] = $D['n'];
                $d['insat'] = count($D['insat']);
                $d['arom'] = !empty($D['arom']);
                if ($d['arom'] || $D['clase'] === 'fus') { $d['mlocs'] = $d['dlocs'] = $d['tlocs'] = []; $d['mn'] = $d['dn'] = 0; }
                if ($fv !== null) { $d['plocs'] = [$lab[$fv]]; $d['princ'] = []; }
                if ($union !== null) { $d['union'] = [$lab[$union]]; }
                if ($D['clase'] !== 'fus') { $cands[] = $d; continue; }
                // ---- hidrógeno indicado (I), añadido (A) e hidro (H) ----
                $Iops = $D['ih'] ? array_map(fn($a) => [$a], array_values(array_intersect($D['ih'], $sat))) : [[]];
                foreach ($Iops as $I) {
                    $Kr = array_values(array_diff(array_unique($K), $I));
                    $libres = array_values(array_diff($sat, $I, $Kr));
                    $mejorTam = null; $As = [];
                    foreach ([0, 1, 2] as $t) {
                        if ($mejorTam !== null) { break; }
                        foreach (self::combinaciones($libres, $t) as $A) {
                            $resto = array_values(array_diff($S['atomos'], $I, $Kr, $A));
                            if (count($resto) % 2 === 0 && MolOrg::emparejamientos($resto, $S['aristas'], 1)) { $As[] = $A; $mejorTam = $t; }
                        }
                    }
                    foreach ($As as $A) {
                        $H = array_values(array_diff($sat, $I, $Kr, $A));
                        $c = $d;
                        $c['ilocs'] = self::ordenarLocs(array_map(fn($a) => $lab[$a], $I));
                        $c['alocs'] = self::ordenarLocs(array_map(fn($a) => $lab[$a], $A));
                        $c['hlocs'] = self::ordenarLocs(array_map(fn($a) => $lab[$a], $H));
                        $cands[] = $c;
                    }
                }
            }
        }
        return $cands;
    }

    private static function combinaciones(array $l, int $t): array
    {
        if ($t === 0) { return [[]]; }
        $r = [];
        foreach ($l as $i => $x) {
            foreach (self::combinaciones(array_slice($l, $i + 1), $t - 1) as $c) { $r[] = array_merge([$x], $c); }
        }
        return $r;
    }

    /** Une partes de un nombre con guion cuando la siguiente empieza por localizador. */
    private static function unir(array $partes): string
    {
        $r = '';
        foreach ($partes as $p) {
            if ($p === '') { continue; }
            $r .= ($r !== '' && preg_match('/^\d/', $p) ? '-' : '') . $p;
        }
        return $r;
    }

    /**
     * Nombre del hidruro de un sistema con su sufijo (sin prefijos sustituyentes):
     * devuelve [hidro, base]; p. ej. ['2,3-dihydro', '1H-inden-1-one'].
     * $sufijo ya incluye los localizadores y el hidrógeno añadido.
     */
    private function nombreSistema(array $d, string $sufijo, bool $citarInsat = true): array
    {
        $en = !$this->es();
        $D = $this->datosSistema($d['sis']);
        $n = $D['n'];
        switch ($D['clase']) {
            case 'mono':
                if (!empty($D['arom'])) { return ['', $this->elidir(($en ? 'benzen§' : 'bencen§') . $sufijo)]; }
                return ['', $this->ensamblar(($en ? 'cyclo' : 'ciclo') . $this->raiz($n), $d['dlocs'], $d['tlocs'], $citarInsat, $sufijo)];
            case 'fus':
                $stem = self::FUSIONADOS[$D['plantilla']][$en ? 0 : 1];
                $ih = $d['ilocs'] ? implode(',', $d['ilocs']) . 'H-' : '';
                $hid = '';
                if ($d['hlocs']) {
                    $mult = self::MULT_SIMPLE[count($d['hlocs'])] . ($en ? 'hydro' : 'hidro');
                    $total = count($D['insat']) === 0 && !$d['ilocs'];
                    $hid = $total ? $mult : implode(',', $d['hlocs']) . '-' . $mult;
                }
                return [$hid, $ih . $this->elidir($stem . '§' . $sufijo)];
            case 'vb':
                return ['', $this->ensamblar(($en ? 'bicyclo' : 'biciclo') . $D['corchete'] . $this->raiz($n), $d['dlocs'], $d['tlocs'], true, $sufijo)];
            case 'espiro':
                return ['', $this->ensamblar(($en ? 'spiro' : 'espiro') . $D['corchete'] . $this->raiz($n), $d['dlocs'], $d['tlocs'], true, $sufijo)];
        }
        throw new NomenclaturaException('Sistema de anillos no reconocido.');
    }

    /** Localizadores de un sufijo con hidrógeno añadido: «-1(2H)-». */
    private static function locsSufijo(array $locs, array $alocs): string
    {
        $t = implode(',', $locs);
        if ($alocs) { $t .= '(' . implode(',', array_map(fn($l) => $l . 'H', $alocs)) . ')'; }
        return '-' . $t . '-';
    }

    /** Anillo o sistema de anillos como sustituyente: fenil, ciclohexil, naftalen-1-il, biciclo[2.2.1]heptan-2-il… */
    private function sustSistema(int $r, int $p, int $ord): array
    {
        $k = $this->anilloDe[$r];
        $D = $this->datosSistema($k);
        $cands = $this->candidatosSistema($k, ['excluir' => [$p => true], 'fv' => $r, 'ordFv' => $ord]);
        foreach ($cands as &$c) { $c = $this->completar($c, false); $c['nombreCmp'] = $c['prefijos']; }
        unset($c);
        usort($cands, fn($a, $b) => $this->comparar($a, $b)[0]);
        $best = $cands[0];
        $en = !$this->es();
        $N = $D['n'];
        if ($D['clase'] === 'mono') {
            if (!empty($D['arom'])) {
                if ($ord === 2) { throw new NomenclaturaException('Anillo aromático unido por doble enlace: no admitido.'); }
                $base = $en ? 'phenyl' : 'fenil';
            } else {
                $raiz = ($en ? 'cyclo' : 'ciclo') . $this->raiz($N);
                $suf = $ord === 2 ? ($en ? 'ylidene' : 'iliden') : $this->yl();
                $insat = (bool)$best['dlocs'] || (bool)$best['tlocs'];
                $base = $insat ? $this->ensamblar($raiz, $best['dlocs'], $best['tlocs'], true, '-1-' . $suf) : $raiz . $suf;
            }
            $nombre = $best['prefijos'] . $base;
        } else {
            $suf = $ord === 2 ? ($en ? 'ylidene' : 'iliden') : $this->yl();
            [$hid, $base] = $this->nombreSistema($best, self::locsSufijo($best['plocs'], $best['alocs'] ?? []) . $suf);
            $nombre = self::unir([$best['prefijos'], $hid, $base]);
            if ($this->estilo === 'trad' && $best['cuenta'] === 0 && $hid === '' && $ord === 1) {
                $tr = $en ? ['naphthalen-1-yl' => '1-naphthyl', 'naphthalen-2-yl' => '2-naphthyl']
                          : ['naftalen-1-il' => '1-naftil', 'naftalen-2-il' => '2-naftil'];
                $nombre = $tr[$nombre] ?? $nombre;
            }
            $base = $nombre;
        }
        if ($best['cuenta'] > 0) { return $this->info($nombre, 'complejo', $N, ['anillo' => true]); }
        return $this->info($nombre, preg_match('/\d/', $base) ? 'compuesto' : 'simple', $N, ['anillo' => true]);
    }

    /** Grupo acilo R-CO- a partir del carbono carbonílico $f (acetil, propanoil, benzoil, formil…). */
    private function acilo(int $f, int $p): array
    {
        $m = $this->mol;
        $en = !$this->es();
        $g = -1;
        foreach ($m->adj[$f] as $x => $o) {
            if ($x !== $p && $m->el[$x] === 'C') { $g = $x; }
        }
        if ($g === -1) { return $this->info($en ? 'formyl' : 'formil', 'simple', 1); }
        if (isset($this->fc[$g])) { throw new NomenclaturaException('Dos grupos carbonilo contiguos en un sustituyente (oxalilo, glioxililo…): no se incluye en esta fase.'); }
        if (isset($this->anilloDe[$g])) {
            $k = $this->anilloDe[$g];
            $D = $this->datosSistema($k);
            $cands = $this->candidatosSistema($k, ['excluir' => [$f => true], 'fv' => $g]);
            foreach ($cands as &$c) { $c = $this->completar($c, false); $c['nombreCmp'] = $c['prefijos']; }
            unset($c);
            usort($cands, fn($a, $b) => $this->comparar($a, $b)[0]);
            $best = $cands[0];
            if ($D['clase'] === 'mono' && !empty($D['arom'])) {
                $nombre = $best['prefijos'] . ($en ? 'benzoyl' : 'benzoil');
            } elseif ($D['clase'] === 'mono') {
                $raiz = ($en ? 'cyclo' : 'ciclo') . $this->raiz($D['n']);
                $citar = $best['cuenta'] > 0 || $best['dlocs'] || $best['tlocs'];
                $nombre = $best['prefijos'] . $this->ensamblar($raiz, $best['dlocs'], $best['tlocs'], true, ($citar ? '-1-' : '') . ($en ? 'carbonyl' : 'carbonil'));
            } else {
                [$hid, $base] = $this->nombreSistema($best, self::locsSufijo($best['plocs'], []) . ($en ? 'carbonyl' : 'carbonil'));
                $nombre = self::unir([$best['prefijos'], $hid, $base]);
            }
            return $this->info($nombre, $best['cuenta'] > 0 ? 'complejo' : (preg_match('/\d/', $nombre) ? 'compuesto' : 'simple'), $D['n'] + 1);
        }
        // cadena: [f, g, ...]
        $atomos = $this->subarbolAciclico($g, $f);
        $perm = array_fill_keys($atomos, true);
        $ext = [$g];
        foreach ($atomos as $a) {
            if ($a === $g) { continue; }
            $n = 0;
            foreach ($m->adj[$a] as $b => $o) { if (isset($perm[$b])) { $n++; } }
            if ($n <= 1) { $ext[] = $a; }
        }
        $cands = [];
        foreach ($ext as $v) {
            $ruta = $v === $g ? [$g] : $this->camino($g, $v, $perm);
            $pos = array_merge([$f], $ruta);
            $d = $this->evaluar($pos, ['excluir' => [$p => true], 'acilo' => $f]);
            if ($d['invalido']) { continue; }
            $d['plocs'] = [1]; $d['princ'] = [];
            $d = $this->completar($d, false);
            $d['nombreCmp'] = $d['prefijos'];
            $cands[] = $d;
        }
        usort($cands, fn($a, $b) => $this->comparar($a, $b)[0]);
        $best = $cands[0];
        $n = $best['len'];
        if ($n === 2) {
            $s = $this->armar($best['subs'], true);
            $nombre = $s['prefijos'] . ($en ? 'acetyl' : 'acetil');
            return $this->info($nombre, $best['cuenta'] > 0 ? 'complejo' : 'simple', 2);
        }
        $nombre = $best['prefijos'] . $this->ensamblar($this->raiz($n), $best['dlocs'], $best['tlocs'], true, $en ? 'oyl' : 'oil');
        return $this->info($nombre, $best['cuenta'] > 0 ? 'complejo' : (preg_match('/\d/', $nombre) ? 'compuesto' : 'simple'), $n);
    }

    /** Prefijo de un carbono funcional no principal: carboxi, ciano, carbamoil, alcoxicarbonil. */
    private function prefijoFC(int $x, int $desde): array
    {
        $en = !$this->es();
        $d = $this->fc[$x] ?? null;
        if ($d === null) {                        // aldehído principal fuera de la estructura principal
            return $this->info($en ? 'formyl' : 'formil', 'simple', 1);
        }
        switch ($d['t']) {
            case 'acido': return $this->info($en ? 'carboxy' : 'carboxi', 'simple', 1);
            case 'nitrilo': return $this->info($en ? 'cyano' : 'ciano', 'simple', 1);
            case 'amida':
                $n = $this->nSustituido($d['n'], [$x], $en ? 'carbamoyl' : 'carbamoil');
                return $n;
            case 'ester':
                $ao = $this->alcoxi($d['r'], $d['o']);
                $nm = ($ao['tipo'] === 'simple' ? $ao['nombre'] : self::encerrar($ao['nombre'])) . ($en ? 'carbonyl' : 'carbonil');
                return $this->info($nm, 'complejo', 1);
        }
        throw new NomenclaturaException('Grupo no reconocido.');
    }

    /** Prefijos amino / carbamoil con sustituyentes en el N: metilamino, dimetilamino, etil(metil)amino… */
    private function nSustituido(int $n, array $excluir, string $base): array
    {
        $en = !$this->es();
        $rs = [];
        foreach ($this->mol->adj[$n] as $c => $o) {
            if (in_array($c, $excluir, true)) { continue; }
            $rs[] = $this->sustituyente($c, $n, 1);
        }
        if (!$rs) { return $this->info($base, 'simple'); }
        if ($base === ($en ? 'amino' : 'amino') && count($rs) === 1 && $rs[0]['nombre'] === ($en ? 'phenyl' : 'fenil')) {
            return $this->info($en ? 'anilino' : 'anilino', 'simple');
        }
        $grupos = [];
        foreach ($rs as $r) { $grupos[$r['nombre']]['info'] = $r; $grupos[$r['nombre']]['n'] = ($grupos[$r['nombre']]['n'] ?? 0) + 1; }
        uasort($grupos, fn($a, $b) => strcmp($a['info']['clave'], $b['info']['clave']) ?: strcmp($a['info']['nombre'], $b['info']['nombre']));
        $txt = ''; $i = 0;
        foreach ($grupos as $nm => $g) {
            $paren = $g['info']['tipo'] !== 'simple';
            $mult = $g['n'] > 1 ? ($g['info']['tipo'] === 'complejo' ? self::MULT_COMPLEJO[$g['n']] : self::MULT_SIMPLE[$g['n']]) : '';
            if ($mult !== '' && !$paren && preg_match('/^(sec|terc|tert)-/', $nm)) { $mult .= '-'; }
            $t = $mult . ($paren ? self::encerrar($nm) : $nm);
            $txt .= ($i > 0 && !$paren && $mult === '') ? '(' . $t . ')' : $t;
            $i++;
        }
        return $this->info($txt . $base, 'complejo');
    }

    private function amino(int $n, int $desde): array
    {
        return $this->nSustituido($n, [$desde], 'amino');
    }

    /** Acilamino (acetamido, benzamido…): N de amida unido a la estructura por el lado N. */
    private function acilamino(int $n, int $desde): array
    {
        $f = null;
        foreach ($this->mol->adj[$n] as $b => $o) { if (isset($this->fc[$b]) && $this->fc[$b]['t'] === 'amida' && $this->fc[$b]['n'] === $n) { $f = $b; } }
        $ac = $this->acilo($f, $n);
        $otros = array_filter(array_keys($this->mol->adj[$n]), fn($b) => $b !== $f && $b !== $desde);
        if ($otros) {                                    // [acetil(metil)amino]
            $rs = [$ac];
            foreach ($otros as $b) { $rs[] = $this->sustituyente($b, $n, 1); }
            usort($rs, fn($a, $b) => strcmp($a['clave'], $b['clave']));
            $txt = '';
            foreach ($rs as $i => $r) {
                $t = $r['tipo'] === 'simple' ? $r['nombre'] : self::encerrar($r['nombre']);
                $txt .= ($i > 0 && $r['tipo'] === 'simple') ? '(' . $t . ')' : $t;
            }
            return $this->info($txt . 'amino', 'complejo');
        }
        $nm = $ac['nombre'];
        if ($this->es()) {
            $nm = preg_replace(['/acetil$/', '/formil$/', '/carbonil$/', '/oil$/'], ['acetamido', 'formamido', 'carboxamido', 'amido'], $nm);
        } else {
            $nm = preg_replace(['/acetyl$/', '/formyl$/', '/carbonyl$/', '/oyl$/'], ['acetamido', 'formamido', 'carboxamido', 'amido'], $nm);
        }
        return $this->info($nm, $ac['tipo']);
    }

    /** Aciloxi (acetiloxi, benzoiloxi…): O de éster unido a la estructura por el lado alcohol. */
    private function aciloxi(int $f, int $o): array
    {
        $ac = $this->acilo($f, $o);
        $nm = ($ac['tipo'] === 'simple' ? $ac['nombre'] : self::encerrar($ac['nombre'])) . ($this->es() ? 'oxi' : 'oxy');
        return $this->info($nm, 'complejo');
    }

    /** Alcoxi / ariloxi (P-63.2.2): metoxi, etoxi, propoxi, butoxi, fenoxi; los demás «R-iloxi». */
    private function alcoxi(int $y, int $o): array
    {
        $r = $this->sustituyente($y, $o, 1);
        $nm = $r['nombre'];
        $en = !$this->es();
        $re = $en ? '/(meth|eth|prop|but|phen)yl$/' : '/(met|et|prop|but|fen)il$/';
        if (preg_match($re, $nm) && !preg_match('/an-\d+-(yl|il)$/', $nm)) {
            $nuevo = preg_replace($en ? '/yl$/' : '/il$/', $en ? 'oxy' : 'oxi', $nm);
            $nuevo = $en ? preg_replace('/phenoxy$/', 'phenoxy', $nuevo) : $nuevo;
            return $this->info($nuevo, $r['tipo'] === 'simple' ? 'simple' : 'complejo');
        }
        if ($r['tipo'] === 'simple') { return $this->info($nm . ($en ? 'oxy' : 'oxi'), 'simple'); }
        return $this->info(self::encerrar($nm) . ($en ? 'oxy' : 'oxi'), 'complejo');
    }

    /** Añade prefijos armados y conteos a un descriptor evaluado. */
    private function completar(array $d, bool $omitir): array
    {
        $a = $this->armar($d['subs'], $omitir);
        return array_merge($d, $a);
    }

    /* ======================= estructura principal y nombre ======================= */

    /** Todos los candidatos a estructura principal (cadenas y anillos, en todas las orientaciones). */
    private function candidatos(): array
    {
        $m = $this->mol;
        $cands = [];
        // sistemas de anillos
        foreach (array_keys($this->sistemas) as $k) {
            foreach ($this->candidatosSistema($k, ['padre' => true]) as $d) { $cands[] = $d; }
        }
        // cadenas
        $esq = [];
        foreach ($m->el as $a => $e) {
            if ($e === 'C' && $this->esEsqueleto($a, true) && !isset($this->anilloDe[$a])) { $esq[$a] = true; }
        }
        $fcP = [];
        if ($this->principal !== null) {
            foreach ($m->el as $a => $e) {
                if ($e === 'C' && $this->tipoFC($a, true) === $this->principal) { $fcP[$a] = $this->anclaFC($a); }
            }
        }
        $ext = [];
        foreach ($esq as $a => $x) {
            $g = 0;
            foreach ($m->adj[$a] as $b => $o) { if (isset($esq[$b])) { $g++; } }
            if ($g <= 1) { $ext[$a] = true; }
        }
        foreach ($fcP as $f => $anc) { if ($anc >= 0 && isset($esq[$anc])) { $ext[$anc] = true; } }
        $ext = array_keys($ext);
        $fcDe = [];
        foreach ($fcP as $f => $anc) { if ($anc >= 0) { $fcDe[$anc][] = $f; } }
        $secuencias = [];
        foreach ($ext as $u) {
            foreach ($ext as $v) {
                if ($v < $u) { continue; }
                $ruta = $u === $v ? [$u] : $this->camino($u, $v, $esq);
                if (!$ruta) { continue; }
                $opU = array_merge([null], $fcDe[$u] ?? []);
                $opV = array_merge([null], $fcDe[$v] ?? []);
                foreach ($opU as $fu) {
                    foreach ($opV as $fv) {
                        if ($fu !== null && $fu === $fv) { continue; }
                        $sec = $ruta;
                        if ($fu !== null) { array_unshift($sec, $fu); }
                        if ($fv !== null) { $sec[] = $fv; }
                        $secuencias[] = $sec;
                    }
                }
            }
        }
        foreach ($fcP as $f => $anc) {
            if ($anc < 0 || !isset($esq[$anc])) {
                if ($anc >= 0 && isset($fcP[$anc])) { if ($f < $anc) { $secuencias[] = [$f, $anc]; } }
                else { $secuencias[] = [$f]; }
            }
        }
        foreach ($secuencias as $sec) {
            foreach ([$sec, array_reverse($sec)] as $k => $s) {
                if ($k === 1 && count($s) === 1) { continue; }
                $d = $this->evaluar($s, ['padre' => true]);
                if ($d['invalido']) { continue; }
                $d['arom'] = false;
                $cands[] = $d;
            }
        }
        if (!$cands) { throw new NomenclaturaException('No se encontró una estructura principal válida (¿triple enlace hacia un sustituyente?).'); }
        // pre-filtro barato por los criterios de selección de la estructura
        $seqSel = array_slice($this->secuencia(), 0, 7);
        usort($cands, fn($a, $b) => $this->comparar($a, $b, $seqSel)[0]);
        $top = $cands[0];
        $finalistas = []; $rivales = [];
        foreach ($cands as $c) {
            if ($this->comparar($top, $c, $seqSel)[0] === 0) { $finalistas[] = $c; }
            else { $rivales[] = $c; }
        }
        foreach ($finalistas as &$c) { $c = $this->completar($c, false); $c['nombreCmp'] = $this->claveEster($c) . $c['prefijos']; }
        unset($c);
        usort($finalistas, fn($a, $b) => $this->comparar($a, $b)[0]);
        return [$finalistas, $rivales];
    }

    /**
     * Nombre completo. Devuelve: nombre, cadena (átomos numerados), anillo (bool), explicacion (array),
     * subs (prefijos con su localizador), radicofuncional (?string).
     */
    public function nombrar(): array
    {
        $m = $this->mol;
        if ($m->numCarbonos() === 0) { throw new NomenclaturaException('La estructura no contiene carbono.'); }
        if ($m->numCarbonos() > 60) { throw new NomenclaturaException('La estructura tiene más de 60 carbonos.'); }
        [$fin, $rivales] = $this->candidatos();
        $best = $fin[0];
        $r = $best['anillo'] ? $this->especial($best) : null;       // ensamblaje o nomenclatura multiplicativa
        if ($r === null) { $r = $this->construirNombre($best); }
        $r['explicacion'] = $this->explicar($best, $fin, $rivales, $r);
        if (!isset($r['cadena'])) { $r['cadena'] = $best['pos']; }
        if (!isset($r['mapa'])) {
            $r['mapa'] = [];
            foreach ($best['pos'] as $i => $a) { $r['mapa'][$a] = isset($best['lab']) ? $best['lab'][$a] : $i + 1; }
        }
        $r['anillo'] = $best['anillo'];
        $r['subs'] = $r['subsNombre'] ?? $best['subs'];
        $r['radicofuncional'] = $this->nombreClaseFuncional();
        return $r;
    }

    /* ---------- ensamblajes de anillos (P-28) y nomenclatura multiplicativa (P-15.3) ---------- */

    private function especial(array $best): ?array
    {
        $k = $best['sis'];
        $f = $this->firmaSistema($k);
        if ($f === null) { return null; }
        $iguales = [];
        foreach (array_keys($this->sistemas) as $j) {
            if ($j !== $k && $this->firmaSistema($j) === $f) { $iguales[] = $j; }
        }
        if (count($iguales) !== 1) { return null; }
        $j = $iguales[0];
        $uS = $uT = null;
        foreach ($this->sistemas[$k]['atomos'] as $a) {
            foreach ($this->mol->adj[$a] as $b => $o) {
                if (($this->anilloDe[$b] ?? -1) === $j) {
                    if ($o !== 1 || $uS !== null) { return null; }
                    $uS = $a; $uT = $b;
                }
            }
        }
        if ($uS !== null) { return $this->ensamblaje($k, $j, $uS, $uT); }
        return $this->multiplicativa($k, $j);
    }

    private static function prima($l): string { return $l . "'"; }

    /** Ensamblaje de dos sistemas idénticos unidos por un enlace sencillo: 1,1'-bifenilo, 1,1'-binaftaleno… */
    private function ensamblaje(int $k, int $j, int $uS, int $uT): ?array
    {
        $DS = $this->datosSistema($k);
        if (!empty($DS['ih'])) { return null; }
        $en = !$this->es();
        $cands = [];
        foreach ([[$k, $uS, $j, $uT], [$j, $uT, $k, $uS]] as [$A, $ua, $B, $ub]) {
            $NA = $this->datosSistema($A)['nums']; $NB = $this->datosSistema($B)['nums'];
            foreach ($NA as $na) {
                foreach ($NB as $nb) {
                    $off = count($na['pos']);
                    $pos = array_merge($na['pos'], $nb['pos']);
                    $locs = array_merge($na['locs'], array_map([self::class, 'prima'], $nb['locs']));
                    $ar = $na['aristas'];
                    foreach ($nb['aristas'] as [$x, $y]) { $ar[] = [$x + $off, $y + $off]; }
                    $d = $this->evaluar($pos, ['padre' => true, 'anillo' => true, 'locs' => $locs, 'aristas' => $ar]);
                    if ($d['invalido']) { continue; }
                    $lab = array_combine($pos, $locs);
                    $d['mlocs'] = $d['dlocs'] = $d['tlocs'] = []; $d['mn'] = $d['dn'] = 0;
                    $d['union'] = self::ordenarLocs([$lab[$ua], $lab[$ub]]);
                    $d['lab'] = $lab; $d['sisA'] = $A;
                    $d = $this->completar($d, false);
                    $d['nombreCmp'] = $d['prefijos'];
                    $cands[] = $d;
                }
            }
        }
        if (!$cands) { return null; }
        $seq = ['union', 'plocs', 'mlocs', 'dlocs', 'cuenta', 'locs', 'citados', 'nombre'];
        usort($cands, fn($a, $b) => $this->comparar($a, $b, $seq)[0]);
        $d = $cands[0];
        $princ = $d['princ'];
        $nk = count($princ);
        $t = $nk ? $princ[0]['t'] : null;
        if ($t === 'ester' && $nk > 1) { return null; }
        $clase = $DS['clase'];
        if ($clase === 'mono' && !empty($DS['arom'])) { $stem = $en ? 'biphenyl' : 'bifenil'; $fin = $en ? '' : 'o'; }
        elseif ($clase === 'mono') { $c = ($en ? 'cyclo' : 'ciclo') . $this->raiz($DS['n']) . ($en ? 'ane' : 'ano'); $stem = 'bi(' . $c . ')'; $fin = ''; }
        elseif ($clase === 'fus') { $stem = 'bi' . self::FUSIONADOS[$DS['plantilla']][$en ? 0 : 1] . ($en ? 'e' : 'o'); $fin = ''; }
        else { return null; }
        $ul = implode(',', $d['union']);
        $arm = $this->armar($d['subs'], false);
        if ($nk > 0) {
            $base = '[' . $ul . '-' . $stem . ']' . self::locsSufijo($d['plocs'], []) . ($nk > 1 ? self::MULT_SIMPLE[$nk] : '') . self::SUF[$this->lang][$t][1];
            $nombre = $arm['prefijos'] . $base;
        } else {
            $nombre = self::unir([$arm['prefijos'], $ul . '-' . $stem . $fin]);
        }
        if (!$en && $t === 'acido') { $nombre = 'ácido ' . $nombre; }
        if ($t === 'ester') { $nombre = $this->componerEster($princ, $nombre); }
        $mapa = [];
        foreach ($d['lab'] as $a => $l) { if (is_int($l)) { $mapa[$a] = $l; } }
        return ['nombre' => $nombre, 'notas' => ['ensamblaje'], 'nsubs' => [], 'omitir' => false, 'grupos' => $arm['grupos'],
            'cadena' => $d['pos'], 'mapa' => $mapa, 'subsNombre' => $d['subs']];
    }

    /** Nombre multiplicativo: 1,1'-metilendibenceno, 4,4'-oxidianilina, 4,4'-(propano-2,2-diil)difenol… */
    private function multiplicativa(int $k, int $j): ?array
    {
        $m = $this->mol;
        $en = !$this->es();
        $enS = array_fill_keys($this->sistemas[$k]['atomos'], true);
        $enT = array_fill_keys($this->sistemas[$j]['atomos'], true);
        // componentes del resto de la molécula
        $comp = []; $idc = 0;
        foreach (array_keys($m->el) as $a) {
            if (isset($enS[$a]) || isset($enT[$a]) || isset($comp[$a])) { continue; }
            $pila = [$a]; $comp[$a] = $idc;
            while ($pila) {
                $x = array_pop($pila);
                foreach ($m->adj[$x] as $y => $o) {
                    if (!isset($enS[$y]) && !isset($enT[$y]) && !isset($comp[$y])) { $comp[$y] = $idc; $pila[] = $y; }
                }
            }
            $idc++;
        }
        $unionS = []; $unionT = [];
        foreach ($comp as $a => $c) {
            foreach ($m->adj[$a] as $b => $o) {
                if (isset($enS[$b])) { $unionS[$c][] = [$b, $a, $o]; }
                if (isset($enT[$b])) { $unionT[$c][] = [$b, $a, $o]; }
            }
        }
        $puente = array_values(array_intersect(array_keys($unionS), array_keys($unionT)));
        if (count($puente) !== 1) { return null; }
        $c = $puente[0];
        if (count($unionS[$c]) !== 1 || count($unionT[$c]) !== 1) { return null; }
        [$sA, $c1, $o1] = $unionS[$c][0]; [$tA, $c2, $o2] = $unionT[$c][0];
        if ($o1 !== 1 || $o2 !== 1) { return null; }
        $C = array_keys(array_filter($comp, fn($x) => $x === $c));
        foreach ($C as $a) { if (isset($this->anilloDe[$a])) { return null; } }
        // el conector no puede contener grupos del grupo principal
        if ($this->principal !== null) {
            foreach ($C as $a) {
                if ($this->tipoFC($a, false) === $this->principal || ($this->principal === 'aldehido' && isset($this->ald[$a]))) { return null; }
                if ($m->el[$a] === 'N' && $this->principal === 'amina' && !$m->esNitroN($a) && !$this->esNAmida($a)) { return null; }
                if ($m->el[$a] === 'O' && count($m->adj[$a]) === 1) {
                    $cc = $m->vecinos($a)[0];
                    if ($m->orden($a, $cc) === 1 && $this->principal === 'alcohol' && !isset($this->fc[$cc])) { return null; }
                    if ($m->orden($a, $cc) === 2 && $this->principal === 'cetona' && !isset($this->fc[$cc]) && !isset($this->ald[$cc])) { return null; }
                }
            }
        }
        $conn = $this->conector($C, $c1, $c2, $sA, $tA);
        if ($conn === null) { return null; }
        // unidades multiplicadas
        $u = [];
        foreach ([[$k, $sA, $c1], [$j, $tA, $c2]] as [$X, $xa, $cx]) {
            if (!in_array($this->principal, [null, 'alcohol', 'amina', 'acido', 'aldehido', 'nitrilo', 'amida', 'cetona'], true)) { return null; }
            $cands = $this->candidatosSistema($X, ['padre' => true, 'excluir' => [$cx => true], 'union' => $xa]);
            foreach ($cands as &$cd) { $cd = $this->completar($cd, false); $cd['nombreCmp'] = $cd['prefijos']; }
            unset($cd);
            if (!$cands) { return null; }
            usort($cands, fn($a, $b) => $this->comparar($a, $b)[0]);
            $bx = $cands[0];
            $nm = $this->construirNombre($bx, true);
            if ($nm['nsubs']) { return null; }
            $u[] = ['d' => $bx, 'nombre' => $nm['nombre'], 'loc' => $bx['union'][0], 'grupos' => $nm['grupos']];
        }
        if ($u[0]['nombre'] !== $u[1]['nombre'] || $u[0]['loc'] !== $u[1]['loc']) { return null; }
        $un = $u[0]['nombre'];
        $acido = false;
        if (!$en && strpos($un, 'ácido ') === 0) { $un = substr($un, strlen('ácido ')); $acido = true; }
        $cuerpo = $u[0]['d']['cuenta'] > 0 ? 'bis' . self::encerrar($un) : 'di' . $un;
        $ct = $conn['tipo'] === 'simple' ? $conn['nombre'] : self::encerrar($conn['nombre']);
        $nombre = $u[0]['loc'] . ',' . self::prima($u[0]['loc']) . '-' . $ct . $cuerpo;
        if ($acido) { $nombre = 'ácido ' . $nombre; }
        $mapa = [];
        foreach ($u[0]['d']['lab'] as $a => $l) { if (is_int($l)) { $mapa[$a] = $l; } }
        return ['nombre' => $nombre, 'notas' => ['multiplicativa'], 'nsubs' => [], 'omitir' => false, 'grupos' => $u[0]['grupos'],
            'cadena' => $u[0]['d']['pos'], 'mapa' => $mapa, 'subsNombre' => $u[0]['d']['subs'], 'conector' => $conn['nombre']];
    }

    /** Sustituyente divalente que une las unidades: oxi, azanodiil, metilen, carbonil, etano-1,2-diil, propano-2,2-diil… */
    private function conector(array $C, int $c1, int $c2, int $sA, int $tA): ?array
    {
        $m = $this->mol;
        $en = !$this->es();
        if (count($C) === 1) {
            $x = $C[0];
            if ($m->el[$x] === 'O') { return $this->info($en ? 'oxy' : 'oxi', 'simple'); }
            if ($m->el[$x] === 'N' && count($m->adj[$x]) === 2) { return $this->info($en ? 'azanediyl' : 'azanodiil', 'simple'); }
        }
        if ($m->el[$c1] !== 'C' || $m->el[$c2] !== 'C' || !$this->esEsqueleto($c1, false) || !$this->esEsqueleto($c2, false)) { return null; }
        $E = [];
        foreach ($C as $a) { if ($m->el[$a] === 'C' && $this->esEsqueleto($a, false)) { $E[$a] = true; } }
        foreach ($C as $a) {                        // heteroátomos del conector: solo como prefijos terminales
            if ($m->el[$a] !== 'C' && count($m->adj[$a]) > 1 && !$m->esNitroN($a)) {
                $puenteo = 0;
                foreach ($m->adj[$a] as $b => $o) { if (isset($E[$b])) { $puenteo++; } }
                if ($puenteo > 1) { return null; }
            }
        }
        $base = $this->camino($c1, $c2, $E);
        if (!$base) { return null; }
        $ext = [$c1 => true, $c2 => true];
        foreach (array_keys($E) as $a) {
            $g = 0;
            foreach ($m->adj[$a] as $b => $o) { if (isset($E[$b])) { $g++; } }
            if ($g <= 1) { $ext[$a] = true; }
        }
        $cands = [];
        foreach (array_keys($ext) as $u) {
            foreach (array_keys($ext) as $v) {
                $ruta = $u === $v ? [$u] : $this->camino($u, $v, $E);
                if (!$ruta || !in_array($c1, $ruta, true) || !in_array($c2, $ruta, true)) { continue; }
                $d = $this->evaluar($ruta, ['excluir' => [$sA => true, $tA => true]]);
                if ($d['invalido']) { continue; }
                $d['plocs'] = self::ordenarLocs([array_search($c1, $ruta, true) + 1, array_search($c2, $ruta, true) + 1]);
                $d['princ'] = [];
                $d = $this->completar($d, count($ruta) === 1);
                $d['nombreCmp'] = $d['prefijos'];
                $cands[] = $d;
            }
        }
        if (!$cands) { return null; }
        usort($cands, fn($a, $b) => $this->comparar($a, $b, ['len', 'mn', 'dn', 'plocs', 'mlocs', 'dlocs', 'cuenta', 'locs', 'citados', 'nombre'])[0]);
        $b = $cands[0];
        if ($b['len'] === 1) {
            if (count($b['subs']) === 1 && $b['subs'][0]['info']['nombre'] === 'oxo') { return $this->info($en ? 'carbonyl' : 'carbonil', 'simple'); }
            $nm = $b['prefijos'] . ($en ? 'methylene' : 'metilen');
            return $this->info($nm, $b['cuenta'] > 0 ? 'complejo' : 'simple');
        }
        $nm = $b['prefijos'] . $this->ensamblar($this->raiz($b['len']), $b['dlocs'], $b['tlocs'], true, '-' . implode(',', $b['plocs']) . '-' . ($en ? 'diyl' : 'diil'));
        return $this->info($nm, $b['cuenta'] > 0 ? 'complejo' : 'compuesto');
    }

    private const SUF = [
        'en' => [
            'acido' => ['oic acid', 'carboxylic acid'], 'ester' => ['oate', 'carboxylate'], 'amida' => ['amide', 'carboxamide'],
            'nitrilo' => ['nitrile', 'carbonitrile'], 'aldehido' => ['al', 'carbaldehyde'], 'cetona' => ['one', 'one'],
            'alcohol' => ['ol', 'ol'], 'amina' => ['amine', 'amine'],
        ],
        'es' => [
            'acido' => ['oico', 'carboxílico'], 'ester' => ['oato', 'carboxilato'], 'amida' => ['amida', 'carboxamida'],
            'nitrilo' => ['nitrilo', 'carbonitrilo'], 'aldehido' => ['al', 'carbaldehído'], 'cetona' => ['ona', 'ona'],
            'alcohol' => ['ol', 'ol'], 'amina' => ['amina', 'amina'],
        ],
    ];
    private const FC_CADENA = ['acido', 'ester', 'amida', 'nitrilo', 'aldehido'];

    private function construirNombre(array $best, bool $unidad = false): array
    {
        $m = $this->mol;
        $en = !$this->es();
        $pos = $best['pos'];
        $N = count($pos);
        $anillo = $best['anillo'];
        $arom = !empty($best['arom']);
        $princ = $best['princ'];
        $k = count($princ);
        $t = $k ? $princ[0]['t'] : null;
        $notas = [];

        $subs = array_values(array_filter($best['subs'], fn($x) => !self::esN($x['loc'])));
        $nsubs = array_values(array_filter($best['subs'], fn($x) => self::esN($x['loc'])));
        $nPref = count($subs) + ($unidad ? 1 : 0);          // en una unidad multiplicativa, el punto de unión cuenta
        $chainFC = !$anillo && in_array($t, self::FC_CADENA, true);

        // ---- sistemas policíclicos: fusionados, von Baeyer, espiro ----
        if ($anillo && ($best['clase'] ?? 'mono') !== 'mono') {
            $arm = $this->armar(array_merge($subs, $nsubs), false);
            $sufijo = '';
            if ($k > 0) {
                $sufijo = self::locsSufijo($best['plocs'], $best['alocs'] ?? []) . ($k > 1 ? self::MULT_SIMPLE[$k] : '') . self::SUF[$this->lang][$t][1];
            }
            [$hid, $base] = $this->nombreSistema($best, $sufijo);
            $nombre = self::unir([$arm['prefijos'], $hid, $base]);
            if (!$en && $t === 'acido') { $nombre = 'ácido ' . $nombre; }
            if ($t === 'ester') { $nombre = $this->componerEster($princ, $nombre); }
            $notas[] = 'sistema_' . $best['clase'];
            if (!empty($best['ilocs'])) { $notas[] = 'h_indicado'; }
            if (!empty($best['alocs'])) { $notas[] = 'h_anadido'; }
            if (!empty($best['hlocs'])) { $notas[] = 'hidro'; }
            return ['nombre' => $nombre, 'notas' => $notas, 'nsubs' => $nsubs, 'omitir' => false, 'grupos' => $arm['grupos']];
        }

        // ---- omisión de localizadores (P-14.3.4) ----
        $omitPref = false; $omitSuf = $chainFC; $omitInsat = false;
        if (!$anillo && $N === 1) { $omitPref = $omitSuf = true; }
        if (!$anillo && $N === 2 && !$chainFC && $nPref + $k === 1) { $omitPref = $omitSuf = $omitInsat = true; }
        if ($anillo && $nPref + $k === 1 && !$best['dlocs'] && !$best['tlocs']) { $omitPref = $omitSuf = true; }
        if ($nPref === 0 && $k === 0 && $best['mn'] === 1 && (($anillo) || (!$anillo && $N <= 3))) { $omitInsat = true; }
        if ($nPref === 0 && $k === 0 && !$anillo && $N === 2) { $omitInsat = true; }
        if ($k === 0 && $nPref > 1) {                     // sustitución total por un mismo prefijo
            $nombres = array_unique(array_map(fn($s) => $s['info']['nombre'], $subs));
            $sinH = true;
            foreach ($pos as $a) { if ($m->hidrogenos($a) > 0) { $sinH = false; } }
            $mono = true;
            foreach ($subs as $s) { if (($s['orden'] ?? 1) === 2 || $s['info']['nombre'] === 'oxo') { $mono = false; } }
            if (count($nombres) === 1 && $sinH && $mono) { $omitPref = true; $omitInsat = true; }
        }

        // ---- nombres retenidos (PIN) ----
        $retenido = null;
        if ($anillo && $arom && $k === 1 && $t !== 'cetona') {
            $retenido = ['alcohol' => ['phenol', 'fenol'], 'amina' => ['aniline', 'anilina'], 'acido' => ['benzoic acid', 'benzoico'],
                'aldehido' => ['benzaldehyde', 'benzaldehído'], 'amida' => ['benzamide', 'benzamida'],
                'nitrilo' => ['benzonitrile', 'benzonitrilo'], 'ester' => ['benzoate', 'benzoato']][$t];
            $notas[] = 'retenido_benceno';
        } elseif ($chainFC && $k === 1 && $N === 2) {
            $retenido = ['acido' => ['acetic acid', 'acético'], 'ester' => ['acetate', 'acetato'], 'amida' => ['acetamide', 'acetamida'],
                'nitrilo' => ['acetonitrile', 'acetonitrilo'], 'aldehido' => ['acetaldehyde', 'acetaldehído']][$t];
            if ($t !== 'amida') { $omitPref = true; }
            $notas[] = 'retenido_acetico';
        } elseif ($chainFC && $k === 1 && $N === 1) {
            $retenido = ['acido' => ['formic acid', 'fórmico'], 'ester' => ['formate', 'formiato'], 'amida' => ['formamide', 'formamida'],
                'nitrilo' => ['formonitrile', 'formonitrilo'], 'aldehido' => ['formaldehyde', 'formaldehído']][$t];
            $notas[] = 'retenido_formico';
        } elseif ($chainFC && $k === 2 && $N === 2 && in_array($t, ['acido', 'ester'], true)) {
            $retenido = ['acido' => ['oxalic acid', 'oxálico'], 'ester' => ['oxalate', 'oxalato']][$t];
            $notas[] = 'retenido_oxalico';
        }

        $todos = array_merge($subs, $nsubs);
        $arm = $this->armar($todos, $omitPref);
        $pref = $arm['prefijos'];

        if ($retenido !== null) {
            $stem = $retenido[$en ? 0 : 1];
            if (!$en && $t === 'acido') { $nombre = 'ácido ' . $pref . $stem; }
            else { $nombre = $pref . $stem; }
        } else {
            $sufijo = '';
            if ($k > 0) {
                $s = self::SUF[$this->lang][$t][$anillo ? 1 : 0];
                $mult = $k > 1 ? self::MULT_SIMPLE[$k] : '';
                $sufijo = ($omitSuf ? '' : '-' . implode(',', $best['plocs']) . '-') . $mult . $s;
            }
            if ($anillo && $arom) {
                $stem = $this->elidir(($en ? 'benzen§' : 'bencen§') . $sufijo);
            } else {
                $raiz = $anillo ? (($en ? 'cyclo' : 'ciclo') . $this->raiz($N)) : $this->raiz($N);
                $stem = $this->ensamblar($raiz, $best['dlocs'], $best['tlocs'], !$omitInsat, $sufijo);
            }
            $nombre = $pref . $stem;
            if (!$en && $t === 'acido') { $nombre = 'ácido ' . $nombre; }
        }
        // éster: nombre del grupo alquilo (R') como palabra aparte
        if ($t === 'ester') {
            $nombre = $this->componerEster($princ, $nombre);
        }
        // hidrocarburos con nombre retenido (P-22.1.3; acetileno P-14.3.4.2; anisol P-63.2)
        $hid = $en
            ? ['ethyne' => 'acetylene', 'methylbenzene' => 'toluene', 'methoxybenzene' => 'anisole',
               '1,2-dimethylbenzene' => '1,2-xylene', '1,3-dimethylbenzene' => '1,3-xylene', '1,4-dimethylbenzene' => '1,4-xylene']
            : ['etino' => 'acetileno', 'metilbenceno' => 'tolueno', 'metoxibenceno' => 'anisol',
               '1,2-dimetilbenceno' => '1,2-xileno', '1,3-dimetilbenceno' => '1,3-xileno', '1,4-dimetilbenceno' => '1,4-xileno'];
        if (isset($hid[$nombre])) { $notas[] = 'retenido_hidrocarburo:' . $nombre; $nombre = $hid[$nombre]; }
        return ['nombre' => $nombre, 'notas' => $notas, 'nsubs' => $nsubs, 'omitir' => $omitPref, 'grupos' => $arm['grupos']];
    }

    /** Nombre de éster: «ethyl propanoate» / «propanoato de etilo»; con grupos distintos, «1-ethyl 5-methyl …» (P-65.6.3.2). */
    private function componerEster(array $princ, string $acido): string
    {
        $en = !$this->es();
        $rs = [];
        foreach ($princ as $p) {
            $d = $this->fc[$p['fc']];
            $info = $this->sustituyente($d['r'], $d['o'], 1);
            $rs[$info['nombre']]['info'] = $info;
            $rs[$info['nombre']]['locs'][] = $p['loc'];
        }
        uasort($rs, fn($a, $b) => strcmp($a['info']['clave'], $b['info']['clave']) ?: strcmp($a['info']['nombre'], $b['info']['nombre']));
        $conLoc = count($rs) > 1;
        $partes = [];
        foreach ($rs as $nm => $g) {
            $base = $en ? $nm : $nm . 'o';
            $paren = $g['info']['tipo'] !== 'simple';
            $n = count($g['locs']);
            $mult = '';
            if ($n > 1) {
                $mult = $g['info']['tipo'] === 'complejo' ? self::MULT_COMPLEJO[$n] : self::MULT_SIMPLE[$n];
                if (!$paren && preg_match('/^(sec|terc|tert)-/', $nm)) { $mult .= '-'; }
            }
            sort($g['locs']);
            $txt = $mult . (($paren && ($n > 1 || $conLoc)) ? self::encerrar($base) : $base);
            $partes[] = ($conLoc ? implode(',', $g['locs']) . '-' : '') . $txt;
        }
        if ($en) { return implode(' ', $partes) . ' ' . $acido; }
        $ult = array_pop($partes);
        $lista = $partes ? implode(', ', $partes) . ' y ' . $ult : $ult;
        return $acido . ' de ' . $lista;
    }

    /** Clave de desempate: localizador más bajo para el grupo R' de éster citado primero. */
    private function claveEster(array $d): string
    {
        $r = [];
        foreach ($d['princ'] as $p) {
            if ($p['t'] !== 'ester') { continue; }
            $f = $this->fc[$p['fc']];
            $r[] = [$this->sustituyente($f['r'], $f['o'], 1)['clave'], $p['loc']];
        }
        sort($r);
        return implode(',', array_map(fn($x) => sprintf('%03d', $x[1]), $r)) . '|';
    }

    /* ======================= nombres de clase funcional ======================= */

    /**
     * Nombre de clase funcional (radicofuncional), aceptado en nomenclatura general (no PIN):
     * haluros de alquilo, alcoholes, éteres, cetonas y alquilaminas con grupos R sencillos.
     */
    public function nombreClaseFuncional(): ?string
    {
        $m = $this->mol;
        $en = !$this->es();
        $cl = $this->clases;
        $otras = fn(array $perm) => array_diff(array_keys($cl), array_merge($perm, ['alqueno', 'alquino']));
        $simple = fn(array $i) => $i['tipo'] === 'simple';
        // R–X
        if (($cl['halogeno'] ?? 0) === 1 && !$otras(['halogeno'])) {
            foreach ($m->el as $x => $e) {
                if (!in_array($e, MolOrg::HALOGENOS, true)) { continue; }
                $c = $m->vecinos($x)[0];
                if (!empty($m->arom[$c])) { return null; }
                $R = $this->sustituyente($c, $x, 1);
                if (!$simple($R)) { return null; }
                $h = self::HALURO[$this->lang][$e];
                return $en ? $R['nombre'] . ' ' . $h : $h . ' de ' . $R['nombre'] . 'o';
            }
        }
        // R–OH
        if (($cl['alcohol'] ?? 0) === 1 && !$otras(['alcohol'])) {
            foreach ($m->el as $o => $e) {
                if ($e !== 'O' || count($m->adj[$o]) !== 1) { continue; }
                $c = $m->vecinos($o)[0];
                if (!empty($m->arom[$c])) { return null; }
                $R = $this->sustituyente($c, $o, 1);
                if (!$simple($R) || !preg_match('/il$|yl$/', $R['nombre'])) { return null; }
                return $en ? $R['nombre'] . ' alcohol' : 'alcohol ' . preg_replace('/il$/', 'ílico', $R['nombre']);
            }
        }
        // R–O–R'
        if (($cl['eter'] ?? 0) === 1 && !$otras(['eter'])) {
            foreach ($m->el as $o => $e) {
                if ($e !== 'O' || count($m->adj[$o]) !== 2) { continue; }
                [$c1, $c2] = $m->vecinos($o);
                $R = [$this->sustituyente($c1, $o, 1), $this->sustituyente($c2, $o, 1)];
                return $this->listaR($R, $en ? 'ether' : 'éter');
            }
        }
        // R–CO–R'
        if (($cl['cetona'] ?? 0) === 1 && !$otras(['cetona'])) {
            foreach ($m->el as $c => $e) {
                if ($e !== 'C' || isset($this->anilloDe[$c])) { continue; }
                $esK = false; $R = [];
                foreach ($m->adj[$c] as $b => $o) {
                    if ($m->el[$b] === 'O' && $o === 2) { $esK = true; }
                    elseif ($m->el[$b] === 'C') { $R[] = $b; }
                }
                if (!$esK || count($R) !== 2 || isset($this->ald[$c])) { continue; }
                $inf = [$this->sustituyente($R[0], $c, 1), $this->sustituyente($R[1], $c, 1)];
                return $this->listaR($inf, $en ? 'ketone' : 'cetona');
            }
        }
        // aminas R-NH2, R2NH, R3N
        if (($cl['amina'] ?? 0) === 1 && !$otras(['amina'])) {
            foreach ($m->el as $n => $e) {
                if ($e !== 'N') { continue; }
                $inf = [];
                foreach ($m->adj[$n] as $c => $o) {
                    if (!empty($m->arom[$c])) { return null; }
                    $inf[] = $this->sustituyente($c, $n, 1);
                }
                foreach ($inf as $i) { if (!$simple($i) || $i['tipo'] === 'compuesto') { return null; } }
                usort($inf, fn($a, $b) => strcmp($a['clave'], $b['clave']));
                $cnt = [];
                foreach ($inf as $i) { $cnt[$i['nombre']] = ($cnt[$i['nombre']] ?? 0) + 1; }
                $txt = ''; $i = 0;
                foreach ($cnt as $nm => $c) {
                    $t = ($c > 1 ? self::MULT_SIMPLE[$c] : '') . $nm;
                    $txt .= ($i++ > 0 && $c === 1) ? '(' . $t . ')' : $t;       // etil(metil)amina (P-62.2)
                }
                return $txt . ($en ? 'amine' : 'amina');
            }
        }
        return null;
    }

    private function listaR(array $R, string $clase): ?string
    {
        foreach ($R as $r) { if ($r['tipo'] === 'complejo') { return null; } }
        usort($R, fn($a, $b) => strcmp($a['clave'], $b['clave']));
        $nm = array_map(fn($r) => $r['tipo'] === 'simple' ? $r['nombre'] : self::encerrar($r['nombre']), $R);
        if ($R[0]['nombre'] === $R[1]['nombre']) {
            return 'di' . ($R[0]['tipo'] === 'simple' ? '' : '') . $nm[0] . ' ' . $clase;
        }
        return $nm[0] . ' ' . $nm[1] . ' ' . $clase;
    }

    /* ======================= explicación didáctica ======================= */

    private static function fmt(array $l): string { return '{' . implode(',', $l) . '}'; }

    private function textoCriterio(string $k, array $a, array $b, bool $numeracion): string
    {
        switch ($k) {
            case 'pn': return 'contiene el mayor número de grupos principales (' . count($a['princ']) . ' frente a ' . count($b['princ']) . ') (P-44.1.1)';
            case 'anillo': return 'a igualdad de grupos principales, el anillo es preferido a la cadena (P-44.1.2.2)';
            case 'len': return 'es la cadena más larga (' . $a['len'] . ' C frente a ' . $b['len'] . ' C) (P-44.3)';
            case 'mn': return 'tiene el mayor número de enlaces múltiples (' . $a['mn'] . ' frente a ' . $b['mn'] . ') (P-44.4.1.1)';
            case 'dn': return 'tiene el mayor número de dobles enlaces (P-44.4.1.1)';
            case 'plocs':
                return $a['princ']
                    ? 'da los localizadores más bajos al grupo principal (sufijo): ' . self::fmt($a['plocs']) . ' frente a ' . self::fmt($b['plocs']) . ' (P-31.1.4.2.1)'
                    : 'da el localizador más bajo a la valencia libre (P-31.1.4.2.1)';
            case 'mlocs': return 'da los localizadores más bajos a los enlaces múltiples (terminaciones -eno/-ino) en conjunto: ' . self::fmt($a['mlocs']) . ' frente a ' . self::fmt($b['mlocs']) . ' (P-31.1.4.2.4)';
            case 'dlocs': return 'a igualdad, da los localizadores más bajos a los dobles enlaces: ' . self::fmt($a['dlocs']) . ' frente a ' . self::fmt($b['dlocs']) . ' (P-31.1.4.2.4)';
            case 'cuenta': return 'tiene el mayor número de sustituyentes (' . $a['cuenta'] . ' frente a ' . $b['cuenta'] . ') (P-45.2.1)';
            case 'locs': return 'da los localizadores más bajos al conjunto de prefijos: ' . self::fmt($a['locs']) . ' frente a ' . self::fmt($b['locs']) . ($numeracion ? ' (P-31.1.4.3.4, primer punto de diferencia)' : ' (P-45.2.2)');
            case 'citados':
                $g = reset($a['grupos']);
                return 'el localizador más bajo corresponde al prefijo citado primero en orden alfanumérico (' . ($g['info']['nombre'] ?? '') . ') (P-31.1.4.3.4, P-45.2.3)';
            case 'nombre': return 'conduce al nombre que aparece primero en orden alfanumérico (P-45.5)';
            case 'nr': return 'tiene más anillos (' . ($a['nr'] ?? 0) . ' frente a ' . ($b['nr'] ?? 0) . ') (P-44.2)';
            case 'insat': return 'tiene más átomos insaturados, es decir, menor grado de hidrogenación (P-44.4.1)';
            case 'ilocs': return 'da el localizador más bajo al hidrógeno indicado: ' . self::fmt($a['ilocs'] ?? []) . ' frente a ' . self::fmt($b['ilocs'] ?? []) . ' (P-31.1.4.2.4)';
            case 'union': return 'da los localizadores más bajos a los puntos de unión: ' . self::fmt($a['union'] ?? []) . ' frente a ' . self::fmt($b['union'] ?? []) . ' (P-28.2, P-15.3)';
            case 'alocs': return 'da el localizador más bajo al hidrógeno añadido: ' . self::fmt($a['alocs'] ?? []) . ' frente a ' . self::fmt($b['alocs'] ?? []) . ' (P-31.1.4)';
            case 'hlocs': return 'da los localizadores más bajos a los prefijos hidro: ' . self::fmt($a['hlocs'] ?? []) . ' frente a ' . self::fmt($b['hlocs'] ?? []) . ' (P-31.1.4.2.4)';
        }
        return '';
    }

    private function explicar(array $best, array $fin, array $rivales, array $res): array
    {
        $exp = [];
        $cl = $this->clases;
        $p = $this->principal;
        $es = true;
        // 1. grupos característicos
        $pres = [];
        foreach (array_merge(self::CLASES, ['eter', 'halogeno', 'nitro']) as $c) {
            if (!empty($cl[$c])) { $pres[] = self::NOMBRE_CLASE[$c] . ($cl[$c] > 1 ? ' (' . $cl[$c] . ')' : ''); }
        }
        if ($p !== null) {
            [$sc, $sa, $pre, $reg] = self::INFO_CLASE[$p];
            $enAnillo = $best['anillo'] && in_array($p, self::FC_CADENA, true);
            $exp[] = 'Grupos característicos: ' . implode(', ', $pres) . '. El de mayor prioridad (orden de clases, P-41) es ' . self::NOMBRE_CLASE[$p]
                . ': es el grupo principal y se expresa como sufijo ' . ($enAnillo ? $sa . ', porque su carbono no forma parte del anillo' : $sc) . ' (' . $reg . ').';
            $otros = [];
            foreach ($cl as $c => $n) {
                if ($c === $p || !isset(self::INFO_CLASE[$c])) { continue; }
                $otros[] = self::NOMBRE_CLASE[$c] . ' → ' . self::INFO_CLASE[$c][2];
            }
            if (!empty($cl['eter'])) { $otros[] = 'éter → alcoxi (P-63.2.2)'; }
            if (!empty($cl['halogeno'])) { $otros[] = 'halógeno → fluoro, cloro, bromo, yodo (P-61.3)'; }
            if (!empty($cl['nitro'])) { $otros[] = 'nitro → nitro (P-61.5)'; }
            if ($otros) { $exp[] = 'Los demás grupos se citan como prefijos: ' . implode('; ', $otros) . '.'; }
        } elseif ($pres) {
            $exp[] = 'Grupos presentes: ' . implode(', ', $pres) . '. Ninguno se expresa como sufijo: los éteres (alcoxi), halógenos y el grupo nitro se nombran siempre como prefijos (P-61, P-63.2), sobre el hidruro progenitor.';
        }
        // ensamblajes y nombres multiplicativos: explicación propia
        if (array_intersect($res['notas'], ['ensamblaje', 'multiplicativa'])) {
            $D = $this->datosSistema($best['sis']);
            $comp = $D['clase'] === 'mono' ? (!empty($D['arom']) ? 'anillos de benceno' : 'anillos de ' . $D['n'] . ' carbonos') : 'sistemas «' . ($D['clase'] === 'fus' ? self::FUSIONADOS[$D['plantilla']][1] . 'o' : ($D['corchete'] ?? '')) . '»';
            if (in_array('ensamblaje', $res['notas'], true)) {
                $exp[] = 'Estructura principal: ensamblaje de dos ' . $comp . ' idénticos, unidos directamente. Un ensamblaje es preferido a cada uno de sus componentes (P-28, P-44.1.2).';
            } else {
                $exp[] = 'Estructura principal: dos unidades idénticas (' . $comp . ', con el mismo grupo principal y los mismos prefijos) unidas por el conector «' . ($res['conector'] ?? '') . '». Se usa nomenclatura multiplicativa en lugar de tratar una unidad como sustituyente de la otra (P-15.3, P-51.3).';
            }
            $todos = $res['subsNombre'] ?? [];
            if ($todos) {
                $t = [];
                foreach ($todos as $sb) { $t[] = $sb['info']['nombre'] . ' en ' . (self::esN($sb['loc']) ? '' : 'C') . $sb['loc']; }
                $exp[] = 'Prefijos ' . (in_array('multiplicativa', $res['notas'], true) ? 'de cada unidad' : 'del ensamblaje') . ': ' . implode('; ', $t) . '.';
            }
            foreach ($res['notas'] as $n) {
                if ($n === 'ensamblaje') { $exp[] = 'Dos sistemas de anillos idénticos unidos por un enlace sencillo forman un ensamblaje: 1,1\'-bifenilo, 1,1\'-binaftaleno, 1,1\'-bi(ciclohexano). El segundo componente se numera con primas y los puntos de unión reciben los localizadores más bajos; después se aplican las reglas habituales al grupo principal y a los prefijos (P-28).'; }
                if ($n === 'multiplicativa') { $exp[] = 'Nombre multiplicativo: localizadores de unión (con prima en la segunda unidad) + conector (metilen, oxi, carbonil, azanodiil, etano-1,2-diil, propano-2,2-diil…) + di (o bis si la unidad lleva prefijos) + nombre de la unidad. En cada unidad, el grupo principal recibe el localizador más bajo y después el punto de unión (P-15.3).'; }
            }
            if (count($res['grupos']) > 1) {
                $exp[] = 'Los prefijos se citan en orden alfanumérico: ' . implode(' < ', array_keys($res['grupos'])) . ' (P-14.5).';
            }
            return array_values(array_filter($exp));
        }
        // 2. estructura principal
        $pos = $best['pos'];
        $N = count($pos);
        if ($best['anillo']) {
            $D = $this->datosSistema($best['sis']);
            switch ($D['clase']) {
                case 'fus':
                    $tipo = 'sistema de dos anillos fusionados (' . $D['n'] . ' C) con el nombre retenido «' . self::FUSIONADOS[$D['plantilla']][1] . 'o»'; break;
                case 'vb':
                    $tipo = 'sistema bicíclico con puente ' . $D['corchete'] . ' (nomenclatura de von Baeyer)'; break;
                case 'espiro':
                    $tipo = 'sistema espiro ' . $D['corchete']; break;
                default:
                    $tipo = !empty($best['arom']) ? 'anillo de benceno' : 'anillo de ' . $N . ' carbonos';
            }
            $txt = 'Estructura principal: ' . $tipo;
        } else {
            $txt = 'Cadena principal de ' . $N . ' carbono' . ($N > 1 ? 's' : '');
        }
        $rival = null; $crit = '';
        $clave = function ($d) { $c = $d['pos']; sort($c); return implode(',', $c); };
        $kb = $clave($best);
        foreach ($fin as $d) { if ($clave($d) !== $kb) { [$s, $crit] = $this->comparar($best, $d); if ($s !== 0) { $rival = $d; break; } } }
        if ($rival === null) {
            foreach ($rivales as $d) {
                if ($clave($d) === $kb) { continue; }
                [$s, $c] = $this->comparar($best, $d);
                if ($rival === null || $this->nivel($c) > $this->nivel($crit)) { $rival = $d; $crit = $c; }
            }
        }
        if ($rival !== null && $crit !== '') {
            $txt .= ': ' . $this->textoCriterio($crit, $best, $rival, false) . '.';
        } else {
            $txt .= '.';
        }
        $exp[] = $txt;
        if ($best['mn'] > 0) {
            $exp[] = 'Los enlaces múltiples de la estructura principal se indican con las terminaciones «-eno» (doble) e «-ino» (triple) y sus localizadores; ' .
                'si hay varios se usan multiplicadores (dieno, triino) y se añade «a» a la raíz (buta-1,3-dieno) (P-31.1.4.2.4).';
        }
        // 3. prefijos
        $todos = $best['subs'];
        if ($todos) {
            $t = [];
            foreach ($todos as $s) { $t[] = $s['info']['nombre'] . (is_int($s['loc']) ? ' en C' . $s['loc'] : ' en ' . $s['loc']); }
            $exp[] = 'Prefijos (' . count($todos) . '): ' . implode('; ', $t) . '.';
        }
        // 4. numeración
        $num = null;
        foreach ($fin as $d) {
            if ($clave($d) === $kb && $d['pos'] !== $pos) {
                [$s, $c] = $this->comparar($best, $d);
                if ($s === 0) { $num = 'la molécula es simétrica: las numeraciones posibles conducen al mismo nombre.'; }
                else { $num = $this->textoCriterio($c, $best, $d, true) . '.'; }
                break;
            }
        }
        if ($num !== null && !$res['omitir']) { $exp[] = 'Numeración: ' . $num; }
        elseif ($res['omitir']) { $exp[] = 'No se escriben localizadores porque no hay ambigüedad (P-14.3.4).'; }
        // 5. orden alfanumérico
        if (count($res['grupos']) > 1) {
            $exp[] = 'Los prefijos se citan en orden alfanumérico, sin considerar multiplicadores (di, tri, bis…) ni sec-/terc-: '
                . implode(' < ', array_keys($res['grupos'])) . ' (P-14.5).';
        } elseif ($res['grupos']) {
            $g = reset($res['grupos']);
            if (count($g['locs']) > 1) { $exp[] = 'Los prefijos idénticos se agrupan con multiplicadores (di, tri…; bis, tris… para prefijos compuestos) (P-16.3).'; }
        }
        // 6. nombres retenidos
        foreach ($res['notas'] as $n) {
            if ($n === 'retenido_benceno') { $exp[] = 'Se usa el nombre retenido del derivado del benceno (fenol, anilina, ácido benzoico, benzaldehído, benzamida, benzonitrilo), que es el nombre IUPAC preferido y admite sustitución; el grupo principal ocupa la posición 1 (P-62, P-63.1, P-65.1, P-66).'; }
            if ($n === 'retenido_acetico') { $exp[] = 'Para dos carbonos se usan los nombres retenidos preferidos: ácido acético, acetato, acetamida, acetonitrilo, acetaldehído (P-65.1, P-66).'; }
            if ($n === 'retenido_formico') { $exp[] = 'Con un solo carbono se usan los nombres retenidos preferidos: ácido fórmico, formiato, formamida, formaldehído, formonitrilo (P-65.1, P-66).'; }
            if ($n === 'retenido_oxalico') { $exp[] = 'HOOC–COOH conserva el nombre retenido preferido ácido oxálico (P-65.1.1).'; }
            if ($n === 'sistema_fus') { $exp[] = 'Los biciclos orto-fusionados con anillos de 5 o más miembros se nombran con el hidruro de fusión retenido (pentaleno, indeno, azuleno, naftaleno, heptaleno), que tiene una numeración fija: los átomos de fusión llevan letra (3a, 4a, 8a…) y solo se elige entre las orientaciones equivalentes (P-25.1, P-25.3).'; }
            if ($n === 'h_indicado') { $exp[] = 'El hidruro progenitor tiene un número impar de átomos: se indica con «xH-» qué carbono del anillo no forma parte de un doble enlace (hidrógeno indicado, P-14.7.1), con el localizador más bajo posible.'; }
            if ($n === 'h_anadido') { $exp[] = 'El sufijo «-ona» exige un carbono saturado que el hidruro progenitor no tiene: el hidrógeno añadido se cita entre paréntesis después del localizador del sufijo, p. ej., naftalen-1(2H)-ona (P-14.7.2).'; }
            if ($n === 'hidro') { $exp[] = 'Los carbonos saturados adicionales se expresan con prefijos «hidro» en número par (dihidro, tetrahidro…), citados justo antes del nombre del progenitor y no alfabetizados con los demás prefijos (P-31.1.4.2.4).'; }
            if ($n === 'sistema_vb') { $exp[] = 'Nomenclatura de von Baeyer: «biciclo» + número de átomos de cada puente en orden decreciente entre corchetes + nombre del alcano con el total de átomos. Se numera desde una cabeza de puente, recorriendo primero el puente más largo, luego el siguiente y al final el puente más corto (P-23.2).'; }
            if ($n === 'sistema_espiro') { $exp[] = 'Compuesto espiro: «espiro» + átomos de cada anillo, sin contar el átomo espiro, en orden creciente. Se numera empezando en el anillo menor, en el átomo contiguo al espiro, y se continúa por el átomo espiro hacia el anillo mayor (P-24.2).'; }
            if (strpos($n, 'retenido_hidrocarburo') === 0) { $exp[] = 'Nombre retenido preferido sin sustitución: tolueno, xileno (P-22.1.3), anisol (P-63.2) o acetileno (P-14.3.4.2).'; }
        }
        if (!empty($cl['halogeno'])) {
            $exp[] = 'Los halógenos solo se expresan como prefijos (fluoro, cloro, bromo, yodo) y se ordenan alfabéticamente con los demás prefijos (P-61.3.1, P-14.5).';
        }
        return array_values(array_filter($exp));
    }

    private function nivel(string $c): int
    {
        return array_search($c, ['', 'pn', 'anillo', 'len', 'mn', 'dn', 'plocs', 'mlocs', 'dlocs', 'cuenta', 'locs', 'citados', 'nombre'], true) ?: 0;
    }
}

/* =====================================================================
 *  Generador aleatorio de ejercicios
 * ===================================================================== */
final class GeneradorOrganico
{
    public const FAMILIAS = ['halogenos', 'alquenos', 'alquinos', 'aromaticos', 'biciclos', 'alcoholes', 'fenoles', 'eteres', 'aminas',
        'aldehidos', 'cetonas', 'acidos', 'esteres', 'amidas', 'nitrilos', 'nitro'];

    /** Esqueletos de dos anillos por nivel (fase 2): fusionados, con puente, espiro, ensamblajes y multiplicativos. */
    private const BICICLOS = [
        1 => ['c1ccc2ccccc2c1', 'C1Cc2ccccc2C1', 'C1CCc2ccccc2C1', 'C1CCC2CCCCC2C1', 'C1CC2CCC1C2', 'c1ccc(cc1)-c1ccccc1',
              'c1ccc(cc1)C1CCCCC1', 'c1ccc(Cc2ccccc2)cc1'],
        2 => ['c1ccc2ccccc2c1', 'C1C=Cc2ccccc21', 'C1CCC2(CC1)CCCC2', 'C1CC2CCC1CC2', 'C1CCC(CC1)C1CCCCC1', 'c1ccc(Oc2ccccc2)cc1',
              'c1ccc(CCc2ccccc2)cc1', 'C1=Cc2ccccc2CC1', 'C1CCC2CC2C1', 'C1CCC2(CC1)CCCC2', 'C1CC2CCC1C2'],
        3 => ['c1ccc2ccccc2c1', 'C1CC2CC(C1)C2', 'C1CC2CCC(C1)C2', 'c1ccc2c(c1)CC2', 'C1CC12CCCCC2', 'C1CCC2CCCC2C1', 'CC(C)(c1ccccc1)c1ccccc1',
              'c1ccc(-c2cccc3ccccc23)cc1', 'C1=CC2CCC1C2', 'C1CCC2(CC1)CCCCC2', 'O=C1CCCc2ccccc12', 'C1CC2CCCCC2C1'],
    ];

    private const RAMAS = [
        1 => ['C', 'C', 'CC'],
        2 => ['C', 'C', 'CC', 'C(C)C', 'CCC', 'CC(C)C', 'C(C)CC', 'C(C)(C)C'],
        3 => ['C', 'CC', 'C(C)C', 'CCC', 'CCCC', 'C(C)CC', 'CC(C)C', 'C(C)(C)C', 'CCC(C)C', 'CC(C)(C)C', 'C(C)CCC', 'C(CC)CC'],
    ];
    /** [cadena mín, máx], [ramas mín, máx], [anillo mín, máx], [grupos mín, máx] */
    private const NIVEL = [
        1 => ['cad' => [2, 6], 'ram' => [0, 1], 'ani' => [3, 6], 'grp' => [1, 1]],
        2 => ['cad' => [3, 8], 'ram' => [0, 2], 'ani' => [4, 7], 'grp' => [1, 2]],
        3 => ['cad' => [4, 12], 'ram' => [1, 3], 'ani' => [5, 8], 'grp' => [2, 3]],
    ];

    public static function generar(string $tipo, int $nivel, array $familias): MolOrg
    {
        $nivel = max(1, min(3, $nivel));
        $familias = array_values(array_intersect($familias, self::FAMILIAS));
        $aromatico = in_array('aromaticos', $familias, true) || in_array('fenoles', $familias, true);
        $cfg = self::NIVEL[$nivel];
        if ($tipo === 'aleatorio') {
            $op = ['lineal', 'ramificado', 'ramificado', 'ciclico', 'ciclico_ramificado'];
            $tipo = $op[random_int(0, count($op) - 1)];
        }
        $m = new MolOrg();
        $anillo = [];
        if (in_array('biciclos', $familias, true)) {
            $l = self::BICICLOS[$nivel];
            $m = MolOrg::desdeSmiles($l[random_int(0, count($l) - 1)]);
            $k = random_int(0, $nivel - 1);
            for ($i = 0; $i < $k; $i++) {
                $c = self::elegir(self::carbonos($m, 'any'));
                if ($c !== null) { self::injertar($m, $c, self::RAMAS[1][random_int(0, 2)]); }
            }
        } elseif ($tipo === 'ciclico' || $tipo === 'ciclico_ramificado') {
            if ($aromatico) {
                for ($i = 0; $i < 6; $i++) { $anillo[] = $m->nuevo('C', true); }
                for ($i = 0; $i < 6; $i++) { $m->enlazar($anillo[$i], $anillo[($i + 1) % 6], MolOrg::AROM); }
            } else {
                $N = random_int($cfg['ani'][0], $cfg['ani'][1]);
                for ($i = 0; $i < $N; $i++) { $anillo[] = $m->nuevo('C'); }
                for ($i = 0; $i < $N; $i++) { $m->enlazar($anillo[$i], $anillo[($i + 1) % $N]); }
            }
            $k = $tipo === 'ciclico' ? random_int(0, 1) : random_int(1, $nivel + 1);
            for ($i = 0; $i < $k; $i++) { self::injertar($m, $anillo[random_int(0, count($anillo) - 1)], self::rama($nivel)); }
        } else {
            $L = $tipo === 'lineal' ? random_int($cfg['cad'][0], $cfg['cad'][1]) : random_int(max(4, $cfg['cad'][0]), max(5, $cfg['cad'][1]));
            $at = [];
            for ($i = 0; $i < $L; $i++) { $at[] = $m->nuevo('C'); if ($i) { $m->enlazar($at[$i - 1], $at[$i]); } }
            if ($tipo === 'ramificado') {
                $k = random_int(max(1, $cfg['ram'][0]), max(1, $cfg['ram'][1]));
                for ($i = 0; $i < $k && $L > 2; $i++) { self::injertar($m, $at[random_int(1, $L - 2)], self::rama($nivel)); }
            }
            if ($aromatico) {                                   // fenilo unido a la cadena
                self::injertar($m, $at[random_int(0, $L - 1)], 'c1ccccc1');
            }
        }
        // grupos funcionales
        $grupos = array_values(array_diff($familias, ['aromaticos', 'biciclos']));
        if ($grupos) {
            [$gmin, $gmax] = $cfg['grp'];
            $n = random_int($gmin, $gmax);
            $elegidos = $grupos;
            shuffle($elegidos);
            $elegidos = array_slice($elegidos, 0, min($n, count($elegidos)));
            while (count($elegidos) < $n) { $elegidos[] = $grupos[random_int(0, count($grupos) - 1)]; }
            foreach ($elegidos as $f) { self::aplicar($m, $f, $nivel); }
        }
        return $m;
    }

    private static function rama(int $nivel): string
    {
        $r = self::RAMAS[$nivel];
        return $r[random_int(0, count($r) - 1)];
    }

    private static function injertar(MolOrg $m, int $ancla, string $smi): ?int
    {
        if ($m->hidrogenos($ancla) < 1) { return null; }
        $sub = MolOrg::desdeSmiles($smi);
        $map = [];
        for ($a = 0; $a < $sub->n; $a++) { $map[$a] = $m->nuevo($sub->el[$a], $sub->arom[$a], $sub->carga[$a]); }
        foreach ($sub->adj as $a => $vs) { foreach ($vs as $b => $o) { if ($a < $b) { $m->enlazar($map[$a], $map[$b], $o); } } }
        $m->enlazar($ancla, $map[0]);
        return $map[0];
    }

    /** Carbonos del esqueleto aptos para recibir un grupo. */
    private static function carbonos(MolOrg $m, string $modo): array
    {
        $res = [];
        $enAnillo = [];
        foreach ($m->sistemas() as $sx) { foreach ($sx['atomos'] as $a) { $enAnillo[$a] = true; } }
        foreach ($m->el as $a => $e) {
            if ($e !== 'C' || $m->hidrogenos($a) < 1) { continue; }
            $het = false; $multiple = false;
            foreach ($m->adj[$a] as $b => $o) {
                if ($m->el[$b] !== 'C') { $het = true; }
                if ($o === 2 || $o === 3) { $multiple = true; }
                // carbonos vecinos de grupos carbonílicos o de nitrilo: se evitan en los niveles bajos
            }
            if ($het) { continue; }
            if ($modo === 'arom' && !$m->arom[$a] && !(($enAnillo[$a] ?? false) && $multiple)) { continue; }
            if ($modo === 'sp3' && ($m->arom[$a] || $multiple)) { continue; }
            if ($modo === 'noarom' && $m->arom[$a]) { continue; }
            // no sobre carbonos unidos a carbonos funcionales (evita casos fuera de dominio)
            $res[] = $a;
        }
        return $res;
    }

    private static function elegir(array $l): ?int { return $l ? $l[random_int(0, count($l) - 1)] : null; }

    private static function aplicar(MolOrg $m, string $f, int $nivel): void
    {
        switch ($f) {
            case 'halogenos':
                $c = self::elegir(self::carbonos($m, 'any'));
                if ($c !== null) { $pool = ['Cl', 'Cl', 'Br', 'Br', 'F', 'I']; self::injertar($m, $c, $pool[random_int(0, 5)]); }
                return;
            case 'alquenos':
            case 'alquinos':
                $tri = $f === 'alquinos';
                $pares = [];
                foreach ($m->adj as $a => $vs) {
                    foreach ($vs as $b => $o) {
                        if ($a >= $b || $o !== 1 || $m->el[$a] !== 'C' || $m->el[$b] !== 'C' || $m->arom[$a] || $m->arom[$b]) { continue; }
                        $minH = $tri ? 2 : 1;
                        if ($m->hidrogenos($a) < $minH || $m->hidrogenos($b) < $minH) { continue; }
                        $ok = true;
                        foreach ([$a, $b] as $x) {
                            foreach ($m->adj[$x] as $y => $oy) { if ($oy >= 2 || $m->el[$y] !== 'C') { $ok = false; } }
                        }
                        if ($tri) {
                            $r = $m->anillos();
                            foreach ($r as $ring) { if (in_array($a, $ring, true) || in_array($b, $ring, true)) { $ok = false; } }
                        }
                        if ($ok) { $pares[] = [$a, $b]; }
                    }
                }
                if ($pares) { [$a, $b] = $pares[random_int(0, count($pares) - 1)]; $m->enlazar($a, $b, $tri ? 3 : 2); }
                return;
            case 'alcoholes':
                $c = self::elegir(self::carbonos($m, 'sp3'));
                if ($c !== null) { self::injertar($m, $c, 'O'); }
                return;
            case 'fenoles':
                $c = self::elegir(self::carbonos($m, 'arom'));
                if ($c !== null) { self::injertar($m, $c, 'O'); }
                return;
            case 'eteres':
                $c = self::elegir(self::carbonos($m, 'any'));
                $r = ['OC', 'OC', 'OCC', 'OC(C)C', 'OCCC'];
                if ($c !== null) { self::injertar($m, $c, $r[random_int(0, $nivel === 1 ? 1 : 4)]); }
                return;
            case 'aminas':
                $c = self::elegir(self::carbonos($m, 'sp3') ?: self::carbonos($m, 'arom'));
                $r = ['N', 'N', 'NC', 'N(C)C', 'NCC'];
                if ($c !== null) { self::injertar($m, $c, $r[random_int(0, $nivel === 1 ? 1 : 4)]); }
                return;
            case 'aldehidos':
                $c = self::elegir(self::carbonos($m, 'any'));
                if ($c !== null) { self::injertar($m, $c, 'C=O'); }
                return;
            case 'cetonas':
                $cand = [];
                foreach (self::carbonos($m, 'sp3') as $a) {
                    $nc = 0;
                    foreach ($m->adj[$a] as $b => $o) { if ($m->el[$b] === 'C') { $nc++; } }
                    if ($nc === 2 && $m->hidrogenos($a) === 2) { $cand[] = $a; }
                }
                $c = self::elegir($cand);
                if ($c !== null) { self::injertar($m, $c, 'O'); $o = $m->n - 1; $m->enlazar($c, $o, 2); }
                else { $c = self::elegir(self::carbonos($m, 'any')); if ($c !== null) { self::injertar($m, $c, 'C(C)=O'); } }
                return;
            case 'acidos':
                $c = self::elegir(self::carbonos($m, 'any'));
                if ($c !== null) { self::injertar($m, $c, 'C(=O)O'); }
                return;
            case 'esteres':
                $c = self::elegir(self::carbonos($m, 'any'));
                $r = ['C(=O)OC', 'C(=O)OCC', 'C(=O)OC(C)C', 'C(=O)OCCC'];
                if ($c !== null) { self::injertar($m, $c, $r[random_int(0, $nivel === 1 ? 1 : 3)]); }
                return;
            case 'amidas':
                $c = self::elegir(self::carbonos($m, 'any'));
                $r = ['C(N)=O', 'C(N)=O', 'C(=O)NC', 'C(=O)N(C)C'];
                if ($c !== null) { self::injertar($m, $c, $r[random_int(0, $nivel === 1 ? 1 : 3)]); }
                return;
            case 'nitrilos':
                $c = self::elegir(self::carbonos($m, 'any'));
                if ($c !== null) { self::injertar($m, $c, 'C#N'); }
                return;
            case 'nitro':
                $c = self::elegir(self::carbonos($m, 'arom') ?: self::carbonos($m, 'sp3'));
                if ($c !== null) {
                    $n = $m->nuevo('N', false, 1); $o1 = $m->nuevo('O'); $o2 = $m->nuevo('O', false, -1);
                    $m->enlazar($c, $n); $m->enlazar($n, $o1, 2); $m->enlazar($n, $o2);
                }
                return;
        }
    }
}
