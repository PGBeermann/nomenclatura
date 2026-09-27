<?php
/**
 * IupacAlcanos.php
 * Motor de generación de estructuras y nomenclatura IUPAC para alcanos
 * acíclicos (lineales y ramificados) y cicloalcanos monocíclicos con
 * sustituyentes alquilo.
 *
 * Fundamento normativo:
 *   Favre, H. A.; Powell, W. H. Nomenclature of Organic Chemistry. IUPAC
 *   Recommendations and Preferred Names 2013. Royal Society of Chemistry, 2014.
 *   DOI: 10.1039/9781849733069
 *     P-14.4  Numeración / regla de los localizadores más bajos
 *     P-14.5  Orden alfanumérico de prefijos
 *     P-16.3  Prefijos multiplicadores (di, tri… / bis, tris…)
 *     P-29    Prefijos sustituyentes (alquil / alcan-x-il); P-29.6 prefijos retenidos
 *     P-31.1.4 Criterios generales de numeración
 *     P-44.1.2.2 Los anillos son preferidos a las cadenas
 *     P-44.3  Selección de la cadena principal (la más larga)
 *     P-45.2  Criterios adicionales: más sustituyentes, localizadores más bajos,
 *             primer prefijo citado en orden alfanumérico
 *     P-45.5  Orden alfanumérico del nombre completo
 *
 * El motor no "memoriza" nombres: los calcula a partir del grafo molecular
 * (átomos de carbono y enlaces sencillos), por lo que es válido para cualquier
 * estructura generada dentro del dominio indicado (árbol o un solo anillo).
 *
 * Requiere PHP >= 7.4.  Autor: generado para UNACHI – Escuela de Química.
 */
declare(strict_types=1);

final class Molecula
{
    /** @var array<int,int[]> */
    public array $adj = [];
    /** @var array<int,string> símbolo del elemento: C, F, Cl, Br, I */
    public array $el = [];
    public int $n = 0;

    public const HALOGENOS = ['F', 'Cl', 'Br', 'I'];

    public function nuevoAtomo(string $el = 'C'): int
    {
        $this->adj[$this->n] = [];
        $this->el[$this->n] = $el;
        return $this->n++;
    }

    public function esCarbono(int $a): bool { return $this->el[$a] === 'C'; }

    public function numCarbonos(): int
    {
        return count(array_filter($this->el, fn($e) => $e === 'C'));
    }

    public function enlazar(int $a, int $b): void
    {
        $this->adj[$a][] = $b;
        $this->adj[$b][] = $a;
    }

    public function grado(int $a): int
    {
        return count($this->adj[$a]);
    }

    public function numEnlaces(): int
    {
        $s = 0;
        foreach ($this->adj as $v) { $s += count($v); }
        return intdiv($s, 2);
    }

    /**
     * Construye una molécula a partir de un SMILES de alcano (solo carbono, enlaces sencillos).
     * Acepta la salida de JSME: átomos C o [CH3], [CH2:4], etc. (los mapas se ignoran),
     * ramas, cierres de anillo (1-9, %nn) y marcas de estereoquímica (se ignoran).
     * Lanza NomenclaturaException con un mensaje didáctico si la estructura no es un alcano.
     */
    public static function desdeSmiles(string $smi): Molecula
    {
        $smi = trim($smi);
        if ($smi === '') { throw new NomenclaturaException('No hay ninguna estructura dibujada.'); }
        if (strlen($smi) > 400) { throw new NomenclaturaException('La estructura es demasiado grande.'); }
        $m = new Molecula();
        $pila = [];
        $prev = -1;
        $anillos = [];
        $len = strlen($smi);
        $agregar = function (string $el) use ($m, &$prev): void {
            $a = $m->nuevoAtomo($el);
            if ($prev >= 0) { $m->enlazar($prev, $a); }
            $prev = $a;
        };
        for ($i = 0; $i < $len; $i++) {
            $ch = $smi[$i];
            $sig = $smi[$i + 1] ?? '';
            if ($ch === 'C' && $sig === 'l') { $agregar('Cl'); $i++; }
            elseif ($ch === 'B' && $sig === 'r') { $agregar('Br'); $i++; }
            elseif ($ch === 'C') { $agregar('C'); }
            elseif ($ch === 'F' || $ch === 'I') { $agregar($ch); }
            elseif ($ch === '[') {
                $j = strpos($smi, ']', $i);
                if ($j === false) { throw new NomenclaturaException('SMILES mal formado.'); }
                $atomo = substr($smi, $i, $j - $i + 1);
                if (preg_match('/^\[C(@{1,2})?(H[0-4]?)?(:\d+)?\]$/', $atomo)) { $agregar('C'); }
                elseif (preg_match('/^\[(F|Cl|Br|I)(:\d+)?\]$/', $atomo, $mm)) { $agregar($mm[1]); }
                else {
                    if (preg_match('/[+\-]/', $atomo)) { throw new NomenclaturaException('La estructura tiene cargas: no es un alcano ni un haloalcano.'); }
                    if (preg_match('/^\[\d/', $atomo)) { throw new NomenclaturaException('No se admiten isótopos.'); }
                    throw new NomenclaturaException('La estructura contiene átomos que no son carbono ni halógeno (' . $atomo . ').');
                }
                $i = $j;
            } elseif ($ch === '(') {
                $pila[] = $prev;
            } elseif ($ch === ')') {
                if (!$pila) { throw new NomenclaturaException('SMILES mal formado.'); }
                $prev = array_pop($pila);
            } elseif (ctype_digit($ch) || $ch === '%') {
                $d = $ch;
                if ($ch === '%') { $d = substr($smi, $i + 1, 2); $i += 2; }
                if (isset($anillos[$d])) {
                    if ($anillos[$d] === $prev || in_array($prev, $m->adj[$anillos[$d]], true)) {
                        throw new NomenclaturaException('SMILES mal formado.');
                    }
                    $m->enlazar($anillos[$d], $prev); unset($anillos[$d]);
                } else { $anillos[$d] = $prev; }
            } elseif ($ch === '-' || $ch === '/' || $ch === '\\') {
                // enlace sencillo explícito o marca estérea: se ignora
            } elseif ($ch === '=' || $ch === '#' || $ch === '$') {
                throw new NomenclaturaException('La estructura tiene enlaces dobles o triples: es un alqueno/alquino, no un alcano.');
            } elseif ($ch === '.') {
                throw new NomenclaturaException('Hay más de una molécula en el editor. Dibuje una sola estructura.');
            } elseif (ctype_lower($ch)) {
                throw new NomenclaturaException('La estructura contiene un anillo aromático: no es un alcano ni un cicloalcano.');
            } else {
                throw new NomenclaturaException('La estructura contiene átomos que no son carbono ni halógeno (' . $ch . ').');
            }
        }
        if ($anillos || $pila) { throw new NomenclaturaException('SMILES mal formado.'); }
        foreach ($m->adj as $a => $v) {
            if ($m->el[$a] === 'C' && count($v) > 4) { throw new NomenclaturaException('Un carbono tiene más de cuatro enlaces.'); }
            if ($m->el[$a] !== 'C') {
                if (count($v) !== 1 || $m->el[$v[0]] !== 'C') {
                    throw new NomenclaturaException('Cada halógeno debe estar unido a un solo átomo de carbono.');
                }
            }
        }
        if ($m->numCarbonos() === 0) { throw new NomenclaturaException('La estructura no contiene carbono.'); }
        return $m;
    }

    /** Longitud (en átomos) del camino más largo del grafo acíclico (diámetro). */
    public function diametro(): int
    {
        if ($this->n <= 1) { return $this->n; }
        $bfs = function (int $s): array {
            $dist = [$s => 1]; $cola = [$s]; $ult = $s;
            while ($cola) {
                $a = array_shift($cola); $ult = $a;
                foreach ($this->adj[$a] as $b) {
                    if (!isset($dist[$b])) { $dist[$b] = $dist[$a] + 1; $cola[] = $b; }
                }
            }
            return [$ult, $dist[$ult]];
        };
        [$u] = $bfs(0);
        return $bfs($u)[1];
    }
}

/** Error de dominio con mensaje apto para mostrar al estudiante. */
final class NomenclaturaException extends RuntimeException {}

final class IupacAlcanos
{
    private const RAIZ = [
        1 => 'met', 2 => 'et', 3 => 'prop', 4 => 'but', 5 => 'pent', 6 => 'hex',
        7 => 'hept', 8 => 'oct', 9 => 'non', 10 => 'dec', 11 => 'undec', 12 => 'dodec',
        13 => 'tridec', 14 => 'tetradec', 15 => 'pentadec', 16 => 'hexadec',
        17 => 'heptadec', 18 => 'octadec', 19 => 'nonadec', 20 => 'icos',
    ];
    private const RAIZ_EN = [
        1 => 'meth', 2 => 'eth', 3 => 'prop', 4 => 'but', 5 => 'pent', 6 => 'hex',
        7 => 'hept', 8 => 'oct', 9 => 'non', 10 => 'dec', 11 => 'undec', 12 => 'dodec',
        13 => 'tridec', 14 => 'tetradec', 15 => 'pentadec', 16 => 'hexadec',
        17 => 'heptadec', 18 => 'octadec', 19 => 'nonadec', 20 => 'icos',
    ];
    private const MULT_SIMPLE = [2 => 'di', 3 => 'tri', 4 => 'tetra', 5 => 'penta', 6 => 'hexa', 7 => 'hepta', 8 => 'octa',
        9 => 'nona', 10 => 'deca', 11 => 'undeca', 12 => 'dodeca', 13 => 'trideca', 14 => 'tetradeca', 15 => 'pentadeca',
        16 => 'hexadeca', 17 => 'heptadeca', 18 => 'octadeca', 19 => 'nonadeca', 20 => 'icosa'];
    private const MULT_COMPLEJO = [2 => 'bis', 3 => 'tris', 4 => 'tetrakis', 5 => 'pentakis', 6 => 'hexakis', 7 => 'heptakis',
        8 => 'octakis', 9 => 'nonakis', 10 => 'decakis'];  // P-16.3.6

    private const HALO = [
        'es' => ['F' => 'fluoro', 'Cl' => 'cloro', 'Br' => 'bromo', 'I' => 'yodo'],
        'en' => ['F' => 'fluoro', 'Cl' => 'chloro', 'Br' => 'bromo', 'I' => 'iodo'],
    ];
    private const HALURO = [
        'es' => ['F' => 'fluoruro', 'Cl' => 'cloruro', 'Br' => 'bromuro', 'I' => 'yoduro'],
        'en' => ['F' => 'fluoride', 'Cl' => 'chloride', 'Br' => 'bromide', 'I' => 'iodide'],
    ];

    private Molecula $mol;  // molécula completa (C + halógenos)
    private Molecula $m;    // esqueleto carbonado (índices propios)
    /** @var int[] índice en el esqueleto → índice en la molécula completa */
    private array $orig = [];
    /** @var array<int,string[]> halógenos unidos a cada carbono del esqueleto */
    private array $hal = [];
    /** Si es true, se ignoran los halógenos (para el nombre radicofuncional R–X). */
    public bool $sinHalogenos = false;
    private string $lang;   // 'es' | 'en'
    private string $estilo; // 'pin' (nombre IUPAC preferido 2013) | 'trad' (prefijos retenidos: isopropil…) | 'sis79' (sistemático 1979: 1-metiletil…)
    private array $memo = [];

    public function __construct(Molecula $mol, string $lang = 'es', string $estilo = 'pin')
    {
        $this->mol = $mol;
        $this->lang = $lang;
        $this->estilo = $estilo;
        // Esqueleto carbonado: la cadena principal y los grupos alquilo se buscan solo entre
        // carbonos; los halógenos son siempre prefijos sustituyentes (P-61.3.1).
        $nuevo = [];
        $esq = new Molecula();
        foreach ($mol->el as $a => $e) {
            if ($e === 'C') { $nuevo[$a] = $esq->nuevoAtomo('C'); $this->orig[$nuevo[$a]] = $a; $this->hal[$nuevo[$a]] = []; }
        }
        foreach ($mol->adj as $a => $vs) {
            foreach ($vs as $b) {
                if ($mol->el[$a] === 'C' && $mol->el[$b] === 'C' && $a < $b) { $esq->enlazar($nuevo[$a], $nuevo[$b]); }
                elseif ($mol->el[$a] === 'C' && $mol->el[$b] !== 'C') { $this->hal[$nuevo[$a]][] = $mol->el[$b]; }
            }
        }
        $this->m = $esq;
    }

    /** Halógenos del carbono $a (vacío si se ignoran). */
    private function halDe(int $a): array
    {
        return $this->sinHalogenos ? [] : $this->hal[$a];
    }

    private function totalHalogenos(): int
    {
        $t = 0;
        foreach ($this->hal as $h) { $t += count($h); }
        return $t;
    }

    private function infoHalogeno(string $x): array
    {
        $n = self::HALO[$this->lang][$x];
        return ['nombre' => $n, 'tipo' => 'simple', 'clave' => self::claveAlfa($n), 'longitud' => 0, 'hal' => true];
    }

    /* ======================= utilidades léxicas ======================= */

    private function raiz(int $n): string
    {
        if ($n > 20) { throw new NomenclaturaException("La cadena principal tiene $n carbonos; el programa admite hasta 20 (icosano)."); }
        return $this->lang === 'es' ? self::RAIZ[$n] : self::RAIZ_EN[$n];
    }

    private function yl(): string { return $this->lang === 'es' ? 'il' : 'yl'; }
    private function ano(): string { return $this->lang === 'es' ? 'ano' : 'ane'; }

    /** Clave de ordenación alfanumérica (P-14.5): se ignoran prefijos en cursiva (sec-, terc-/tert-),
     *  localizadores, signos y multiplicadores externos; se conservan los internos de prefijos complejos. */
    private static function claveAlfa(string $nombre): string
    {
        $s = preg_replace('/\b(tert|terc|sec)-/', '', $nombre);
        return preg_replace('/[^a-z]/', '', strtolower($s));
    }

    /* ======================= topología ======================= */

    /** Devuelve los átomos del ciclo en orden (vacío si la molécula es acíclica). */
    public function ciclo(): array
    {
        $m = $this->m;
        $r = $m->numEnlaces() - $m->n + 1;
        if ($r === 0) { return []; }
        if ($r > 1) { throw new NomenclaturaException('La estructura tiene más de un anillo. El programa nombra alcanos acíclicos y cicloalcanos con un solo anillo.'); }
        $grado = [];
        foreach ($m->adj as $a => $v) { $grado[$a] = count($v); }
        $vivo = array_fill(0, $m->n, true);
        $cola = [];
        foreach ($grado as $a => $g) { if ($g <= 1) { $cola[] = $a; } }
        while ($cola) {
            $a = array_pop($cola);
            if (!$vivo[$a]) { continue; }
            $vivo[$a] = false;
            foreach ($m->adj[$a] as $b) {
                if ($vivo[$b] && --$grado[$b] === 1) { $cola[] = $b; }
            }
        }
        $inicio = array_search(true, $vivo, true);
        $orden = [$inicio];
        $prev = -1; $act = $inicio;
        while (true) {
            $sig = -1;
            foreach ($m->adj[$act] as $b) {
                if ($vivo[$b] && $b !== $prev) { $sig = $b; break; }
            }
            if ($sig === $inicio || $sig === -1) { break; }
            $orden[] = $sig; $prev = $act; $act = $sig;
        }
        return $orden;
    }

    /** Átomos del subárbol que cuelga de $r cuando se llega desde $p. */
    private function subarbol(int $r, int $p): array
    {
        $vis = [$r => true, $p => true];
        $pila = [$r];
        $res = [];
        while ($pila) {
            $a = array_pop($pila);
            $res[] = $a;
            foreach ($this->m->adj[$a] as $b) {
                if (!isset($vis[$b])) { $vis[$b] = true; $pila[] = $b; }
            }
        }
        return $res;
    }

    /** Camino único entre u y v dentro de un árbol (restringido al conjunto $permitidos). */
    private function camino(int $u, int $v, array $permitidos): array
    {
        $padre = [$u => -1];
        $cola = [$u];
        while ($cola) {
            $a = array_shift($cola);
            if ($a === $v) { break; }
            foreach ($this->m->adj[$a] as $b) {
                if (isset($permitidos[$b]) && !array_key_exists($b, $padre)) {
                    $padre[$b] = $a; $cola[] = $b;
                }
            }
        }
        $ruta = [];
        for ($x = $v; $x !== -1; $x = $padre[$x]) { $ruta[] = $x; }
        return array_reverse($ruta);
    }

    /* ======================= sustituyentes ======================= */

    /**
     * Nombra el sustituyente que nace en $r (unido a $p).
     * Devuelve ['nombre','tipo' => simple|compuesto|complejo,'clave','atomos'=>n]
     */
    public function sustituyente(int $r, int $p): array
    {
        $k = "$r|$p";
        if (isset($this->memo[$k])) { return $this->memo[$k]; }
        $atomos = $this->subarbol($r, $p);
        $n = count($atomos);
        $yl = $this->yl();

        $conHal = false;
        foreach ($atomos as $a) { if ($this->halDe($a)) { $conHal = true; break; } }

        if ($n === 1) {
            if (!$conHal) { return $this->memo[$k] = $this->info($this->raiz(1) . $yl, 'simple', 1); }
            // metilo halogenado: clorometil, trifluorometil…
            $d = $this->evaluarCadena([$r], $p, false, true);
            return $this->memo[$k] = $this->info($d['prefijos'] . $this->raiz(1) . $yl, 'complejo', 1);
        }
        // terc-butilo: C(CH3)3. Prefijo retenido y preferido en 2013 (P-29.6), solo sin sustituir;
        // en el estilo sistemático de 1979 se escribe 1,1-dimetiletil.
        if ($this->estilo !== 'sis79' && !$conHal && $n === 4 && count(array_diff($this->m->adj[$r], [$p])) === 3) {
            return $this->memo[$k] = $this->info(($this->lang === 'es' ? 'terc-but' : 'tert-but') . $yl, 'simple', 3);
        }

        // PIN 2013 (P-29.2): cadena más larga que contenga el átomo de unión; la valencia
        // libre recibe el localizador más bajo posible (propan-2-il, pentan-3-il…).
        // IUPAC 1979 (A-2.6): cadena más larga que EMPIECE en el átomo de unión (C1): 1-metiletil…
        $desdeUnion = $this->estilo !== 'pin';
        $perm = array_fill_keys($atomos, true);
        $extremos = [$r];
        foreach ($atomos as $a) { if ($a !== $r && $this->m->grado($a) === 1) { $extremos[] = $a; } }

        $cand = [];
        $maxLen = 0;
        foreach ($extremos as $u) {
            if ($desdeUnion && $u !== $r) { continue; }
            foreach ($extremos as $v) {
                if ($u === $v) { continue; }
                $ruta = $this->camino($u, $v, $perm);
                $pos = array_search($r, $ruta, true);
                if ($pos === false) { continue; }
                $L = count($ruta);
                if ($L < $maxLen) { continue; }
                if ($L > $maxLen) { $maxLen = $L; $cand = []; }
                $d = $this->evaluarCadena($ruta, $p, false, true);
                $d['fv'] = $pos + 1;
                $cand[] = $d;
            }
        }
        usort($cand, function ($a, $b) {
            if ($a['fv'] !== $b['fv']) { return $a['fv'] <=> $b['fv']; } // valencia libre: localizador más bajo
            return self::comparar($a, $b)[0];
        });
        $best = $cand[0];
        $raiz = $this->raiz($maxLen);
        $base = $best['fv'] === 1 ? $raiz . $yl : $raiz . 'an-' . $best['fv'] . '-' . $yl;
        $nombre = $best['prefijos'] . $base;

        // Prefijos retenidos no sustituidos (IUPAC 1979, A-2.25; 1993, R-9.1)
        if ($this->estilo === 'trad') {
            $ret = $this->lang === 'es'
                ? ['1-metiletil' => 'isopropil', '1-metilpropil' => 'sec-butil', '2-metilpropil' => 'isobutil',
                   '3-metilbutil' => 'isopentil', '2,2-dimetilpropil' => 'neopentil', '1,1-dimetilpropil' => 'terc-pentil']
                : ['1-methylethyl' => 'isopropyl', '1-methylpropyl' => 'sec-butyl', '2-methylpropyl' => 'isobutyl',
                   '3-methylbutyl' => 'isopentyl', '2,2-dimethylpropyl' => 'neopentyl', '1,1-dimethylpropyl' => 'tert-pentyl'];
            if (isset($ret[$nombre])) { return $this->memo[$k] = $this->info($ret[$nombre], 'simple', $maxLen); }
        }
        if ($best['cuenta'] === 0) {
            $tipo = $best['fv'] === 1 ? 'simple' : 'compuesto';
            return $this->memo[$k] = $this->info($nombre, $tipo, $maxLen);
        }
        return $this->memo[$k] = $this->info($nombre, 'complejo', $maxLen);
    }

    private function info(string $nombre, string $tipo, int $longitud): array
    {
        return ['nombre' => $nombre, 'tipo' => $tipo, 'clave' => self::claveAlfa($nombre), 'longitud' => $longitud];
    }

    /* ======================= evaluación de cadenas ======================= */

    /**
     * Evalúa una cadena orientada (lista de átomos; el índice 0 recibe el localizador 1).
     * $excluir: átomo que no debe considerarse sustituyente (unión al padre), o -1.
     */
    private function evaluarCadena(array $cadena, int $excluir = -1, bool $esAnillo = false, bool $esSust = false): array
    {
        $enCadena = array_fill_keys($cadena, true);
        $subs = [];
        foreach ($cadena as $i => $a) {
            foreach ($this->m->adj[$a] as $b) {
                if ($b === $excluir || isset($enCadena[$b])) { continue; }
                $s = $this->sustituyente($b, $a);
                $subs[] = ['loc' => $i + 1, 'info' => $s, 'atomo' => $b, 'desde' => $a];
            }
            // halógenos: siempre prefijos (fluoro, cloro, bromo, yodo)
            foreach ($this->halDe($a) as $x) {
                $subs[] = ['loc' => $i + 1, 'info' => $this->infoHalogeno($x), 'atomo' => -1, 'desde' => $a];
            }
        }
        return $this->armar($subs, $cadena, $esAnillo, $esSust);
    }

    private function armar(array $subs, array $cadena, bool $esAnillo, bool $esSust = false): array
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
        $partes = [];
        $citados = [];
        $locs = [];
        // Omisión de localizadores (P-14.3.4):
        //  · cicloalcano monosustituido (metilciclohexano)
        //  · metano / grupo metilo (diclorometano, trifluorometil)
        //  · etano monosustituido (cloroetano)
        //  · sustitución total por un mismo prefijo (hexacloroetano)
        $nC = count($cadena);
        $hDisp = $esAnillo ? 2 * $nC : 2 * $nC + 2;
        $nGrupos = count(array_unique(array_map(fn($s) => $s['info']['nombre'], $subs)));
        $omitirLoc = ($esAnillo && count($subs) === 1)
            || $nC === 1
            || (!$esSust && !$esAnillo && $nC === 2 && count($subs) === 1)
            || (!$esSust && count($subs) === $hDisp && $nGrupos === 1);
        foreach ($grupos as $nm => $g) {
            sort($g['locs']);
            $cnt = count($g['locs']);
            $info = $g['info'];
            foreach ($g['locs'] as $l) { $citados[] = $l; $locs[] = $l; }
            $paren = $info['tipo'] !== 'simple';
            if ($cnt === 1) { $mult = ''; }
            elseif ($info['tipo'] === 'complejo') { $mult = self::MULT_COMPLEJO[$cnt]; }
            else { $mult = self::MULT_SIMPLE[$cnt]; }
            // di + prefijo simple que comienza en "i" (p. ej. diisopropil): se escribe sin elisión.
            // Ante prefijos en cursiva se usa guion: di-terc-butil, di-sec-butil.
            if ($mult !== '' && !$paren && preg_match('/^(sec|terc|tert)-/', $nm)) { $mult .= '-'; }
            // signos de inclusión anidados: ( ) → [ ] → { }  (P-16.5.4)
            if ($paren) {
                if (strpos($nm, '[') !== false)      { $env = '{' . $nm . '}'; }
                elseif (strpos($nm, '(') !== false)  { $env = '[' . $nm . ']'; }
                else                                 { $env = '(' . $nm . ')'; }
            }
            $txt = $mult . ($paren ? $env : $nm);
            $partes[] = ($omitirLoc ? '' : implode(',', $g['locs']) . '-') . $txt;
        }
        sort($locs);
        $pref = implode('-', $partes);
        return [
            'prefijos' => $pref,
            'cuenta'   => count($subs),
            'locs'     => $locs,
            'citados'  => $citados,
            'cadena'   => $cadena,
            'subs'     => $subs,
            'grupos'   => $grupos,
        ];
    }

    /**
     * Compara dos candidatos según P-45.2 / P-31.1.4 / P-45.5.
     * Devuelve [signo, nivel]; nivel 1=nº sustituyentes, 2=conjunto de localizadores,
     * 3=primer prefijo citado, 4=orden alfanumérico, 0=idénticos.
     */
    private static function comparar(array $a, array $b): array
    {
        if ($a['cuenta'] !== $b['cuenta']) { return [$b['cuenta'] <=> $a['cuenta'], 1]; }
        $c = self::cmpLocs($a['locs'], $b['locs']);
        if ($c !== 0) { return [$c, 2]; }
        $c = self::cmpLocs($a['citados'], $b['citados']);
        if ($c !== 0) { return [$c, 3]; }
        $c = strcmp($a['nombreCmp'] ?? $a['prefijos'], $b['nombreCmp'] ?? $b['prefijos']);
        if ($c !== 0) { return [$c < 0 ? -1 : 1, 4]; }
        return [0, 0];
    }

    private static function cmpLocs(array $x, array $y): int
    {
        $n = min(count($x), count($y));
        for ($i = 0; $i < $n; $i++) {
            if ($x[$i] !== $y[$i]) { return $x[$i] <=> $y[$i]; }
        }
        return count($x) <=> count($y);
    }

    /* ======================= nombre completo ======================= */

    /**
     * Devuelve:
     *  nombre, cadena (átomos en orden de localizador), anillo(bool), explicacion(array de strings),
     *  subs (sustituyentes con localizador)
     */
    public function nombrar(): array
    {
        $r = $this->nombrarEsqueleto();
        $nHal = $this->totalHalogenos();
        if ($nHal > 0 && !$this->sinHalogenos) {
            $r['explicacion'][] = 'Los halógenos se nombran únicamente como prefijos (fluoro, cloro, bromo, yodo), no como sufijos; '
                . 'cuentan como sustituyentes al elegir la cadena y al numerar, y se ordenan alfabéticamente junto con los grupos alquilo (P-61.3.1, P-14.5).';
        }
        $r['radicofuncional'] = $nHal === 1 ? $this->nombreRadicofuncional() : null;
        // índices del esqueleto → índices de la molécula completa
        $r['cadena'] = array_map(fn($a) => $this->orig[$a], $r['cadena']);
        return $r;
    }

    /**
     * Nombre de clase funcional (radicofuncional) R–X para monohaloalcanos:
     * «cloruro de isopropilo», «isopropyl chloride». Aceptado en nomenclatura general,
     * no es nombre preferido (P-61.3.2). Devuelve null si no es aplicable.
     */
    public function nombreRadicofuncional(): ?string
    {
        $c = null; $x = null;
        foreach ($this->hal as $a => $h) { if ($h) { $c = $a; $x = $h[0]; } }
        if ($c === null) { return null; }
        $sin = new IupacAlcanos($this->mol, $this->lang, $this->estilo);
        $sin->sinHalogenos = true;
        if ($sin->ciclo()) {
            // solo cicloalquilo sin otros sustituyentes: cloruro de ciclohexilo
            if ($this->m->n !== count($sin->ciclo())) { return null; }
            $R = ($this->lang === 'es' ? 'ciclo' : 'cyclo') . $this->raiz($this->m->n) . $this->yl();
        } else {
            $info = $sin->sustituyente($c, -1);
            // solo para grupos R sencillos (metilo, isopropilo, sec-butilo, propan-2-ilo…);
            // con grupos R complejos el nombre radicofuncional pierde utilidad didáctica
            if ($info['tipo'] === 'complejo') { return null; }
            $R = $info['nombre'];
        }
        $haluro = self::HALURO[$this->lang][$x];
        return $this->lang === 'es' ? $haluro . ' de ' . $R . 'o' : $R . ' ' . $haluro;
    }

    private function nombrarEsqueleto(): array
    {
        $m = $this->m;
        if ($m->n === 0) { throw new NomenclaturaException('No hay ninguna estructura dibujada.'); }
        if ($m->n > 60) { throw new NomenclaturaException('La estructura tiene más de 60 carbonos.'); }
        $ano = $this->ano();
        $es = $this->lang === 'es';
        $ciclo = $this->ciclo();
        $exp = [];

        if ($ciclo) {
            $N = count($ciclo);
            $padre = ($es ? 'ciclo' : 'cyclo') . $this->raiz($N) . $ano;
            $cands = [];
            for ($s = 0; $s < $N; $s++) {
                foreach ([1, -1] as $dir) {
                    $orden = [];
                    for ($i = 0; $i < $N; $i++) { $orden[] = $ciclo[(($s + $dir * $i) % $N + $N) % $N]; }
                    $d = $this->evaluarCadena($orden, -1, true);
                    $d['nombreCmp'] = $d['prefijos'] . $padre;
                    $cands[] = $d;
                }
            }
            usort($cands, fn($a, $b) => self::comparar($a, $b)[0]);
            $best = $cands[0];
            $nombre = $best['prefijos'] . $padre;
            $total = $m->n;
            $exp[] = "Hidruro progenitor: anillo de $N carbonos ($padre). Según IUPAC 2013 (P-44.1.2.2) el anillo es preferido a las cadenas como estructura principal.";
            $larga = 0; $subLarga = '';
            foreach ($best['subs'] as $s) {
                if ($s['info']['longitud'] > $larga) { $larga = $s['info']['longitud']; $subLarga = $s['info']['nombre']; }
            }
            if ($larga > $N) {
                $exp[] = "Nota: el sustituyente $subLarga contiene una cadena de $larga carbonos, más larga que el anillo ($N C). Las reglas IUPAC de 1979/1993 tomarían la cadena como progenitor (el anillo sería un sustituyente cicloalquil); las recomendaciones de 2013 prefieren siempre el anillo.";
            }
            if ($best['cuenta'] > 0) {
                $exp[] = $this->describirSustituyentes($best);
                if ($best['cuenta'] === 1) {
                    $exp[] = 'Cicloalcano monosustituido: el localizador 1 se omite (P-14.3.4).';
                } elseif ($best['prefijos'] !== '' && strpos($best['prefijos'], '-') === false) {
                    $exp[] = 'Todos los hidrógenos están sustituidos por el mismo prefijo: se omiten los localizadores (P-14.3.4).';
                } else {
                    $exp[] = $this->describirNumeracion($best, $cands, true);
                }
            } else {
                $exp[] = 'Cicloalcano sin sustituyentes.';
            }
            $exp[] = $this->describirOrdenAlfa($best);
            return ['nombre' => $nombre, 'cadena' => $best['cadena'], 'anillo' => true, 'explicacion' => array_values(array_filter($exp)), 'subs' => $best['subs']];
        }

        // ---------- acíclico ----------
        if ($m->n === 1) {
            $d = $this->evaluarCadena([0]);
            $exp = ['Un solo átomo de carbono: hidruro progenitor ' . $this->raiz(1) . $ano . '.'];
            if ($d['cuenta'] > 0) {
                $exp[] = $this->describirSustituyentes($d);
                $exp[] = 'En el metano todas las posiciones son equivalentes: no se escriben localizadores (P-14.3.4).';
                $exp[] = $this->describirOrdenAlfa($d);
            }
            return ['nombre' => $d['prefijos'] . $this->raiz(1) . $ano, 'cadena' => [0], 'anillo' => false,
                'explicacion' => array_values(array_filter($exp)), 'subs' => $d['subs']];
        }
        $hojas = [];
        foreach ($m->adj as $a => $v) { if (count($v) === 1) { $hojas[] = $a; } }
        $todos = array_fill_keys(range(0, $m->n - 1), true);
        $cands = [];
        $maxLen = 0;
        foreach ($hojas as $u) {
            foreach ($hojas as $v) {
                if ($u === $v) { continue; }
                $ruta = $this->camino($u, $v, $todos);
                $L = count($ruta);
                if ($L < $maxLen) { continue; }
                if ($L > $maxLen) { $maxLen = $L; $cands = []; }
                $cands[] = $ruta;
            }
        }
        $padre = $this->raiz($maxLen) . $ano;
        $datos = [];
        foreach ($cands as $ruta) {
            $d = $this->evaluarCadena($ruta);
            $d['nombreCmp'] = $d['prefijos'] . $padre;
            $datos[] = $d;
        }
        usort($datos, fn($a, $b) => self::comparar($a, $b)[0]);
        $best = $datos[0];
        $nombre = $best['prefijos'] . $padre;

        // número de cadenas distintas (como conjunto de átomos) de longitud máxima
        $conj = [];
        foreach ($datos as $d) { $c = $d['cadena']; sort($c); $conj[implode(',', $c)] = true; }
        $nCad = count($conj);
        if ($best['cuenta'] === 0) {
            $exp[] = "Alcano lineal (no ramificado) de $maxLen carbonos: $padre.";
        } else {
            $exp[] = "Cadena principal: la cadena continua más larga, de $maxLen carbonos → $padre (P-44.3).";
            if ($nCad > 1) {
                $exp[] = $this->describirEleccionCadena($best, $datos, $nCad);
            }
            $exp[] = $this->describirSustituyentes($best);
            if (strpos($best['prefijos'], '-') === false) {
                $exp[] = $maxLen === 2 && $best['cuenta'] === 1
                    ? 'Etano monosustituido: ambas posiciones son equivalentes y el localizador se omite (P-14.3.4).'
                    : 'Todos los hidrógenos están sustituidos por el mismo prefijo: se omiten los localizadores (P-14.3.4).';
            } else {
                $exp[] = $this->describirNumeracion($best, $datos, false);
            }
            $exp[] = $this->describirOrdenAlfa($best);
        }
        return ['nombre' => $nombre, 'cadena' => $best['cadena'], 'anillo' => false, 'explicacion' => array_values(array_filter($exp)), 'subs' => $best['subs']];
    }

    /* ======================= explicación didáctica ======================= */

    private static function fmtLocs(array $l): string { return '{' . implode(',', $l) . '}'; }

    private function describirEleccionCadena(array $best, array $datos, int $nCad): string
    {
        $cb = $best['cadena']; sort($cb); $kb = implode(',', $cb);
        // rival más cercano: cadena distinta cuyo nombre difiere y que empata en más criterios
        // ($datos ya está ordenado: la primera aparición de cada cadena es su mejor numeración)
        $d = null; $nivel = 0; $vistos = [$kb => true];
        foreach ($datos as $x) {
            $c = $x['cadena']; sort($c); $kc = implode(',', $c);
            if (isset($vistos[$kc])) { continue; }
            $vistos[$kc] = true;
            [, $nv] = self::comparar($best, $x);
            if ($nv === 0) { continue; }
            if ($d === null || $nv > $nivel) { $d = $x; $nivel = $nv; }
        }
        {
            switch ($nivel) {
                case 1: return "Hay $nCad cadenas de igual longitud; se elige la que tiene el mayor número de sustituyentes ({$best['cuenta']} frente a {$d['cuenta']}) (P-45.2.1).";
                case 2: return "Hay $nCad cadenas de igual longitud y con igual número de sustituyentes; se elige la que tiene los localizadores más bajos: " . self::fmtLocs($best['locs']) . ' frente a ' . self::fmtLocs($d['locs']) . ' (P-45.2.2).';
                case 3: return "Hay $nCad cadenas equivalentes en longitud, número y localizadores; decide el localizador más bajo para el prefijo citado primero en orden alfanumérico (P-45.2.3).";
                case 4: return "Hay $nCad cadenas equivalentes; se elige la que conduce al nombre que aparece primero en orden alfanumérico (P-45.5).";
                default: return "Hay $nCad cadenas de longitud máxima, pero todas conducen al mismo nombre (son equivalentes por simetría).";
            }
        }
        return '';
    }

    private function describirSustituyentes(array $best): string
    {
        $t = [];
        foreach ($best['subs'] as $s) { $t[] = $s['info']['nombre'] . ' en C' . $s['loc']; }
        return 'Sustituyentes (' . $best['cuenta'] . '): ' . implode('; ', $t) . '.';
    }

    private function describirNumeracion(array $best, array $datos, bool $anillo): string
    {
        $mismaCad = function ($d) use ($best, $anillo) {
            if ($anillo) { return true; }
            return $d['cadena'] === array_reverse($best['cadena']);
        };
        $rival = null; $nivelMax = -1;
        foreach ($datos as $d) {
            if ($d === $best || !$mismaCad($d)) { continue; }
            if ($d['cadena'] === $best['cadena']) { continue; }
            [, $nivel] = self::comparar($best, $d);
            // el rival más "cercano" es el que empata en más criterios
            if ($nivel === 0) { $rival = $d; $nivelMax = 0; break; }
            if ($nivel > $nivelMax) { $nivelMax = $nivel; $rival = $d; }
        }
        $pref = $anillo ? 'Numeración del anillo' : 'Numeración de la cadena';
        if ($rival === null) { return ''; }
        switch ($nivelMax) {
            case 0: return "$pref: la molécula es simétrica; ambas numeraciones conducen al mismo nombre.";
            case 2: return "$pref: se elige el conjunto de localizadores más bajo " . self::fmtLocs($best['locs']) . ' frente a ' . self::fmtLocs($rival['locs']) . ' (primer punto de diferencia, P-14.4 y P-31.1.4).';
            case 3:
                $g = reset($best['grupos']);
                return "$pref: los conjuntos de localizadores coinciden " . self::fmtLocs($best['locs']) . "; el localizador más bajo se asigna al sustituyente citado primero en orden alfanumérico ({$g['info']['nombre']}) (P-31.1.4, P-14.5).";
            case 4: return "$pref: se elige la numeración que da el nombre primero en orden alfanumérico (P-45.5).";
            default: return "$pref: localizadores más bajos " . self::fmtLocs($best['locs']) . ' (P-31.1.4).';
        }
    }

    private function describirOrdenAlfa(array $best): string
    {
        if (count($best['grupos']) < 2) {
            $g = $best['grupos'] ? reset($best['grupos']) : null;
            if ($g && count($g['locs']) > 1) {
                return 'Sustituyentes idénticos se agrupan con prefijos multiplicadores (di, tri, tetra…; bis, tris… para prefijos complejos) (P-16.3).';
            }
            return '';
        }
        $n = [];
        foreach ($best['grupos'] as $nm => $g) { $n[] = $nm; }
        return 'Los prefijos se citan en orden alfanumérico, sin considerar multiplicadores (di, tri…) ni prefijos en cursiva (sec-, terc-): ' . implode(' < ', $n) . ' (P-14.5).';
    }

    /* ======================= salida SMILES ======================= */

    /** SMILES con mapas de átomo opcionales: $mapas[atomo] = número a mostrar. */
    public static function smiles(Molecula $m, array $mapas = [], int $inicio = 0): string
    {
        if ($m->n === 1) { return isset($mapas[0]) ? '[CH4:' . $mapas[0] . ']' : 'C'; }
        // 1ª pasada: DFS para detectar cierres de anillo
        $vis = []; $cierres = []; $padre = [$inicio => -1];
        $dfs = function (int $a) use (&$dfs, &$vis, &$cierres, &$padre, $m) {
            $vis[$a] = true;
            foreach ($m->adj[$a] as $b) {
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
        $digitos = []; $d = 1;
        foreach ($cierres as $k => [$x, $y]) { $digitos[$x][] = $d; $digitos[$y][] = $d; $d++; }
        $esCierre = [];
        foreach ($cierres as $k => $v) { $esCierre[$k] = true; }

        $emit = function (int $a, int $p) use (&$emit, $m, $mapas, $digitos, $esCierre): string {
            if ($m->el[$a] !== 'C') {
                $s = $m->el[$a];
            } elseif (isset($mapas[$a])) {
                $h = 4 - $m->grado($a);
                $s = '[C' . ($h > 0 ? 'H' . ($h > 1 ? $h : '') : '') . ':' . $mapas[$a] . ']';
            } else {
                $s = 'C';
            }
            foreach ($digitos[$a] ?? [] as $dg) { $s .= $dg > 9 ? '%' . $dg : (string)$dg; }
            $hijos = [];
            foreach ($m->adj[$a] as $b) {
                if ($b === $p) { continue; }
                $k = min($a, $b) . '-' . max($a, $b);
                if (isset($esCierre[$k])) { continue; }
                $hijos[] = $b;
            }
            $nh = count($hijos);
            foreach ($hijos as $i => $b) {
                $sub = $emit($b, $a);
                $s .= ($i < $nh - 1) ? '(' . $sub . ')' : $sub;
            }
            return $s;
        };
        return $emit($inicio, -1);
    }

    /** Fórmula molecular en orden de Hill: C, H y luego los demás elementos alfabéticamente. */
    public static function formula(Molecula $m): string
    {
        $cuenta = array_count_values($m->el);
        $nC = $cuenta['C'] ?? 0;
        $r = $m->numEnlaces() - $m->n + 1;
        $nX = $m->n - $nC;
        $h = 2 * $nC + 2 - 2 * $r - $nX;
        $f = 'C' . ($nC > 1 ? $nC : '') . ($h > 0 ? 'H' . ($h > 1 ? $h : '') : '');
        foreach (['Br', 'Cl', 'F', 'I'] as $x) {
            if (!empty($cuenta[$x])) { $f .= $x . ($cuenta[$x] > 1 ? $cuenta[$x] : ''); }
        }
        return $f;
    }
}

/* ======================= generador aleatorio ======================= */

final class GeneradorAlcanos
{
    /** Ramas disponibles por nivel (SMILES del sustituyente; el primer átomo es el de unión). */
    private const RAMAS = [
        1 => ['C', 'C', 'C', 'CC'],
        2 => ['C', 'C', 'CC', 'CC', 'CCC', 'C(C)C', 'CCCC', 'C(C)CC', 'CC(C)C', 'C(C)(C)C'],
        3 => [
            'C', 'CC', 'C(C)C', 'CCCC', 'C(C)CC', 'CC(C)C', 'C(C)(C)C',   // metil … terc-butil
            'CCCCC',                // pentil
            'CCC(C)C',              // 3-metilbutil (isopentil)
            'CC(C)(C)C',            // 2,2-dimetilpropil (neopentil)
            'C(C)(C)CC',            // 2-metilbutan-2-il (terc-pentil)
            'C(C)CCC',              // pentan-2-il
            'C(CC)CC',              // pentan-3-il
            'CC(C)CC',              // 2-metilbutil
            'C(C)C(C)C',            // 3-metilbutan-2-il
            'C(C)C(C)(C)C',         // 3,3-dimetilbutan-2-il
            'CC(C)C(C)C',           // 2,3-dimetilbutil
            'C(CC)C(C)C',           // 2-metilpentan-3-il
            'CCC(C)(C)C',           // 3,3-dimetilbutil
            'C(C)(C)C(C)C',         // 2,3-dimetilbutan-2-il
        ],
    ];

    /** Rango de la cadena principal y número de ramificaciones por nivel. */
    private const CADENA = [1 => [4, 7], 2 => [5, 10], 3 => [8, 20]];
    private const NRAMAS = [1 => [1, 2], 2 => [1, 3], 3 => [2, 5]];
    private const ANILLO = [1 => [3, 6], 2 => [3, 8], 3 => [5, 12]];
    private const MAX_C = 40;

    /** Halógenos por nivel: [mínimo, máximo] y reparto (Cl y Br más frecuentes). */
    private const NHAL = [1 => [1, 1], 2 => [1, 2], 3 => [1, 4]];
    private const HAL_POOL = ['Cl', 'Cl', 'Cl', 'Br', 'Br', 'Br', 'F', 'I'];

    public static function generar(string $tipo, int $nivel, bool $halogenos = false): Molecula
    {
        $m = self::esqueleto($tipo, $nivel, $halogenos);
        if ($halogenos) { self::halogenar($m, $nivel); }
        return $m;
    }

    /**
     * Coloca halógenos sobre carbonos con hidrógenos disponibles. En el nivel básico
     * no se colocan dos halógenos en el mismo carbono; en el avanzado se permiten
     * grupos CX2/CX3 (p. ej. triclorometil) y halógenos distintos.
     */
    private static function halogenar(Molecula $m, int $nivel): void
    {
        [$a, $b] = self::NHAL[$nivel];
        $k = random_int($a, $b);
        $carbonos = array_keys(array_filter($m->el, fn($e) => $e === 'C'));
        $mismo = $nivel < 3 || random_int(0, 1) === 1; // un solo tipo de halógeno
        $x0 = self::HAL_POOL[random_int(0, count(self::HAL_POOL) - 1)];
        for ($i = 0, $int = 0; $i < $k && $int < 50; $int++) {
            $c = $carbonos[random_int(0, count($carbonos) - 1)];
            if ($m->grado($c) >= 4) { continue; }
            $yaHal = count(array_filter($m->adj[$c], fn($v) => $m->el[$v] !== 'C'));
            if ($nivel === 1 && $yaHal > 0) { continue; }
            $x = $mismo ? $x0 : self::HAL_POOL[random_int(0, count(self::HAL_POOL) - 1)];
            $h = $m->nuevoAtomo($x);
            $m->enlazar($c, $h);
            $i++;
        }
    }

    private static function esqueleto(string $tipo, int $nivel, bool $halogenos): Molecula
    {
        $nivel = max(1, min(3, $nivel));
        if ($tipo === 'aleatorio') {
            $tipos = ['lineal', 'ramificado', 'ramificado', 'ciclico', 'ciclico_ramificado'];
            $tipo = $tipos[random_int(0, count($tipos) - 1)];
        }
        switch ($tipo) {
            case 'lineal':
                // con halógenos se incluyen metano y etano (diclorometano, 1,2-dibromoetano…)
                return self::cadena(random_int($halogenos ? 1 : 3, [1 => 8, 2 => 12, 3 => 20][$nivel]))[0];
            case 'ciclico':
                return self::anillo(random_int(3, [1 => 8, 2 => 10, 3 => 12][$nivel]), random_int(0, 1), $nivel);
            case 'ciclico_ramificado':
                [$a, $b] = self::ANILLO[$nivel];
                return self::anillo(random_int($a, $b), random_int(max(1, $nivel - 1), $nivel + 1), $nivel);
            case 'ramificado':
            default:
                return self::ramificado($nivel);
        }
    }

    /** @return array{0:Molecula,1:int[]} */
    private static function cadena(int $n): array
    {
        $m = new Molecula();
        $at = [];
        for ($i = 0; $i < $n; $i++) {
            $at[] = $m->nuevoAtomo();
            if ($i > 0) { $m->enlazar($at[$i - 1], $at[$i]); }
        }
        return [$m, $at];
    }

    private static function injertar(Molecula $m, int $ancla, string $rama): int
    {
        $sub = Molecula::desdeSmiles($rama);
        $map = [];
        foreach (range(0, $sub->n - 1) as $a) { $map[$a] = $m->nuevoAtomo(); }
        foreach ($sub->adj as $a => $vs) {
            foreach ($vs as $b) { if ($a < $b) { $m->enlazar($map[$a], $map[$b]); } }
        }
        $m->enlazar($ancla, $map[0]);
        return $sub->n;
    }

    private static function copiar(Molecula $m): Molecula
    {
        $c = new Molecula();
        $c->adj = $m->adj;
        $c->el = $m->el;
        $c->n = $m->n;
        return $c;
    }

    private static function nC(string $rama): int
    {
        return substr_count($rama, 'C');
    }

    private static function ramificado(int $nivel): Molecula
    {
        [$min, $max] = self::CADENA[$nivel];
        $L = random_int($min, $max);
        [$m, $at] = self::cadena($L);
        [$kmin, $kmax] = self::NRAMAS[$nivel];
        $k = random_int($kmin, min($kmax, max(1, intdiv($L - 2, 2) + 1)));
        $ramas = self::RAMAS[$nivel];
        $puestas = 0; $intentos = 0;
        while ($puestas < $k && $intentos < 60) {
            $intentos++;
            $pos = random_int(1, $L - 2);            // nunca en los extremos
            $a = $at[$pos];
            if ($m->grado($a) >= 4) { continue; }
            $rama = $ramas[random_int(0, count($ramas) - 1)];
            if ($m->n + self::nC($rama) > self::MAX_C) { continue; }
            $prueba = self::copiar($m);
            self::injertar($prueba, $a, $rama);
            if ($prueba->diametro() > 20) { continue; } // la cadena principal no debe exceder 20 C
            $m = $prueba;
            $puestas++;
        }
        return $m;
    }

    private static function anillo(int $N, int $k, int $nivel): Molecula
    {
        $m = new Molecula();
        $at = [];
        for ($i = 0; $i < $N; $i++) {
            $at[] = $m->nuevoAtomo();
            if ($i > 0) { $m->enlazar($at[$i - 1], $at[$i]); }
        }
        $m->enlazar($at[$N - 1], $at[0]);
        $ramas = self::RAMAS[$nivel];
        $k = min($k, $N);
        $puestas = 0; $intentos = 0;
        while ($puestas < $k && $intentos < 60) {
            $intentos++;
            $a = $at[random_int(0, $N - 1)];
            if ($m->grado($a) >= 4) { continue; }
            // geminales poco frecuentes en nivel 1
            if ($nivel === 1 && $m->grado($a) >= 3) { continue; }
            $rama = $ramas[random_int(0, count($ramas) - 1)];
            // sustituyente siempre más corto que el anillo: el anillo es el progenitor
            // tanto por IUPAC 2013 como por la regla tradicional (1979/1993)
            if (self::nC($rama) >= $N) { continue; }
            self::injertar($m, $a, $rama);
            $puestas++;
        }
        return $m;
    }
}
