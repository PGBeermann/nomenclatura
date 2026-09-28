<?php
declare(strict_types=1);
session_set_cookie_params([
    'lifetime' => 0, 'path' => '/',
    'secure'   => !empty($_SERVER['HTTPS']) && $_SERVER['HTTPS'] !== 'off',
    'httponly' => true, 'samesite' => 'Strict',
]);
session_start();
if (empty($_SESSION['csrf'])) { $_SESSION['csrf'] = bin2hex(random_bytes(32)); }
header('X-Content-Type-Options: nosniff');
header('Referrer-Policy: same-origin');

/** Casillas de familias: valor => [etiqueta, ejemplo] */
$familias = [
    'halogenos'  => ['Halógenos (F, Cl, Br, I)', '2-cloropropano'],
    'alquenos'   => ['Alquenos (C=C)', 'but-2-eno'],
    'alquinos'   => ['Alquinos (C≡C)', 'prop-1-ino'],
    'aromaticos' => ['Aromáticos (benceno)', 'etilbenceno'],
    'biciclos'   => ['Dos anillos (fusionados, puente, espiro, bifenilo…)', 'naftaleno, biciclo[2.2.1]heptano'],
    'heterociclos' => ['Heterociclos (piridina, furano, pirrolidina, indol…)', 'piridin-3-ol, oxolano, 1H-indol'],
    'alcoholes'  => ['Alcoholes', 'propan-2-ol'],
    'fenoles'    => ['Fenoles', '4-metilfenol'],
    'eteres'     => ['Éteres', 'metoxietano'],
    'aminas'     => ['Aminas', 'N-metiletanamina'],
    'aldehidos'  => ['Aldehídos', 'butanal'],
    'cetonas'    => ['Cetonas', 'pentan-3-ona'],
    'acidos'     => ['Ácidos carboxílicos', 'ácido propanoico'],
    'esteres'    => ['Ésteres', 'acetato de etilo'],
    'amidas'     => ['Amidas', 'N-metilacetamida'],
    'nitrilos'   => ['Nitrilos', 'butanonitrilo'],
    'nitro'      => ['Nitro', 'nitrobenceno'],
    'azufre'     => ['Azufre (tioles, sulfuros, sulfóxidos, sulfonas, ácidos sulfónicos)', 'etanotiol, (metilsulfanil)benceno'],
];
?><!doctype html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="csrf-token" content="<?= htmlspecialchars($_SESSION['csrf'], ENT_QUOTES) ?>">
<title>Nomenclatura IUPAC de compuestos orgánicos</title>
<link rel="stylesheet" href="assets/estilos.css?v=4.2">
</head>
<body>
<header class="cabecera">
  <div class="cabecera__inner cabecera__marca">
    <a href="https://www.unachi.ac.pa/" target="_blank" rel="noopener" title="Universidad Autónoma de Chiriquí">
      <img class="cabecera__logo" src="assets/img/unachi-logo.png" alt="Logo de la Universidad Autónoma de Chiriquí (UNACHI)" width="84" height="84" decoding="async">
    </a>
    <div class="cabecera__texto">
      <p class="cabecera__inst">Universidad Autónoma de Chiriquí · Facultad de Ciencias Naturales y Exactas · Escuela de Química</p>
      <h1>Nomenclatura IUPAC de compuestos orgánicos</h1>
      <p class="cabecera__sub">Hidrocarburos, compuestos aromáticos y grupos funcionales. Genere o dibuje una estructura, escriba su nombre y compruébelo con las Recomendaciones IUPAC 2013.</p>
    </div>
  </div>
</header>

<main class="contenedor">
  <section class="panel controles" aria-label="Opciones del ejercicio">
    <label>Esqueleto
      <select id="tipo">
        <option value="aleatorio">Aleatorio</option>
        <option value="lineal">Cadena lineal</option>
        <option value="ramificado" selected>Cadena ramificada</option>
        <option value="ciclico">Anillo simple</option>
        <option value="ciclico_ramificado">Anillo con ramificaciones</option>
      </select>
    </label>
    <label>Nivel
      <select id="nivel">
        <option value="1">Básico (1 grupo)</option>
        <option value="2" selected>Intermedio (1–2 grupos)</option>
        <option value="3">Avanzado (2–3 grupos)</option>
      </select>
    </label>
    <div class="botones">
      <button id="btnGenerar" class="btn btn--primario" type="button">Generar</button>
      <button id="btnDesplegar" class="btn btn--secundario" type="button" disabled>Desplegar</button>
      <button id="btnDibujar" class="btn btn--linea" type="button" title="Dibuje su propia estructura en el editor y verifique su nombre">Dibujar</button>
    </div>
    <p class="marcador" aria-live="polite">Ejercicios: <b id="cntEj">0</b> · Aciertos: <b id="cntOk">0</b></p>

    <fieldset class="familias" id="familias">
      <legend>Tipos de compuestos que se incluirán en los ejercicios
        <span class="familias__acciones">
          <button type="button" class="enlace" id="btnTodas">marcar todas</button> ·
          <button type="button" class="enlace" id="btnNinguna">solo alcanos</button>
        </span>
      </legend>
      <?php foreach ($familias as $v => [$etq, $ej]): ?>
        <label class="check check--familia" title="Ejemplo: <?= htmlspecialchars($ej, ENT_QUOTES) ?>">
          <input type="checkbox" name="familia" value="<?= $v ?>"> <?= htmlspecialchars($etq) ?>
        </label>
      <?php endforeach; ?>
      <p class="familias__nota">Sin casillas marcadas se generan alcanos y cicloalcanos. Con varias marcadas, en los niveles intermedio y avanzado se combinan grupos para practicar el orden de prioridad (el grupo principal va como sufijo y los demás como prefijos).</p>
    </fieldset>
  </section>

  <div class="rejilla">
    <section class="panel" aria-label="Estructura">
      <div class="panel__titulo">
        <h2>Estructura</h2>
        <span id="formula" class="formula" title="Fórmula molecular"></span>
      </div>
      <p id="avisoDibujo" class="aviso" hidden>Modo dibujo: dibuje una molécula con C, H, N, O, S y halógenos (cadenas; carbociclos y heterociclos simples, fusionados, con puente o espiro de hasta dos anillos; enlaces sencillos, dobles o triples), escriba su nombre y pulse <b>Verificar</b> o <b>Desplegar</b>.</p>
      <div id="jsme_container" class="jsme"><p class="cargando">Cargando editor JSME…</p></div>
      <div class="opciones">
        <label class="check"><input type="checkbox" id="chkNumeracion" disabled> Mostrar numeración de la estructura principal</label>
        <label class="check"><input type="checkbox" id="chkEditar"> Permitir edición (modo editor JSME)</label>
        <span id="origen" class="origen"></span>
      </div>
    </section>

    <section class="panel" aria-label="Respuesta">
      <h2>Su respuesta</h2>
      <form id="frmRespuesta" autocomplete="off">
        <input id="respuesta" type="text" placeholder="p. ej. ácido 3-hidroxibutanoico" disabled spellcheck="false">
        <button class="btn btn--linea" type="submit" id="btnVerificar" disabled>Verificar</button>
      </form>
      <p id="retro" class="retro" aria-live="polite"></p>

      <div id="solucion" class="solucion" hidden>
        <h3>Nombre IUPAC preferido (2013)</h3>
        <p class="nombre" id="nombre"></p>
        <div id="bloqueAlt" class="alternativas" hidden>
          <h3>Otras formas aceptadas</h3>
          <div id="listaAlt"></div>
          <div class="tabla-eq" id="bloqueEq" hidden>
            <table>
              <caption>Equivalencia de los sustituyentes</caption>
              <thead><tr><th>Preferido 2013</th><th>Retenido / tradicional</th><th>Sistemático 1979</th></tr></thead>
              <tbody id="tablaEq"></tbody>
            </table>
          </div>
        </div>
        <dl class="alternativos">
          <div><dt>Inglés (PIN)</dt><dd id="nombreEn"></dd></div>
          <div id="filaEnAlt" hidden><dt>Inglés (alternativos)</dt><dd id="nombreEnAlt"></dd></div>
        </dl>
        <h3>Justificación</h3>
        <ol id="explicacion" class="explicacion"></ol>
      </div>
      <p id="error" class="error" role="alert" hidden></p>
    </section>
  </div>

  <details class="panel reglas">
    <summary>Resumen de reglas aplicadas (IUPAC 2013)</summary>
    <ol>
      <li><b>Grupo principal</b>: entre los grupos presentes se elige el de mayor prioridad (P-41): ácido carboxílico &gt; éster &gt; amida &gt; nitrilo &gt; aldehído &gt; cetona &gt; alcohol y fenol &gt; amina. Se expresa como <b>sufijo</b> (-oico, -oato, -amida, -nitrilo, -al, -ona, -ol, -amina); los demás grupos se citan como <b>prefijos</b> (carboxi, alcoxicarbonil, carbamoil, ciano, oxo/formil, hidroxi, amino). Éteres (alcoxi), halógenos y nitro son siempre prefijos.</li>
      <li><b>Estructura principal</b>: la que tiene el máximo número de grupos principales (P-44.1.1); a igualdad, un anillo es preferido a una cadena (P-44.1.2.2).</li>
      <li><b>Cadena principal</b> (2013): la más larga (P-44.3); luego la de más enlaces múltiples y, después, más dobles enlaces (P-44.4.1). <i>Cambio respecto de 1979/1993</i>, que anteponían la insaturación a la longitud: por eso ahora aparecen prefijos como «metiliden» (3-metilidenhexano, antes 2-etilpent-1-eno).</li>
      <li><b>Numeración</b> (P-31.1.4): localizadores más bajos, por este orden, para el grupo principal, los enlaces múltiples en conjunto («-eno/-ino»), los dobles enlaces, todos los prefijos juntos y, finalmente, el prefijo citado primero en orden alfanumérico.</li>
      <li><b>Sufijos en anillos</b>: si el carbono del grupo no pertenece al anillo se usan -carboxílico, -carboxilato, -carboxamida, -carbonitrilo y -carbaldehído (ácido ciclohexanocarboxílico).</li>
      <li><b>Nombres retenidos preferidos</b>: fenol, anilina, ácido benzoico, benzaldehído, benzamida, benzonitrilo (admiten sustitución en el anillo); ácido fórmico, ácido acético, ácido oxálico, formaldehído, acetaldehído, acetamida, acetonitrilo; tolueno, xileno y anisol solo sin sustituir; acetileno. La acetona se llama propan-2-ona.</li>
      <li><b>Orden alfanumérico</b> de los prefijos sin considerar di-, tri-, sec-, terc-; sí se considera iso- y los multiplicadores dentro de un prefijo compuesto (P-14.5). <b>Multiplicadores</b>: di, tri… para prefijos simples; bis, tris… para prefijos compuestos (P-16.3).</li>
      <li><b>Prefijos sustituyentes</b> preferidos (P-29): propan-2-il, butan-2-il, etenil, prop-2-en-1-il, metiliden; se retienen terc-butil, fenil y bencil. Las formas isopropil, sec-butil, vinil, alil e isopropiliden se aceptan en nomenclatura general.</li>
      <li><b>Ésteres, aminas y amidas</b>: «propanoato de etilo»; los sustituyentes del nitrógeno llevan el localizador N (N,N-dimetiletanamina, N-fenilacetamida). Se muestran también los nombres de clase funcional aceptados (etil metil éter, alcohol isopropílico, cloruro de vinilo, etil(metil)amina).</li>
      <li><b>Dos anillos</b> (fase 2): biciclos fusionados con nombre retenido (naftaleno, indeno, azuleno, pentaleno, heptaleno) y su numeración fija (4a, 8a); hidrógeno indicado (1H-indeno), prefijos hidro (2,3-dihidro-1H-indeno, decahidronaftaleno) e hidrógeno añadido (3,4-dihidronaftalen-1(2H)-ona) (P-25, P-31.1.4.2.4, P-14.7); biciclos con puente de von Baeyer (biciclo[2.2.1]heptano, P-23.2); espiro (espiro[4.5]decano, P-24.2); ensamblajes (1,1'-bifenilo, P-28); nomenclatura multiplicativa (1,1'-metilendibenceno, 4,4'-(propano-2,2-diil)difenol, P-15.3); antigüedad de anillos: más anillos, más átomos, menos hidrogenado (P-44.2, P-44.4).</li>
      <li><b>Heterociclos</b> (fase 3A): nombres retenidos (pirrol, furano, tiofeno, imidazol, pirazol, piridina, pirimidina, pirazina, pirano; pirrolidina, piperidina, piperazina, morfolina; indol, quinolina, isoquinolina, 1-benzofurano, 1H-bencimidazol, purina…) y de Hantzsch-Widman (oxirano, azetidina, oxolano, 1,3-oxazol, 1,4-dioxano, 1,3,5-triazina) con los localizadores más bajos para los heteroátomos, primero O, luego S, luego N (P-22.2, P-31.1.4.2.2); prefijos de reemplazo en puentes y espiro (1-azabiciclo[2.2.2]octano, 1,4-dioxaespiro[4.5]decano, P-15.4); lactonas y lactamas como cetonas del heterociclo (oxolan-2-ona, pirrolidin-2-ona); hidrógeno añadido (piridin-2(1H)-ona, pirimidina-2,4(1H,3H)-diona). Los heterociclos con N son preferidos a los demás anillos (P-44.2).</li>
      <li><b>Azufre</b>: tioles (etanotiol, sufijo -tiol, prefijo sulfanil), sulfuros ((metilsulfanil)metano), sulfóxidos y sulfonas ((metanosulfinil)metano, (metanosulfonil)metano) y ácidos sulfónicos (ácido bencenosulfónico, prefijo sulfo) (P-63, P-65.3).</li>
      <li><b>Omisión de localizadores</b> (P-14.3.4): etanol, cloroetano, ciclohexanol, propeno, ciclohexeno; pero 3-metilciclohex-1-eno, propan-2-ona, butan-2-ona.</li>
    </ol>
    <p><b>Dr. Pedro González Beermann · UNACHI 2026</b></p>
    <p class="ref">Ref.: Favre, H. A.; Powell, W. H. <i>Nomenclature of Organic Chemistry. IUPAC Recommendations and Preferred Names 2013</i>. RSC, 2014. doi:10.1039/9781849733069 · Editor molecular: Bienfait, B.; Ertl, P. <i>J. Cheminform.</i> 2013, 5, 24. doi:10.1186/1758-2946-5-24</p>
  </details>
</main>

<script src="assets/app.js?v=4.2"></script>
<script src="jsme/jsme.nocache.js"></script>
</body>
</html>
