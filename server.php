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
?><!doctype html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="csrf-token" content="<?= htmlspecialchars($_SESSION['csrf'], ENT_QUOTES) ?>">
<title>Nomenclatura IUPAC de alcanos</title>
<link rel="stylesheet" href="assets/estilos.css?v=3.0">
<style>
  /* Logo institucional en el encabezado (esquina superior izquierda) */
  .cabecera__marca { display: flex; align-items: center; gap: 1.1rem; }
  .cabecera__logo  { flex: 0 0 auto; height: 84px; width: auto; display: block; }
  .cabecera__texto { min-width: 0; }
  @media (max-width: 640px) {
    .cabecera__marca { gap: .75rem; align-items: flex-start; }
    .cabecera__logo  { height: 56px; }
  }
</style>
</head>
<body>
<header class="cabecera">
  <div class="cabecera__inner cabecera__marca">
    <a href="https://www.unachi.ac.pa/" target="_blank" rel="noopener" title="Universidad Autónoma de Chiriquí">
      <img class="cabecera__logo" src="assets/img/logo-unachi.png" alt="Logo de la Universidad Autónoma de Chiriquí (UNACHI)" width="84" height="84" decoding="async">
    </a>
    <div class="cabecera__texto">
      <p class="cabecera__inst">Universidad Autónoma de Chiriquí · Facultad de Ciencias Naturales y Exactas · Escuela de Química</p>
      <h1>Nomenclatura IUPAC de alcanos, cicloalcanos y haloalcanos</h1>
      <p class="cabecera__sub">Genere o dibuje una estructura, escriba su nombre y compruébelo con las Recomendaciones IUPAC 2013.</p>
    </div>
  </div>
</header>

<main class="contenedor">
  <section class="panel controles" aria-label="Opciones del ejercicio">
    <label>Tipo de estructura
      <select id="tipo">
        <option value="aleatorio">Aleatorio</option>
        <option value="lineal">Alcano lineal (no ramificado)</option>
        <option value="ramificado" selected>Alcano ramificado</option>
        <option value="ciclico">Cicloalcano simple</option>
        <option value="ciclico_ramificado">Cicloalcano con ramificaciones</option>
      </select>
    </label>
    <label>Nivel
      <select id="nivel">
        <option value="1">Básico</option>
        <option value="2" selected>Intermedio</option>
        <option value="3">Avanzado</option>
      </select>
    </label>
    <label class="check check--control"><input type="checkbox" id="chkHalogenos"> Incluir halógenos (F, Cl, Br, I)</label>
    <div class="botones">
      <button id="btnGenerar" class="btn btn--primario" type="button">Generar</button>
      <button id="btnDesplegar" class="btn btn--secundario" type="button" disabled>Desplegar</button>
      <button id="btnDibujar" class="btn btn--linea" type="button" title="Dibuje su propio alcano en el editor y verifique su nombre">Dibujar</button>
    </div>
    <p class="marcador" aria-live="polite">Ejercicios: <b id="cntEj">0</b> · Aciertos: <b id="cntOk">0</b></p>
  </section>

  <div class="rejilla">
    <section class="panel" aria-label="Estructura">
      <div class="panel__titulo">
        <h2>Estructura</h2>
        <span id="formula" class="formula" title="Fórmula molecular"></span>
      </div>
      <p id="avisoDibujo" class="aviso" hidden>Modo dibujo: dibuje un alcano, cicloalcano o haloalcano (carbonos, F, Cl, Br, I y solo enlaces sencillos), escriba su nombre y pulse <b>Verificar</b> o <b>Desplegar</b>.</p>
      <div id="jsme_container" class="jsme"><p class="cargando">Cargando editor JSME…</p></div>
      <div class="opciones">
        <label class="check"><input type="checkbox" id="chkNumeracion" disabled> Mostrar numeración de la cadena/anillo principal</label>
        <label class="check"><input type="checkbox" id="chkEditar"> Permitir edición (modo editor JSME)</label>
        <span id="origen" class="origen"></span>
      </div>
    </section>

    <section class="panel" aria-label="Respuesta">
      <h2>Su respuesta</h2>
      <form id="frmRespuesta" autocomplete="off">
        <input id="respuesta" type="text" placeholder="p. ej. 4-etil-2,2-dimetilheptano" disabled spellcheck="false">
        <button class="btn btn--linea" type="submit" id="btnVerificar" disabled>Verificar</button>
      </form>
      <p id="retro" class="retro" aria-live="polite"></p>

      <div id="solucion" class="solucion" hidden>
        <h3>Nombre IUPAC preferido (2013)</h3>
        <p class="nombre" id="nombre"></p>
        <div id="bloqueAlt" class="alternativas" hidden>
          <h3>Otras formas aceptadas</h3>
          <div class="alt" id="filaTrad" hidden>
            <span class="alt__etq">Con prefijos tradicionales retenidos (IUPAC 1979/1993)</span>
            <span class="alt__nom" id="nombreTrad"></span>
          </div>
          <div class="alt" id="filaSis" hidden>
            <span class="alt__etq">Con prefijos sistemáticos numerados desde el punto de unión (IUPAC 1979)</span>
            <span class="alt__nom" id="nombreSis"></span>
          </div>
          <div class="alt" id="filaRadico" hidden>
            <span class="alt__etq">Nombre radicofuncional o de clase funcional (aceptado en nomenclatura general, no preferido; P-61.3)</span>
            <span class="alt__nom" id="nombreRadico"></span>
          </div>
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
      <li><b>Cadena principal</b>: la cadena continua más larga (P-44.3). En cicloalcanos, el anillo es el hidruro progenitor (P-44.1.2.2).</li>
      <li>Si hay varias cadenas de igual longitud: la de <b>más sustituyentes</b>; luego la de <b>localizadores más bajos</b>; luego la que da el localizador más bajo al prefijo citado primero en orden alfanumérico (P-45.2).</li>
      <li><b>Numeración</b>: conjunto de localizadores más bajo, comparado término a término en el primer punto de diferencia (P-14.4, P-31.1.4).</li>
      <li><b>Orden alfanumérico</b> de los prefijos sin considerar di-, tri-, sec-, terc-; sí se considera iso- y los multiplicadores dentro de un prefijo complejo (P-14.5).</li>
      <li><b>Multiplicadores</b>: di, tri, tetra… para prefijos simples; bis, tris, tetrakis… para prefijos complejos, entre paréntesis (P-16.3).</li>
      <li><b>Prefijos sustituyentes</b> preferidos (P-29): la cadena más larga que contiene el átomo de unión, con la valencia libre en el localizador más bajo: propan-2-il, butan-2-il, 2-metilpropil, 3-metilbutan-2-il… Se retiene terc-butil (P-29.6).</li>
      <li><b>Halógenos</b>: se nombran solo como prefijos (fluoro, cloro, bromo, yodo); cuentan como sustituyentes para elegir la cadena y numerar, y se alfabetizan con los grupos alquilo (P-61.3.1, P-14.5). Ej.: 1-bromo-4-cloro-2-metilpentano. Para los monohaloalcanos se da también el nombre radicofuncional (bromuro de isopropilo).</li>
      <li><b>Prefijos sec- y neo-</b>: solo existen como prefijos retenidos no sustituidos en <i>sec</i>-butil y neopentil (IUPAC 1979, A-2.25). No se usan «sec-pentil» ni «neohexil» porque son ambiguos. <i>sec</i>- (en cursiva) se ignora al alfabetizar; neo- e iso- sí cuentan.</li>
      <li><b>Formas alternativas</b> aceptadas en nomenclatura general: prefijos retenidos no sustituidos (isopropil, sec-butil, isobutil, terc-butil, isopentil, neopentil, terc-pentil) y prefijos numerados desde el átomo de unión (1-metiletil, 1,2-dimetilpropil), según las reglas IUPAC 1979 (A-2.25, A-2.6).</li>
      <li>En un cicloalcano monosustituido no se escribe el localizador 1 (P-14.3.4).</li>
    </ol>
    <p><b>Dr. Pedro González Beermann Unachi-2026</b></p>
    <p class="ref">Ref.: Favre, H. A.; Powell, W. H. <i>Nomenclature of Organic Chemistry. IUPAC Recommendations and Preferred Names 2013</i>. RSC, 2014. doi:10.1039/9781849733069 · Editor molecular: Bienfait, B.; Ertl, P. <i>J. Cheminform.</i> 2013, 5, 24. doi:10.1186/1758-2946-5-24</p>
  </details>
</main>

<script src="assets/app.js?v=3.0"></script>
<script src="jsme/jsme.nocache.js"></script>
</body>
</html>