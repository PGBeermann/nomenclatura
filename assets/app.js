/* app.js – Interfaz del generador de alcanos y haloalcanos (JSME + api.php)  v3.0 */
'use strict';

const estado = {
  jsme: null,
  id: null,            // ejercicio activo en el servidor
  smiles: '',          // SMILES de la estructura (sin numeración)
  smilesNum: '',       // SMILES con los localizadores como mapas de átomo
  smilesCargado: '',   // lo que JSME devolvía justo después de la última carga programática
  longitud: 8,         // carbonos de la cadena/anillo principal (para escalar el dibujo)
  propio: false,       // true si la estructura fue dibujada por el estudiante
  desplegado: false, acertado: false, ej: 0, ok: 0
};
const $ = (id) => document.getElementById(id);
const CSRF = document.querySelector('meta[name="csrf-token"]').content;

/* ---------- JSME: jsme.nocache.js invoca esta función al cargar ---------- */
window.jsmeOnLoad = function () {
  const cont = $('jsme_container');
  cont.innerHTML = '';
  const w = Math.max(280, cont.clientWidth - 2);
  estado.jsme = new JSApplet.JSME('jsme_container', w + 'px', '400px', {
    options: 'depict,nosearchinchiKey,noquery'
  });
  window.addEventListener('resize', ajustarTamano);
  generar();
};

function ajustarTamano() {
  if (!estado.jsme) return;
  const w = Math.max(280, $('jsme_container').clientWidth - 2);
  estado.jsme.setSize(w + 'px', '400px');
}

function smilesActual() {
  try { return (estado.jsme && estado.jsme.smiles()) || ''; } catch (e) { return ''; }
}

function mostrarMolecula(smiles) {
  if (!estado.jsme) return;
  estado.jsme.readGenericMolecularInput(smiles);
  // escala (debe aplicarse DESPUÉS de leer): ~21,5 px por enlace a escala 1
  const ancho = $('jsme_container').clientWidth;
  const L = estado.longitud || 8;
  const escala = Math.max(0.75, Math.min(2.4, (0.82 * ancho) / (21.5 * (L + 1))));
  const a = estado.jsme;
  if (typeof a.setMolecularAreaScale === 'function') a.setMolecularAreaScale(escala);
  if (typeof a.setMolecularAreaLineWidth === 'function') a.setMolecularAreaLineWidth(1);
  if (typeof a.setAtomMolecularAreaFontSize === 'function') a.setAtomMolecularAreaFontSize(11);
  estado.smilesCargado = smilesActual();
}

/** true si el estudiante cambió la estructura en el editor desde la última carga. */
function estructuraModificada() {
  return estado.jsme !== null && smilesActual() !== estado.smilesCargado;
}

function modoEdicion(activo) {
  $('chkEditar').checked = activo;
  if (estado.jsme) estado.jsme.options(activo ? 'nodepict' : 'depict');
}

/* ---------- comunicación con el servidor ---------- */
async function api(accion, datos = {}) {
  const r = await fetch('api.php', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json', 'X-CSRF-Token': CSRF },
    credentials: 'same-origin',
    body: JSON.stringify({ accion, ...datos })
  });
  const j = await r.json().catch(() => ({ error: 'Respuesta no válida del servidor' }));
  if (!r.ok || j.error) throw new Error(j.error || ('HTTP ' + r.status));
  return j;
}

function mostrarError(msg) {
  const e = $('error');
  e.textContent = msg;
  e.hidden = !msg;
}

function limpiarSolucion() {
  $('solucion').hidden = true;
  $('retro').textContent = '';
  $('retro').className = 'retro';
  $('chkNumeracion').checked = false;
  $('chkNumeracion').disabled = true;
  estado.desplegado = false;
  estado.acertado = false;
}

function ponerFormula(f) {
  $('formula').innerHTML = f ? f.replace(/(\d+)/g, '<sub>$1</sub>') : '';
}

function habilitarRespuesta() {
  $('respuesta').disabled = false;
  $('btnVerificar').disabled = false;
  $('btnDesplegar').disabled = false;
}

/* ---------- <Generar> ---------- */
async function generar() {
  mostrarError('');
  $('btnGenerar').disabled = true;
  try {
    const d = await api('generar', { tipo: $('tipo').value, nivel: +$('nivel').value, halogenos: $('chkHalogenos').checked });
    Object.assign(estado, { id: d.id, smiles: d.smiles, smilesNum: '', longitud: d.longitud, propio: false });
    $('avisoDibujo').hidden = true;
    mostrarMolecula(d.smiles);
    ponerFormula(d.formula);
    $('origen').textContent = 'Estructura generada';
    limpiarSolucion();
    $('respuesta').value = '';
    habilitarRespuesta();
    $('cntEj').textContent = ++estado.ej;
    $('respuesta').focus({ preventScroll: true });
  } catch (err) {
    mostrarError(err.message);
  } finally {
    $('btnGenerar').disabled = false;
  }
}

/* ---------- <Dibujar>: lienzo vacío en modo editor ---------- */
function dibujar() {
  if (!estado.jsme) return;
  mostrarError('');
  modoEdicion(true);
  estado.jsme.reset();
  Object.assign(estado, { id: null, smiles: '', smilesNum: '', smilesCargado: '', propio: true });
  ponerFormula('');
  $('origen').textContent = 'Estructura dibujada por el estudiante';
  $('avisoDibujo').hidden = false;
  limpiarSolucion();
  $('respuesta').value = '';
  habilitarRespuesta();
}

/**
 * Si la estructura del editor cambió (o se dibujó desde cero), la envía al
 * servidor para crear un nuevo ejercicio con ella. Devuelve true si hay ejercicio.
 */
async function sincronizarEstructura() {
  if (estado.id && !estructuraModificada()) return true;
  const smi = smilesActual();
  if (!smi) { mostrarError('Dibuje una estructura en el editor primero.'); return false; }
  let d;
  try {
    d = await api('analizar', { smiles: smi });
  } catch (err) {
    estado.id = null;                     // la estructura del editor no es válida
    limpiarSolucion();
    throw err;
  }
  const eraGenerada = !estado.propio;
  Object.assign(estado, { id: d.id, smiles: smi, smilesNum: '', longitud: d.longitud, propio: true });
  estado.smilesCargado = smi;             // el dibujo del estudiante se conserva tal cual
  ponerFormula(d.formula);
  $('origen').textContent = eraGenerada ? 'Estructura modificada por el estudiante' : 'Estructura dibujada por el estudiante';
  limpiarSolucion();
  return true;
}

/* ---------- <Desplegar> ---------- */
async function desplegar() {
  mostrarError('');
  try {
    if (!(await sincronizarEstructura())) return;
    const d = await api('desplegar', { id: estado.id });
    estado.desplegado = true;
    estado.smilesNum = d.smiles_num;

    $('nombre').textContent = d.nombre;
    $('filaTrad').hidden = !d.nombre_trad;
    $('nombreTrad').textContent = d.nombre_trad || '';
    $('filaSis').hidden = !d.nombre_sis;
    $('nombreSis').textContent = d.nombre_sis || '';

    const radico = [d.radicofuncional, d.radicofuncional_pin].filter(Boolean);
    $('filaRadico').hidden = radico.length === 0;
    $('nombreRadico').textContent = radico.join('  ·  ');
    const tb = $('tablaEq');
    tb.innerHTML = '';
    (d.equivalencias || []).forEach((e) => {
      const tr = document.createElement('tr');
      [e.pin, e.trad, e.sis].forEach((t) => { const td = document.createElement('td'); td.textContent = t; tr.appendChild(td); });
      tb.appendChild(tr);
    });
    $('bloqueEq').hidden = !(d.equivalencias && d.equivalencias.length);
    $('bloqueAlt').hidden = !d.nombre_trad && !d.nombre_sis && $('bloqueEq').hidden && radico.length === 0;

    $('nombreEn').textContent = d.nombre_en;
    const enAlt = [d.nombre_en_trad, d.nombre_en_sis, d.radicofuncional_en].filter(Boolean);
    $('filaEnAlt').hidden = enAlt.length === 0;
    $('nombreEnAlt').textContent = enAlt.join('  ·  ');

    const ol = $('explicacion');
    ol.innerHTML = '';
    d.explicacion.forEach((t) => { const li = document.createElement('li'); li.textContent = t; ol.appendChild(li); });
    $('solucion').hidden = false;
    $('chkNumeracion').disabled = false;
    $('chkNumeracion').checked = true;
    mostrarMolecula(estado.smilesNum);
  } catch (err) {
    mostrarError(err.message);
  }
}

/* ---------- Verificar respuesta del estudiante ---------- */
async function verificar(ev) {
  ev.preventDefault();
  mostrarError('');
  try {
    if (!(await sincronizarEstructura())) return;
    const d = await api('verificar', { id: estado.id, respuesta: $('respuesta').value });
    const r = $('retro');
    r.textContent = d.mensaje;
    r.className = 'retro ' + (d.correcto ? 'ok' : (d.nivel === 'casi' ? 'casi' : 'mal'));
    if (d.correcto && !estado.acertado && !estado.desplegado && !estado.propio) {
      estado.acertado = true;
      $('cntOk').textContent = ++estado.ok;
    }
  } catch (err) {
    mostrarError(err.message);
  }
}

/* ---------- eventos ---------- */
$('btnGenerar').addEventListener('click', generar);
$('btnDesplegar').addEventListener('click', desplegar);
$('btnDibujar').addEventListener('click', dibujar);
$('frmRespuesta').addEventListener('submit', verificar);
$('chkNumeracion').addEventListener('change', (e) => {
  if (estructuraModificada()) { e.target.checked = false; return; }
  mostrarMolecula(e.target.checked && estado.smilesNum ? estado.smilesNum : estado.smiles);
});
$('chkEditar').addEventListener('change', (e) => modoEdicion(e.target.checked));
