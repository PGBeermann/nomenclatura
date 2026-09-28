# Fase 3A: heterociclos (retenidos, Hantzsch-Widman, fusionados, reemplazo) y compuestos de azufre:
# round-trip con OPSIN e invariancia ante el orden de los átomos.
import random, subprocess, sys
from rdkit import Chem
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
from validar_org import validar, DESPLAZ
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
BASES = ['c1ccncc1', 'c1cc[nH]c1', 'c1ccoc1', 'c1ccsc1', 'c1cnc[nH]1', 'c1cn[nH]c1', 'c1cncnc1', 'c1cnccn1', 'c1ccnnc1',
         'C1CCNC1', 'C1CCNCC1', 'C1COCCN1', 'C1CCOC1', 'C1CCOCC1', 'C1CO1', 'C1CN1', 'C1CS1', 'C1COC1', 'C1CNC1', 'C1COCO1',
         'C1COCCO1', 'c1cocn1', 'c1cscn1', 'c1conc1', 'c1ncncn1', 'c1nc[nH]n1', 'C1CCCNCC1', 'C1CCSC1', 'C1CNCCN1', 'C1=COCC1',
         'C1=CCCOC1', 'C1=CCNC1', 'c1ccc2[nH]ccc2c1', 'c1ccc2ncccc2c1', 'c1ccc2cnccc2c1', 'c1ccc2occc2c1', 'c1ccc2sccc2c1',
         'c1ccc2[nH]cnc2c1', 'c1ccc2[nH]ncc2c1', 'c1ccc2ocnc2c1', 'c1ccc2scnc2c1', 'c1ccc2ncncc2c1', 'c1ccc2nccnc2c1',
         'C1Cc2ccccc2N1', 'C1CCc2ccccc2N1', 'C1=Cc2ccccc2OC1', 'C1CCc2ccccc2O1', 'c1ncc2nc[nH]c2n1', 'C1CN2CCC1CC2',
         'C1CC2CCC1N2', 'C1CCC2(CC1)OCCO2', 'C1CC2CCC1O2', 'C1CC2(C1)CNC2', 'c1ccc(cc1)-c1ccccn1', 'c1ccnc(c1)-c1ccccn1',
         'C1CCC(CC1)N1CCCC1', 'c1ccc(Oc2ccncc2)cc1', 'CSc1ccccc1', 'CCSCC', 'CS(=O)C', 'CCS(=O)(=O)CC', 'OS(=O)(=O)c1ccccc1',
         'SCCCC', 'CC(C)S', 'c1ccc(Sc2ccccc2)cc1']
GR = ['O', 'N', 'NC', '=O', 'C=O', 'C(=O)O', 'C(=O)OC', 'C(N)=O', 'C#N', 'OC', 'Cl', 'Br', 'F', 'C', 'C', 'CC', 'C(C)C', 'S', 'SC', 'S(C)=O', 'S(C)(=O)=O', 'S(=O)(=O)O', 'c1ccccc1']
def aleatoria():
    s = random.choice(BASES)
    for _ in range(random.choice([0, 1, 1, 2, 2, 3])):
        mm = Chem.MolFromSmiles(s)
        cand = [a.GetIdx() for a in mm.GetAtoms() if a.GetTotalNumHs() >= 1 and (a.GetAtomicNum() == 6 or (a.GetAtomicNum() == 7 and a.IsInRing()))]
        if not cand: break
        g = random.choice(GR); a = random.choice(cand)
        if g == '=O' and (mm.GetAtomWithIdx(a).GetTotalNumHs() < 2): continue
        rw = Chem.RWMol(mm); frag = Chem.MolFromSmiles(g[1:] if g.startswith('=') else g)
        off = rw.GetNumAtoms(); combo = Chem.RWMol(Chem.CombineMols(rw, frag))
        combo.AddBond(a, off, Chem.BondType.DOUBLE if g.startswith('=') else Chem.BondType.SINGLE)
        try:
            mol = combo.GetMol(); Chem.SanitizeMol(mol); s = Chem.MolToSmiles(mol)
        except Exception: pass
    return s
N = int(sys.argv[1]); smis = list({aleatoria() for _ in range(N)})
res = subprocess.run(['php', 'nombrar_org.php'], input='\n'.join(smis) + '\n', capture_output=True, text=True).stdout.rstrip('\n').split('\n')
filas = []; errs = {}
for s, r in zip(smis, res):
    f = r.split('\t')
    if f[0] == 'ERR': errs[f[1]] = errs.get(f[1], 0) + 1; continue
    filas.append([s, f[0], f[2], f[3], f[4], f[1]])
n, bad = validar(filas, [1, 2, 3, 4], {1: 'pin', 2: 'trad', 3: 'sis79', 4: 'clase'})
print('estructuras', len(smis), 'nombradas', len(filas), 'nombres', n, 'discrepancias', bad, '| desplazamiento de enlaces:', DESPLAZ[0])
for e, c in sorted(errs.items(), key=lambda x: -x[1]): print('  fuera de dominio', c, e)
inp = []; exp = []
for f in filas:
    m = Chem.MolFromSmiles(f[0])
    for _ in range(3):
        inp.append(Chem.MolToSmiles(m, doRandom=True, canonical=False, kekuleSmiles=random.random() < 0.5)); exp.append((f[1], f[5]))
out = subprocess.run(['php', 'nombrar_org.php'], input='\n'.join(inp) + '\n', capture_output=True, text=True).stdout.rstrip('\n').split('\n')
inv = 0
for e, g, s in zip(exp, out, inp):
    gg = g.split('\t')
    if (gg[0], gg[1] if len(gg) > 1 else '') != e: inv += 1; print('NO-INVARIANTE', s, e, gg[:2])
print('ordenaciones aleatorias', len(inp), 'no invariantes', inv)
