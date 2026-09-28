# Fase 2: estructuras aleatorias con dos anillos (fusionados, puente, espiro, ensamblajes, multiplicativas):
# round-trip con OPSIN e invariancia ante el orden de los átomos.
import random, subprocess, sys
from rdkit import Chem
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
from validar_org import validar, DESPLAZ
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
BASES = ['c1ccc2ccccc2c1', 'C1Cc2ccccc2C1', 'C1CCc2ccccc2C1', 'C1CCC2CCCCC2C1', 'C1C=Cc2ccccc21', 'c1ccc2cccccc2c1'.replace('c1ccc2cccccc2c1', 'C1=CC=C2C=CC=CC=C2C=C1'),
         'C1CC2CCC1C2', 'C1CC2CCC1CC2', 'C1CC2CC(C1)C2', 'C1CC2CCC(C1)C2', 'C1CCC2CC2C1', 'C1CC2CC12', 'c1ccc2c(c1)CC2',
         'C1CCC2(CC1)CCCC2', 'C1CC12CCCCC2', 'C1CCC2(CC1)CCCCC2', 'C1CC2(C1)CCC2', 'C1CCC2CCCC2C1', 'C1CC2=CC=CC=C2C1'.replace('C1CC2=CC=CC=C2C1', 'C1CC2CC=CC=C2C1'),
         'c1ccc(cc1)-c1ccccc1', 'c1ccc(Cc2ccccc2)cc1', 'c1ccc(Oc2ccccc2)cc1', 'C1CCC(CC1)C1CCCCC1', 'c1ccc(cc1)C1CCCCC1',
         'c1ccc(CCc2ccccc2)cc1', 'CC(C)(c1ccccc1)c1ccccc1', 'c1ccc(C(=O)c2ccccc2)cc1', 'c1ccc(-c2cccc3ccccc23)cc1', 'C1CCC(C1)C1CCCCC1',
         'c1ccc(Nc2ccccc2)cc1', 'C1=CC2C=CC1C2', 'C1=CC2CC1C=C2', 'O=C1CCCc2ccccc12', 'C1=CCC2CCCCC2=C1']
GR = ['O', 'O', 'N', 'NC', '=O', 'C=O', 'C(=O)O', 'C(=O)OC', 'C(N)=O', 'C#N', 'OC', 'Cl', 'Br', 'F', 'C', 'C', 'CC', 'C(C)C', '[N+](=O)[O-]']
def aleatoria():
    s = random.choice(BASES)
    for _ in range(random.choice([0, 1, 1, 2, 2, 3])):
        mm = Chem.MolFromSmiles(s)
        cand = [a.GetIdx() for a in mm.GetAtoms() if a.GetAtomicNum() == 6 and a.GetTotalNumHs() >= 1]
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
print('estructuras', len(smis), 'nombradas', len(filas), 'nombres', n, 'discrepancias', bad, '| desplazamiento de enlaces (pentaleno/heptaleno):', DESPLAZ[0])
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
