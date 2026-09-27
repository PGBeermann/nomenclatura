# Estructuras aleatorias arbitrarias (como las que podría dibujar un estudiante):
# árboles de 2-26 C y monociclos con sustituyentes arbitrarios. Round-trip con OPSIN (3 estilos).
import random, subprocess, os, sys
from rdkit import Chem
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
JAR=os.environ.get('OPSIN_JAR','opsin_cli.jar'); random.seed(int(sys.argv[2]) if len(sys.argv)>2 else 1)
def aleatoria():
    n=random.randint(1,26); rw=Chem.RWMol(); deg=[]
    ring = random.random()<0.4
    if ring:
        k=random.randint(3,10)
        for i in range(k): rw.AddAtom(Chem.Atom(6)); deg.append(0)
        for i in range(k): rw.AddBond(i,(i+1)%k,Chem.BondType.SINGLE); deg[i]+=1; deg[(i+1)%k]+=1
    else:
        rw.AddAtom(Chem.Atom(6)); deg.append(0)
    while rw.GetNumAtoms()<n+(k if ring else 0) - (0 if ring else 0) and rw.GetNumAtoms()<30:
        cand=[i for i in range(rw.GetNumAtoms()) if deg[i]<4]
        a=random.choice(cand); b=rw.AddAtom(Chem.Atom(6)); deg.append(0)
        rw.AddBond(a,b,Chem.BondType.SINGLE); deg[a]+=1; deg[b]+=1
    if random.random()<0.6:
        for _ in range(random.randint(1,5)):
            cand=[i for i in range(rw.GetNumAtoms()) if rw.GetAtomWithIdx(i).GetAtomicNum()==6 and deg[i]<4]
            if not cand: break
            a=random.choice(cand); x=rw.AddAtom(Chem.Atom(random.choice([9,17,35,53]))); deg.append(1)
            rw.AddBond(a,x,Chem.BondType.SINGLE); deg[a]+=1
    m=rw.GetMol(); Chem.SanitizeMol(m); return Chem.MolToSmiles(m)
N=int(sys.argv[1]); smis=list({aleatoria() for _ in range(N)})
res=subprocess.run(['php','nombrar_smiles.php'],input='\n'.join(smis)+'\n',capture_output=True,text=True).stdout.strip().split('\n')
ok=[(s,r.split('\t')) for s,r in zip(smis,res) if not r.startswith('ERR')]
errs=[(s,r) for s,r in zip(smis,res) if r.startswith('ERR')]
names=[f[k] for k in (0,2,3) for s,f in ok]
out=subprocess.run(['java','-jar',JAR,'-osmi'],input='\n'.join(names)+'\n',capture_output=True,text=True).stdout.split('\n')
out=[o for o in out if not o.startswith(('Picked','Run the'))]
can=lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s)) if s and Chem.MolFromSmiles(s) else None
bad=0
for j,k in enumerate((0,2,3)):
    for i,(s,f) in enumerate(ok):
        if can(out[j*len(ok)+i])!=can(s): bad+=1; print('MISMATCH',s,f[k],out[j*len(ok)+i])
print('estructuras',len(smis),'nombradas',len(ok),'fuera de dominio',len(errs),'nombres',3*len(ok),'discrepancias',bad)
for s,r in errs[:3]: print('  ',s,r)
