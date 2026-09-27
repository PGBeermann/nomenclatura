import subprocess, sys, random
from rdkit import Chem
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
import os
JAR=os.environ.get('OPSIN_JAR','opsin_cli.jar')  # OPSIN CLI (p. ej. incluido en 'pip install pyopsin')
N=int(sys.argv[1])
rows=[l.split('\t') for l in subprocess.run(['php','volcar.php',str(N)],capture_output=True,text=True).stdout.strip().split('\n')]
K=3
names=[r[1+k] for k in range(K) for r in rows]
out=subprocess.run(['java','-jar',JAR,'-osmi'],input='\n'.join(names)+'\n',capture_output=True,text=True).stdout.split('\n')
out=[o for o in out if not o.startswith(('Picked','Run the'))]
out=out[:K*len(rows)]
assert len(out)==K*len(rows),(len(out),len(rows))
can=lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s)) if s and Chem.MolFromSmiles(s) else None
bad=0
for i,r in enumerate(rows):
    ref=can(r[0])
    for k in range(K):
        nm=r[1+k]
        got=can(out[k*len(rows)+i])
        if got!=ref: bad+=1; print('MISMATCH',r[0],nm,out[k*len(rows)+i])
print('structures',len(rows),'names checked',K*len(rows),'mismatches',bad)
# invariance: random atom orderings must give identical names
uniq=list({can(r[0]):r for r in rows}.values())
inp=[];exp=[]
for r in uniq:
    m=Chem.MolFromSmiles(r[0])
    for _ in range(4):
        inp.append(Chem.MolToSmiles(m,doRandom=True,canonical=False)); exp.append((r[1],r[4]))
res=subprocess.run(['php','nombrar_smiles.php'],input='\n'.join(inp)+'\n',capture_output=True,text=True).stdout.strip().split('\n')
inv=sum(1 for e,g in zip(exp,res) if tuple(g.split('\t')[:2])!=e)
for e,g,s in zip(exp,res,inp):
    if tuple(g.split('\t')[:2])!=e: print('NONINVARIANT',s,e,g)
print('unique',len(uniq),'random orderings',len(inp),'non-invariant',inv)
