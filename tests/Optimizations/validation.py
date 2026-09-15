import argparse,hashlib,json,shutil,subprocess,time,zipfile
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description='Run from a real Windows console: GRAMM uses console cursor APIs.')
    p.add_argument('--baseline',type=Path,required=True);p.add_argument('--patched',type=Path,required=True)
    p.add_argument('--sample-computation',type=Path,required=True,help='GUI/SampleProjects/AscendingBridge/Computation')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    root=a.output.resolve();root.mkdir(parents=True,exist_ok=False);reports=[]
    for custom in [False,True]:
        pair=[]
        for label,dll in [('baseline',a.baseline.resolve()),('patched',a.patched.resolve())]:
            folder=root/f'{label}_{custom}';folder.mkdir()
            for name in ['ggeom.asc','landuse.asc','GRAMM.geb','IIN.dat']:shutil.copy2(a.sample_computation/name,folder/name)
            lines=(folder/'IIN.dat').read_text(encoding='utf-8-sig').splitlines()
            for line,value in {2:5,3:60,4:1,13:37.46}.items():lines[line]=lines[line].split(':')[0]+': '+str(value)+' ! regression fixture'
            (folder/'IIN.dat').write_text('\n'.join(lines)+'\n',encoding='utf-8')
            (folder/'Max_Proc.txt').write_text('1\n');(folder/'GRAMMin.dat').write_text('Version 17.01\ny\n0.2\n1,0\nyes\n0\n')
            met=['27,2,4,0.25','9,1,7,0.25','18,3,2,0.5'] if custom else ['27,2,4,1']
            (folder/'meteopgt.all').write_text('10,1,10\ndirection,speed,stability,frequency\n'+'\n'.join(met)+'\n')
            if custom:
                extra=['1000,267.65,270.15,274.4,-1.5,0.85,250,-0.0065,0.012,-0.004','2000,302.9,303.65,293.15,1.5,0.55,1000,-0.005,0.009,-0.006','0,285.65,0,0,0,0.92,0,0,0,0']
                (folder/'CustomInit.txt').write_text('GRAMM custom initial conditions; values at sea level\nheader\n'+'\n'.join(x+','+y for x,y in zip(met,extra))+'\n')
            start=time.perf_counter();subprocess.run(['dotnet',str(dll),str(folder),'1',str(len(met))],cwd=folder,check=True,timeout=180)
            elapsed=time.perf_counter()-start;hashes={}
            for f in sorted(folder.glob('000*')):
                if zipfile.is_zipfile(f):
                    with zipfile.ZipFile(f) as z:
                        for n in z.namelist():hashes[f.name+'/'+n]=hashlib.sha256(z.read(n)).hexdigest()
                else:hashes[f.name]=hashlib.sha256(f.read_bytes()).hexdigest()
            assert len(hashes)==len(met)*5,'missing wind/scalar/steady-state outputs'
            pair.append(hashes);reports.append(dict(custom=custom,label=label,seconds=elapsed,sha256=hashlib.sha256(dll.read_bytes()).hexdigest(),payloads=hashes))
        assert pair[0]==pair[1],'numerical output regression (ZIP timestamps are ignored)'
    (root/'validation.json').write_text(json.dumps(dict(status='pass',runs=4,numerical_payloads=20,reports=reports),indent=2));print('PASS GRAMM 20 numerical payloads')
if __name__=='__main__':main()
