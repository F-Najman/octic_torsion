"""Fresh standalone native reproduction; writes only a NEW run directory.

Run from the repository root:
 python3 C2x30/reproduce_C2x30.py NEW_RUN_DIRECTORY
Expected runtime on this machine is about four CPU minutes, <=8GB.
No earlier-degree result or retained session is used.
"""
import argparse,json,subprocess,sys
from pathlib import Path

P=Path(__file__).resolve().parent
ROOT=P.parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('run',type=Path);ap.add_argument('--magma',default='magma');args=ap.parse_args()
    out=args.run.resolve();out.mkdir(parents=True,exist_ok=True)
    stages=[('classes256','A3_D8_RR_PROBE_DONE'),('quartic_filter','X130_D4_FILTER_DONE'),('quartic_allcusps','X130_D4_ALLCUSPS_DONE'),('exceptional_fibres','EXCEPTIONAL_FIBRES_DONE'),('three_quotients','THREE_QUOTIENTS_DONE')]
    for stage,marker in stages:
        wrapper=out/(stage+'_driver.m')
        inputs={'exceptional_fibres':'exceptional_pencil_inputs.m','three_quotients':'exceptional_fibres_data.m'}
        pre=''
        if stage in inputs:pre='load '+json.dumps(str(out/inputs[stage]))+';\n'
        wrapper.write_text('C230Out:='+json.dumps(str(out))+';\n'+pre+'load '+json.dumps(str(P/'src'/(stage+'.m')))+';\n')
        console=out/(stage+'_console.log')
        with console.open('w') as f:
            proc=subprocess.run([args.magma,'-b',str(wrapper)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,timeout=1800)
        text=console.read_text()
        assert proc.returncode==0 and marker in text and 'Runtime error' not in text and 'User error' not in text and 'syntax error' not in text,(stage,text[-3000:])
        print('C230_NATIVE_STAGE_PASS',stage,flush=True)
        if stage=='classes256':
            with (out/'series_analysis.log').open('w') as f:
                subprocess.run([sys.executable,str(P/'analyze_classes.py'),str(out)],cwd=ROOT,check=True,stdout=f,stderr=subprocess.STDOUT)
            print('C230_SERIES_CLASSIFICATION_PASS',flush=True)
    with (out/'integer_replay.log').open('w') as f:
        subprocess.run([sys.executable,str(P/'replay_arithmetic.py'),str(out)],cwd=ROOT,check=True,stdout=f,stderr=subprocess.STDOUT)
    print('C230_FRESH_REPRODUCTION_PASS run='+str(out),flush=True)

if __name__=='__main__':main()
