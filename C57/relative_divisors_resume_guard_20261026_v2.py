from pathlib import Path
import os,signal,time
monitor=1315713;worker=1315714
expected=[b"/opt/magma-2.30-1/magma.exe",b"-b",b"work/degree8/large_levels/C57_followup/relative_divisors_resume_20261026.m"]
limit=3600;hz=os.sysconf(os.sysconf_names['SC_CLK_TCK'])
while True:
 try:st=Path(f"/proc/{worker}/stat").read_text().rsplit(")",1)[1].split()
 except FileNotFoundError:
  print("Worker already reaped; no signal sent.",flush=True);break
 state=st[0];elapsed=float(Path('/proc/uptime').read_text().split()[0])-int(st[19])/hz
 if state=='Z':
  print("NATIVE_MAGMA_EXIT_CODE",int(st[49]),"elapsed",elapsed,flush=True)
  if Path(f"/proc/{monitor}/cmdline").exists():os.kill(monitor,signal.SIGCONT)
  break
 args=Path(f"/proc/{worker}/cmdline").read_bytes().split(b"\0");assert args[:3]==expected,args
 if elapsed>=limit:
  print("Extended3600secondwallbound reached; stopping only verified owned worker.",flush=True)
  os.kill(worker,signal.SIGTERM)
  if Path(f"/proc/{monitor}/cmdline").exists():os.kill(monitor,signal.SIGCONT)
  break
 time.sleep(10)
