import subprocess
from pathlib import Path
import os

# mkdir example_result
if os.path.exists('example_results'):
    pass
else:
    subprocess.run(['mkdir example_results'], shell=True)

for files in os.listdir('test_images'):
    inp = Path('test_images', files)
    out   = Path('example_results', files)
    if os.path.exists(out):
        pass
    else:
        subprocess.run(['mkdir %s' %out], shell=True)

    subprocess.run(['python3', 'run.py', '-input', inp, '--outdir', out, '--show', 'True'])
      
