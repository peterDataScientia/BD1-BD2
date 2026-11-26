#!/bin/bash
# Run on login1; source this script, not as a subprocess.
module load python/3.9.6
python -m venv ~/brd4_protonation_env
source ~/brd4_protonation_env/bin/activate
python -m pip install 'propka==3.5.1' 'pdb2pqr==3.6.2' 'urllib3==1.26.20'
python -m pip check
python -c 'import ssl,urllib3,requests; print(ssl.OPENSSL_VERSION,urllib3.__version__,requests.__version__)'
propka3 --version
pdb2pqr --version
