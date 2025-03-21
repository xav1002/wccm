import os
import subprocess

import pandas as pd

substrates = ['Acetate','Acetaldehyde','NADH','NAD+']
SMILES = [
    'CC(=O)[O-]',
    'CC=O',
    'C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O',
    'C1=CC(=C[N+](=C1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)([O-])OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N'
    ]
sequences = [
    'maenhdyereinrlfelqkknvvrlrtssideriaklkklkeyiwenkekiqeavyndlrkppeevllteiypvvseirhviknlkkwtkpkkvrtpislfgaksyyrfeakgvvliispwnypfelsigplitaiaagnavvlkpselsphtsgyikklvadifdesevavvegdavvaqkllemgfnhifftgstkvakavlkkasetlssvtlelggkspviidgkfdieeaakkitwgkylnagqtciapdyvfvkkellgdfvshlkhyikkyyysdgsgrcsnycgiinerhfnrlknvfevtvkegakvcegglfvenecyisptvltdvgrdsyimeeeifgpilpvltyekiddvieyinskpaplvlyvfsrdrkfyrhvinnvisgdclindviahfanprlpfgghnasgigkshgyygfrefshlrsimiqpkrtmlqllyppygefvkkliewstkyf',
    'maenhdyereinrlfelqkknvvrlrtssideriaklkklkeyiwenkekiqeavyndlrkppeevllteiypvvseirhviknlkkwtkpkkvrtpislfgaksyyrfeakgvvliispwnypfelsigplitaiaagnavvlkpselsphtsgyikklvadifdesevavvegdavvaqkllemgfnhifftgstkvakavlkkasetlssvtlelggkspviidgkfdieeaakkitwgkylnagqtciapdyvfvkkellgdfvshlkhyikkyyysdgsgrcsnycgiinerhfnrlknvfevtvkegakvcegglfvenecyisptvltdvgrdsyimeeeifgpilpvltyekiddvieyinskpaplvlyvfsrdrkfyrhvinnvisgdclindviahfanprlpfgghnasgigkshgyygfrefshlrsimiqpkrtmlqllyppygefvkkliewstkyf',
    'maenhdyereinrlfelqkknvvrlrtssideriaklkklkeyiwenkekiqeavyndlrkppeevllteiypvvseirhviknlkkwtkpkkvrtpislfgaksyyrfeakgvvliispwnypfelsigplitaiaagnavvlkpselsphtsgyikklvadifdesevavvegdavvaqkllemgfnhifftgstkvakavlkkasetlssvtlelggkspviidgkfdieeaakkitwgkylnagqtciapdyvfvkkellgdfvshlkhyikkyyysdgsgrcsnycgiinerhfnrlknvfevtvkegakvcegglfvenecyisptvltdvgrdsyimeeeifgpilpvltyekiddvieyinskpaplvlyvfsrdrkfyrhvinnvisgdclindviahfanprlpfgghnasgigkshgyygfrefshlrsimiqpkrtmlqllyppygefvkkliewstkyf',
    'maenhdyereinrlfelqkknvvrlrtssideriaklkklkeyiwenkekiqeavyndlrkppeevllteiypvvseirhviknlkkwtkpkkvrtpislfgaksyyrfeakgvvliispwnypfelsigplitaiaagnavvlkpselsphtsgyikklvadifdesevavvegdavvaqkllemgfnhifftgstkvakavlkkasetlssvtlelggkspviidgkfdieeaakkitwgkylnagqtciapdyvfvkkellgdfvshlkhyikkyyysdgsgrcsnycgiinerhfnrlknvfevtvkegakvcegglfvenecyisptvltdvgrdsyimeeeifgpilpvltyekiddvieyinskpaplvlyvfsrdrkfyrhvinnvisgdclindviahfanprlpfgghnasgigkshgyygfrefshlrsimiqpkrtmlqllyppygefvkkliewstkyf'
    ]
sequences = list(map(lambda x: x.upper(),sequences))
pdb_paths = ['seq10','seq10','seq10','seq10']

test = pd.DataFrame({'substrate':substrates,'SMILES':SMILES,'sequence':sequences,'pdbpath':pdb_paths})

print(test)

test.to_csv('./src/catpred_pipeline/pipeline_entry_kcat.csv')
test.to_csv('./src/catpred_pipeline/pipeline_entry_km.csv')

# cmd = 'conda activate catpred'
# subprocess.call(cmd, shell=True, executable='/bin/bash')

cmd = 'cd src/catpred_pipeline/catpred;'
os.system(cmd)
subprocess.run(cmd,shell=True)

cmd = 'cd src/catpred_pipeline/catpred; conda run -n catpred python3 ./demo_run.py --parameter kcat --input_file ../pipeline_entry_kcat.csv --checkpoint_dir ../data/pretrained/production/kcat/ --use_gpu'
subprocess.call(cmd,shell=True,executable='/bin/bash')