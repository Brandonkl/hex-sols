import os
import numpy as np
import yaml
import copy as cp

dir = os.sys.argv[1] # requires input of directory of input directories (default "data")

for f in os.walk(dir):

    if f[0]== dir: continue #os walk checks the parent directory as the first on the search tree, we only care about subdirectories.

    dat = f[0].split('/')[1].split('-')
    n, d= dat[0], dat[2]
    ang = open(f[0]+"/angle.txt", "w")
    ang.write("acos( "+n+" / "+d+" ) = "+str(round(np.degrees(np.arccos(int(n)/int(d))), 3))) #grabs the respective twist angle acos((3p^2 - q^2)/(3p^2 + q^2))
    ang.close()

    y = yaml.load(open(f[0]+"/layers.yaml", "r+"), Loader=yaml.FullLoader)
    lnew = [] # To be MoS2 Layers
    s = 0
    for l in y["layer"]:

        """If anyone has the misfortune of reading this I appologize, the essence of this repeated spaghetti code is takintg each 2 atom frac-lattice for graphene,
        and duplicating one of the atoms to be the sulfur atoms, while the other is the molybdenum. Due to the geometry of MoS2, this requires 3 rsp2 layers per actual layer"""

        lnew.append({'frac-lattice': cp.deepcopy(l['frac-lattice']), 'frac-sites': [cp.deepcopy(l['frac-sites'][1-s%2])], 'repeat': cp.deepcopy(l['repeat']), 'elements': ['S']})
        lnew.append({'frac-lattice': cp.deepcopy(l['frac-lattice']), 'frac-sites': [cp.deepcopy(l['frac-sites'][0+s%2])], 'repeat': cp.deepcopy(l['repeat']), 'elements': ['Mo']})
        lnew.append({'frac-lattice': cp.deepcopy(l['frac-lattice']), 'frac-sites': [cp.deepcopy(l['frac-sites'][1-s%2])], 'repeat': cp.deepcopy(l['repeat']), 'elements': ['S']})
        s += 1

    y["layer"] = lnew
    yaml.dump(y, open(f[0]+"/layers.yaml", "w"))
