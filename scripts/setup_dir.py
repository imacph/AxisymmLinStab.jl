import os
import sys
dn = os.path.dirname(os.path.realpath(__file__))

file_id = int(sys.argv[1])
if not os.path.exists(dn+'/folder_{:d}'.format(file_id)):
    os.makedirs(dn+'/folder_{:d}'.format(file_id))


