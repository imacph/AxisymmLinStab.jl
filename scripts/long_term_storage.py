import pandas as pd
import os
import sys
import sqlite3
import uuid
from subprocess import call

def cp_dir(source, target):

    print(source,target)
    #call(['cp', '-a', source, target]) # Linux

def main():
    conn = sqlite3.connect('spring_2025_AxisymmLinStab.db')
    dn = os.path.dirname(os.path.realpath(__file__))

    task_id = int(sys.argv[1])

    folder_uuid = str(uuid.uuid4())
    df = pd.read_csv(dn+"/array_config.csv",index_col='task_id',dtype=object)
    

    row = df.loc[[task_id]].copy()

    print(folder_uuid)

    row['folder_uuid'] = folder_uuid
    
    print(row)

    row.to_sql('run_info', conn, if_exists='append',index=False)

    
    conn.close()

    os.mkdir(dn+'/'+folder_uuid)
    cp_dir(dn+"/folder_1/.",dn+'/'+folder_uuid)
    return 0


main()