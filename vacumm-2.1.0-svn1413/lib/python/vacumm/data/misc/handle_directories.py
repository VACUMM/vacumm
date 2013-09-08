# -*- coding: utf8 -*-
"""Variables and utilities about colors and color maps


    

"""

import os



def make_directories(SCRIPT_DIR,  WK_DIR,  dir1, dir2):
    # ########################
    #  debut I 
    # changement de repertoire et creation eventuelle de
    # nouveaux repertoires
    # les fichiers modeles sont rapatries dans vacumm/work/MODEL/MARS
    
    # repositionnement a vacumm/work
    #os.chdir(SCRIPT_DIR)
    os.chdir(WK_DIR)
    WORKDIR=os.getcwd()
    
    if os.path.isdir(WORKDIR+'/'+dir1)==False:
        os.mkdir(dir1)
    
    # repositionnement a vacumm/work/dir1
    os.chdir(dir1)
    
    # creation du repertoire vacumm/work/dir1/dir2
    if os.path.isdir(WORKDIR+'/'+dir1+'/'+dir2)==False:
        os.mkdir(dir2)
    
    # repositionnement a vacumm/work/dir1/dir2
    os.chdir(dir2)

    
