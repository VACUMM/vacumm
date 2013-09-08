#!/usr/bin/env python
# -*- coding: utf8 -*-
"""Travail avec fichiers distants (:mod:`vacumm.misc.remote`) :


- Fichiers d'entreé : télécharger automatique les fichiers nécessaires
  en fonction de leur disponibilité locale et de leur date.
- Fichier de sortie : envoyer un fichier de sortie après son écriture.


"""


from vcmq import os, InputWorkFiles, OutputWorkFile

host = 'caparmor-sftp'

# Materiel
if not os.path.exists('local'): os.mkdir('local')
if not os.path.exists('remote'): os.mkdir('remote')
rfiles = ['remote/data1.txt', 'remote/data2.txt']
for rfile in rfiles:
    if not os.path.exists(rfile):
        f = open(rfile, 'w')
        f.write('data\n')
        f.close()
wdir = os.getcwd()


# Fichiers d'entrée
# - init et recuperation
wfi = InputWorkFiles('(sftp://%(host)s%(wdir)s/remote>local)/data?.txt'%vars())
print wfi.remote_files()
print wfi.local_files()
# - update
wfi.get()                                   # ESSAYER APRES UN TOUCH, PUIS PARAM CHECK=..


# Fichier de sortie
wfo = OutputWorkFile('(local>sftp://%(host)s%(wdir)s/remote)/out.dat'%vars())
# - ecriture
f = open(wfo.local_file, 'w')
f.write('out\n')
f.close()
# - envoi
wfo.put()
