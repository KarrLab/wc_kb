""" Migration of WC-KB-encoded files

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-10-10
:Copyright: 2019, Karr Lab
:License: MIT
"""

import os.path
import sys
sys.path.insert(0, 'migrations')
import migration_2020_03_09 as migration

base_dir = os.path.expanduser('~/Documents')

paths = [
    # wc_kb
    {'taxon': 'pro', 'path': 'wc_kb/tests/fixtures/prokaryote_core.xlsx'},
    {'taxon': 'eu', 'path': 'wc_kb/tests/fixtures/eukaryote_core.xlsx'},    

    # h1_hesc
    {'taxon': 'eu', 'path': 'h1_hesc/h1_hesc/kb_gen/core.xlsx'},
    {'taxon': 'eu', 'path': 'h1_hesc/tests/code/fixtures/eukaryote_core.xlsx'},
    {'taxon': 'eu', 'path': 'h1_hesc/h1_hesc/scaled_down_model/kb_core.xlsx'},

    # mycoplasma_pneumoniae
    {'taxon': 'pro', 'path': 'mycoplasma_pneumoniae/mycoplasma_pneumoniae/kb/mycoplasma_pneumoniae_kb.xlsx'},
    {'taxon': 'pro', 'path': 'mycoplasma_pneumoniae/mycoplasma_pneumoniae/kb/mycoplasma_pneumoniae_kb_old_version.xlsx'},
    {'taxon': 'pro', 'path': 'mycoplasma_pneumoniae/mycoplasma_pneumoniae/network model/mycoplasma_pneumoniae_kb_wFunctions.xlsx'},
    {'taxon': 'pro', 'path': 'mycoplasma_pneumoniae/mycoplasma_pneumoniae/kb/mycoplasma_kb_wProcess.xlsx'},
    {'taxon': 'pro', 'path': 'mycoplasma_pneumoniae/mycoplasma_pneumoniae/kb/mycoplasmaPneumonia_KB.xlsx'},
    {'taxon': 'pro', 'path': 'mycoplasma_pneumoniae/tests/fixtures/kb/min_model_kb.xlsx'},

    # wc_model_gen
    {'taxon': 'pro', 'path': 'wc_model_gen/tests/fixtures/min_model_kb.xlsx'},
]


for i_path, path in enumerate(paths):
    print('Migrating path {} of {}: {}'.format(i_path + 1, len(paths), path['path']))

    abs_path = os.path.join(base_dir, path['path'])

    # migrate
    migration.transform(abs_path, path['taxon'])
