cmd='python3.6 migrations/2019_09_20.py'
cmd='wc-kb validate'

# wc_kb
${cmd} pro ~/Documents/wc_kb/tests/fixtures/prokaryote_core.xlsx
${cmd} eu ~/Documents/wc_kb/tests/fixtures/eukaryote_core.xlsx

# h1_hesc
${cmd} eu ~/Documents/h1_hesc/h1_hesc/kb_gen/core.xlsx
${cmd} eu ~/Documents/h1_hesc/tests/code/fixtures/eukaryote_core.xlsx

# mycoplasma_pneumoniae
${cmd} pro ~/Documents/mycoplasma_pneumoniae/mycoplasma_pneumoniae/kb/mycoplasma_pneumoniae_kb.xlsx
${cmd} pro ~/Documents/mycoplasma_pneumoniae/mycoplasma_pneumoniae/kb/mycoplasma_pneumoniae_kb_old_version.xlsx
${cmd} pro ~/Documents/mycoplasma_pneumoniae/mycoplasma_pneumoniae/network\ model/mycoplasma_pneumoniae_kb_wFunctions.xlsx
${cmd} pro ~/Documents/mycoplasma_pneumoniae/mycoplasma_pneumoniae/kb/mycoplasma_kb_wProcess.xlsx
${cmd} pro ~/Documents/mycoplasma_pneumoniae/mycoplasma_pneumoniae/kb/mycoplasmaPneumonia_KB.xlsx
${cmd} pro ~/Documents/mycoplasma_pneumoniae/tests/fixtures/kb/min_model_kb.xlsx

# wc_model_gen
${cmd} pro ~/Documents/wc_model_gen/tests/fixtures/min_model_kb.xlsx
