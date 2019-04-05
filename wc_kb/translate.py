from wc_kb import kb_translater
import wc_kb
import pdb

core_path = '/media/sf_vm_share/wc/wc_kb/kbs/original_core.xlsx'
kb_path = '/media/sf_vm_share/wc/wc_kb/kbs/mycoplasmaPneumonia_KB_backup.xlsx'

translator = kb_translater.KbTranslater(core_path = core_path, kb_path = kb_path)

#translator.translateChromosomeFeatures() # included - need t remove _MPNXXX from DnaBox types
#translator.translateTranscriptionUnits() # inclued check after genes for ob mapping
#translator.translateGenes() # included
#translator.translateReferences() # included
translator.translateMetabolites()

translator.kb.save('/media/sf_vm_share/wc/wc_kb/kbs/mycoplasmaPneumonia_KB2.xlsx')
kb_reload = wc_kb.io.Reader().run(core_path = '/media/sf_vm_share/wc/wc_kb/kbs/mycoplasmaPneumonia_KB2.xlsx',
                                seq_path = '/media/sf_vm_share/wc/wc_kb/kbs/seq.fna')

#kb_reload = wc_kb.io.Reader().run(core_path = '/media/sf_vm_share/wc/wc_kb/kbs/mycoplasmaPneumonia_KB.xlsx',
#                                  seq_path = '/media/sf_vm_share/wc/wc_kb/kbs/seq.fna')
