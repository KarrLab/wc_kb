from wc_kb import kb_translater
import pdb

core_path = '/media/sf_vm_share/wc/wc_kb/kbs/original_core.xlsx'
kb_path = '/media/sf_vm_share/wc/wc_kb/kbs/kb_core_empty.xlsx'

translator = kb_translater.KbTranslater(core_path = core_path, kb_path = kb_path)

#translator.translateChromosomeFeatures
#translator.translateTranscriptionUnits()
#translator.translateGenes()
#translator.translateReferences()
translator.translateMetabolites()

translator.kb.save('/media/sf_vm_share/wc/result_kb.xlsx')
kb_reload = wc_kb.io.Reader().run(core_path = '/media/sf_vm_share/wc/result_kb.xlsx',
                                seq_path = '/media/sf_vm_share/wc/wc_kb/kbs/seq.fna')
