""" Initial pickle creation """
import wc_kb

kb = wc_kb.io.Reader().run(core_path='/home/balazs/Desktop/mycoplasma_pneumoniae/mycoplasma_pneumoniae/kb/core_wComplexes.xlsx',
                              seq_path= '/home/balazs/Desktop/mycoplasma_pneumoniae/mycoplasma_pneumoniae/kb/seq.fna')



""" Example script on how to build wc_lang models from wc_kb using model generator
import wc_kb
import wc_lang
import pickle

kb = pickle.load(open('kb.pickle','rb'))
model = wc_lang.model_gen.ModelGenerator(knowledge_base=kb, version='0.1').run()

mycoplasma_pneumoniae.model_gen.CompartmentsGenerator(kb,model).run()
mycoplasma_pneumoniae.model_gen.ParametersGenerator(kb,model).run()
mycoplasma_pneumoniae.model_gen.MetaboliteSpeciesGenerator(kb,model).run()

mycoplasma_pneumoniae.model_gen.TranscriptionSubmodelGenerator(kb,model).generate_species()
mycoplasma_pneumoniae.model_gen.TranscriptionSubmodelGenerator(kb,model).generate_reactions()


wc_lang.io.Writer().run(model,'/home/balazs/Desktop/test_wc_lang2.xlsx')
"""
