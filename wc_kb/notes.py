sudo docker run -v /media/sf_vm_share/wc/:/media/sf_vm_share/wc/ -it karrlab/wc_env_dependencies

pip install -e /media/sf_vm_share/wc/wc_utils/
pip install -e /media/sf_vm_share/wc/obj_model/
pip install -e /media/sf_vm_share/wc/wc_lang/
pip install -e /media/sf_vm_share/wc/wc_kb/
