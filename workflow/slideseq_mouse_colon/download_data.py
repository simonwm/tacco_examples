import os

files_to_download = [ value.split('/')[-1] for key, value in snakemake.output.items() ]
for key, value in snakemake.output.items():
    folder_to_put_them = os.path.abspath('./' + '/'.join(value.split('/')[:-1]))
    break

msg = f'As it is currently not possible to automatically download data from the Single Cell Portal without user interaction, please download the following files from the https://singlecell.broadinstitute.org/single_cell/study/SCP2038 and put them into {folder_to_put_them!r} [the Single Cell Portal can create a personalized bulk download command for the commandline for you]: {files_to_download!r}'

print(msg)
raise RuntimeError(msg)
