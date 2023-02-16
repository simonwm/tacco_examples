import urllib.request
import io
import gzip
import shutil
import tarfile

def get_slideseq_urls(replicate,slide):

    overview_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169012/miniml/GSE169012_family.xml.tgz'
    with urllib.request.urlopen(overview_url) as f:
        response = f.read()
    with io.BytesIO(response) as f:
        with tarfile.open(fileobj=f, mode="r:gz") as tar:
            with tar.extractfile(tar.getnames()[0]) as xmlfile:
                content = xmlfile.read().decode('utf-8')
    available_files = [c.replace('ftp://','https://') for c in content.split('\n') if c.startswith('ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5173nnn/')]
    
    slide = int(slide)
    replicate = int(replicate)

    location_url = [f for f in available_files if f'OB{replicate}_{slide:02d}_BeadLocationsForR.csv.gz' in f]
    expression_url = [f for f in available_files if f'OB{replicate}_{slide:02d}_expression_matrix.mtx.gz' in f]
    rds_url = [f for f in available_files if f'OB{replicate}_Slide{slide}.rds.gz' in f]
    
    if len(location_url) == 0 or len(expression_url) == 0 or len(rds_url) == 0:
        raise ValueError(f'There is no data available for replicate {replicate!r} and slide {slide!r}!')
    
    if len(location_url) > 1 or len(expression_url) > 1 or len(rds_url) > 1:
        raise ValueError(f'There is no unique data available for replicate {replicate!r} and slide {slide!r}!')
    
    return {
        'location_url': location_url[0],
        'expression_url': expression_url[0],
        'rds_url': rds_url[0],
    }

urls = get_slideseq_urls(snakemake.wildcards['replicate'],snakemake.wildcards['slide'])

# there are inconsistencies in puck 1_2 between the gene names (only available in the rds) and the count matrix from mtx - so extract everything from the rds
# 
#with urllib.request.urlopen(urls['location_url']) as f:
#    with open(snakemake.output['location'], 'wb') as g:
#        shutil.copyfileobj(f, g)
#
#with urllib.request.urlopen(urls['expression_url']) as f:
#    with open(snakemake.output['expression'], 'wb') as g:
#        shutil.copyfileobj(f, g)
#
with urllib.request.urlopen(urls['rds_url']) as f:
    with gzip.GzipFile(mode='r', fileobj=f) as gz:
        with open(snakemake.output['rds'], 'wb') as g:
            shutil.copyfileobj(gz, g)

