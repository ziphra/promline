# Promline

## recquirements 
- a proper installation of miniconda3 (`~/miniconda3`)
- [Dorado](https://github.com/nanoporetech/dorado) should be in the PATH.
- Docker should run without sudo. To do so: 
  ```
  sudo groupadd docker
  sudo usermod -aG docker username
  su - username
  docker run hello-world
  ``` 
## install 
- git clone this repository
    ```
    git clone https://github.com/ziphra/promline
    ```

- Create a `conda env` from `promline.yml`
    ```
    cd promline
    conda env create -f promline.yml
    conda activate promline 
    mv dorado_models ${CONDA_PREFIX}/bin/
    mv clair3_models ${CONDA_PREFIX}/bin/
    ```
    The folders dorado_models and clair3_models contain some available basecalling models for Dorado and variant calling models for Clair3. The most up-to-date models can be downloaded from the Dorado GitHub page for Dorado models or from the Rerio GitHub page for Clair3 models.    
- Make the script executable:
    ```
    chmod +x pipeline.sh
    ``` 

## run the pipeline
```
USAGE: pipeline2.sh [flags] args
flags:
  -w,--base:  working directory (default: '.')
  -s,--sample:  sample name (default: 'JohnDoe')
  -R,--ref:  reference genome (default: '')
  -f,--fast5:  directory containing raw fast5 (default: '')
  -p,--pod5:  directory for  fast5 to pod5 conversion output or already containing pod5 (default: '')
  -q,--fastqs:  basecalling directory from guppy or dorado, in case of real time basecalling during sequencing. This pipeline will only use pass reads. The sequencing summary should be in this directory. If one fastq file is provided, instead of a guppy basecalling directory, a copy of the sequencing summary should be placed in the base directory (-w)
                (default: '')
  -S,--summary:  path to sequencing_summary.txt generated during real time sequencing. Provide a sequencing_summary.txt if you start the pipeline with a bam (default: '')
  -k,--bam:  path to aligned BAM. File has to be indexed. (default: '')
  -c,--snp_caller:  snp caller: either pmdv, clair3 or all (default: 'clair3')
  -B,--[no]basecalling:  basecalling (default: false)
  -M,--[no]modified:  modified bases calling (default: false)
  -A,--[no]alignment:  alignment (default: true)
  -r,--flowcell:  flowcell and basecalling model: either r9 or r10 (default: 'r10')
  -b,--bps:  bases called per second (default: '400')
  -a,--acc:  basecalling accuracy (default: 'sup')
  -d,--[no]duplex:  duplex basecalling (default: false)
  -t,--threads:  max number of threads to use (default: '')
  -v,--[no]version:  tools versions (default: false)
  -h,--help:  show this help (default: false)
```

Run as: 
```
pipeline2.sh -w ./workingdirectory \
    -s samplename \
    -R ref.fa \
    -f fast5 \
    -p pod5 \
    -M \
    --bps 260 \
    --acc fast
    -B \
    -A \
    -c clair3 \
    -r r10 \
    -t 32 2>&1 | tee log.txt
```

`2>&1 | tee log.txt` allow to store the pipeline log to `log.txt`.

Some tools will use GPUs when available.

## Output
A folder structured as follow: 
- dorado
  - unaligned bam file
- mmi
  - aligned bam file and its index
- vc
  - clair3
    - merge_output.vcf.gz
    - phased_output.bam
    - ...
  - sniffles
    - sniffles_BND.vcf
    - noQC_snifles.vcf
- QC.html
- sequencing_summary.txt


## Promline's steps
### 1. Conversion to pod5 
The pod5 file format is used for storing nanopore sequencing data. It requires less space than fast5 and offers improved read/write performance. The first step of this pipeline is to convert fast5 files to pod5 format. This conversion is performed if the pod5 folder doesn't exist yet or if it is empty. To initiate the conversion, use the -f and -p options. By providing `-f` and `-p`, the first step of this pipeline will be to convert fast5 to pod5, *if* the pod5 folder doesn't exist yet or if it is empty.


### 2. basecalling with [`dorado`](https://github.com/nanoporetech/dorado)
During the basecalling process, you can choose the basecalling speed, basecalling accuracy, and whether you want modified basecalling. The pipeline now supports stereo duplex basecalling as well. However, please note that stereo duplex and modified basecalling cannot be used simultaneously.  


### 3. alignment 
[`minimap2`](https://github.com/lh3/minimap2)
The minimap2 tool is used for the alignment step. If fastqs are provided, they can be used with the -q option. For example, if basecalling was done directly on the sequencer.

### 4. QC 
[`pycoQC`](https://github.com/a-slide/pycoQC)   
For quality control, the pycoQC tool is utilized. If fastqs are provided, the pipeline will search for a sequencing summary in the fastqs directory to perform the quality check. You can also provide a sequencing summary using the -S option.

### 5. Structural variants calling
[`Sniffles`](https://github.com/fritzsedlazeck/Sniffles) is the structural variant caller recommended by nanopore.
Sniffles will be run twice, with and without the quality filtering of variants.    

The recommended structural variant caller is Sniffles, which will be run twice, once with quality filtering of variants and once without. Additionally, a tandem-repeats file is provided to aid Sniffles in improving variant calling in repetitive regions. The tandem-repeats file can be found [here](https://github.com/fritzsedlazeck/Sniffles/tree/master/annotations) and is also available in the current directory. After SV calling, BND variants will be duplicated to create lines in the VCF with their BND mate's coordinates, allowing both breakpoints to be represented and equally annotated in subsequent steps.


### 6. Small variants calling
Small variants calling can be done using either Clair3, the caller recommended by nanopore, or PEPPER-Margin-DeepVariant, or both.

### 7. CNV 
CNV calling is performed using cnvpytor.

