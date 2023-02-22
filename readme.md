# Promline

## recquirements 
- a proper installation of miniconda3 (`~/miniconda3`)
- [Dorado](https://github.com/nanoporetech/dorado) should be in the PATH.
- [PEPPER-Margin-DeepVariant](https://github.com/kishwarshafin/pepper) Docker install, r0.8 or r0.8-gpu if GPU available and if PMDV calling wanted. Otherwise, calling with Clair3 is available trough the conda environment.
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
    ```
- Make the script executable:
    ```
    chmod +x pipeline.sh
    ``` 

## run the pipeline
```
USAGE: pipeline.sh [flags] args
flags:
  -w,--base:  working directory (default: '.')
  -s,--sample:  sample name (default: 'JohnDoe')
  -R,--ref:  reference genome (default: '')
  -f,--fast5:  directory containing raw fast5 (default: '')
  -p,--pod5:  directory for  fast5 to pod5 conversion output or already containing pod5 (default: '')
  -q,--fastqs:  basecalling directory from guppy, in case of real time basecalling during sequencing. This pipeline will only use pass reads. The sequencing summary should be in this directory. If one fastq file is provided, instead of a guppy basecalling directory, a copy of the sequencing summary should be placed in the base directory (-w)
                (default: '')
  -S,--summary:  path to sequencing_summary.txt generated during real time sequencing. Provide a sequencing_summary.txt if you start the pipeline with a bam
                 (default: '')
  -k,--bam:  path to aligned BAM. File has to be indexed. (default: '')
  -c,--snp_caller:  snp caller: either pmdv, clair3 or all (default: 'clair3')
  -B,--[no]basecalling:  basecalling (default: false)
  -M,--[no]modified:  modified bases calling (default: false)
  -A,--[no]alignment:  alignment (default: false)
  -r,--flowcell:  flowcell and basecalling model: either r9 or r10 (default: 'r9')
  -b,--bps:  bases called per second (default: '400')
  -a,--acc:  basecalling accuracy (default: 'sup')
  -d,--[no]duplex:  duplex basecalling (default: false)
  -t,--threads:  max number of threads to use (default: '')
  -v,--[no]version:  tools versions (default: false)
  -h,--help:  show this help (default: false)

```

Run as: 
```
pipeline.sh -w ./workingdirectory \
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
  - sniffles
    - sniffles_BND.vcf
    - noQC_snifles.vcf
- QC.html


## Promline's steps
### 1. Conversion to pod5 
pod5 is the last file format for storing nanopore sequencing data. It takes less space than fast5 and improves read/write performance.
By providing `-f` and `-p`, the first step of this pipeline is to convert fast5 to pod5, if the pod5 folder doesn't exist yet or if it is empty.


### 2. basecalling with [`dorado`](https://github.com/nanoporetech/dorado)
You can select basecalling speed during the run, basecalling accuracy, and if modified basecalling is wanted.    
Stereo duplex basecalling is now also supported.    
Stereo duplex and modified basecalling at the same time is not supported yet.   


### 3. alignment 
[`minimap2`](https://github.com/lh3/minimap2)
Fastqs can be provided (with `-q`), if basecalling was done irt on the sequencer for example. 

### 4. QC 
[`pycoQC`](https://github.com/a-slide/pycoQC)   
If fastqs are provided, the pipeline will look into the fastqs directory for a sequencing summary to do a quality check.

A sequencing summary can also be provided with `-S`.

Otherwise, the pipeline will try to reconstitute one if an unaligned dorado bam file was generated before.


### 5. Structural variants calling
[`Sniffles`](https://github.com/fritzsedlazeck/Sniffles) is the structural variant caller recommended by nanopore.
Sniffles will be run twice, with and without the quality filtering of variants.    

the `tandem-repeats` file is a tandem repeat annotations file that can be used by `Sniffles` to improve variant calling in repetitive regions. This file can be found [here](https://github.com/fritzsedlazeck/Sniffles/tree/master/annotations) and is also in this directory.   

After SV calling, BND variants will be duplicate to create lines in the VCF with their BND mates coordinates, so both breakpoints can be represented and egally annotated in following steps.

### 6. Small variants calling
Small variants calling can either be done with `Clair3`, the caller recommended by nanopore, or with PEPPER-Margin-DeepVariant, or both.

#### [`PMDV`](https://github.com/kishwarshafin/pepper)
