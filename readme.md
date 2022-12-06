# Promline

## recquirements 
- a proper installation of miniconda3
- Guppy 6
- Dorado, only if basecalling with Dorado wanted.
- PEPPER-Margin-DeepVariant Docker install, r0.8 or r0.8-gpu if GPU available and if PMDV calling wanted. Otherwise, calling with Clair3 is available trough the conda environment.
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
    conda env create -f environment.yml
    ```
- Make the script executable:
    ```
    chmod +x promline.sh
    ```

## run the pipeline
```
USAGE: promline/pipeline.sh [flags] args
flags:
  -w,--base:  working directory (default: '.')
  -s,--sample:  sample name (default: 'JohnDoe')
  -r,--ref:  reference genome (default: '')
  -f,--fast5:  directory containing raw fast5 (default: '')
  -p,--pod5:  directory for  fast5 to pod5 conversion output or already containing pod5 (default: '')
  -q,--fastqs:  basecalling directory from guppy, in case of real time basecalling during sequencing. This pipeline will only use pass reads. The sequencing summary should be in this directory. If one fastq file is provided, instead of a guppy basecalling directory, a copy of the sequencing summary should be placed in the base directory (-w)
                (default: '')
  -c,--snp_caller:  snp caller: either pmdv, clair3 or all (meaning both) (default: 'clair3')
  -b,--basecalling:  basecaller: either 'guppy', 'dorado', or 'all' (meaning both) or 'none' (meaning no basecalling, and fastqs are provided in the path defined by -q 
                     (default: 'guppy')
  -m,--model:  flowcell and basecalling model: either r9 or r10 (default: 'r9')
  -t,--trf:  tandem repeat annotation file for sniffles (default: '')
  -h,--help:  show this help (default: false)

```

Run as: 
```
promline.sh -w ./workingdirectory \
    -s samplename \
    -r ref.fa \
    -f fast5 \
    -p pod5 \
    -c all \
    -b all \
    -m r9 \
    -t human_GRCh38_no_alt_analysis_set.trf.bed 2>&1 | tee log.txt
```

`2>&1 | tee log.txt` allow to store the pipeline log to `log.txt`.

Some tools are using GPUs when available.

## Promline

### Conversion to pod5 `-f -p`
pod5 is the last file format for storing nanopore sequencing data. It takes less space tahn fast5 and improves read/write performance.
By providing `-f` and `-p`, the first step of this pipeline would be to convert fast5 to pod5.

```
pod5-convert-from-fast5 $FAST5 pod5/
```

### Basecalling `b`
If asked, basecalling will be performed with `guppy`, `dorado`, or both. 
If one alignment will be made per basecalling tools, only `guppy` data will be used for variant calling and QC.

#### [`guppy`](https://nanoporetech.com/)
```
guppy_basecaller \
    -i $FAST5 \
    -s $GUPPY \
    --records_per_fastq 0 \
    -r \
    -c ${MODEL_GUPPY} \
    --device 'auto' \
    --compress_fastq \
    --num_callers $FLAGS_threads \
    --chunk_size 1000 \
    --gpu_runners_per_device 4 \
    --chunks_per_runner 512 \
    --disable_pings
    
cat $GUPPY/pass/* > $FASTQ
```

#### [`dorado`](https://github.com/nanoporetech/dorado)
```
${dorado}/bin/dorado basecaller -b 256 ${dorado}/${MODEL_DORADO} pod5/ | samtools view -Sh -@ 6 - > $DORADOBAM
```

### alignment 
[`minimap2`](https://github.com/lh3/minimap2)

```
minimap2 -t $FLAGS_threads \
    -ax map-ont \
    $REFMMI \
    --MD \
    -Y \
    $FASTQ | samtools sort -o $BAM 
```

### QC 
[`pycoQC`](https://github.com/a-slide/pycoQC)

```
pycoQC -f $GUPPY/sequencing_summary.txt -a $BAM -o output.html
```

### Structural variants calling
[`Sniffles`](https://github.com/fritzsedlazeck/Sniffles) is the structural variant caller recommended by nanopore.
Sniffles will be run twice, with and without the quality filtering of variants.

```
sniffles -i $BAM \
	--vcf $VCF_SNF \
	--tandem-repeats $0/human_GRCh38_no_alt_analysis_set.trf.bed \
	--reference $REF \
    --long-del-coverage 5 \
    --long-dup-coverage 0.5 \
	-t $FLAGS_threads 
    --no-qc
```

the `tandem-repeats` file is a tandem repeat annotations file that can be used by `Sniffles` to improve variant calling in repetitive regions. This file can be found [here](https://github.com/fritzsedlazeck/Sniffles/tree/master/annotations) and is also in this directory.

After calling, BND variants will be duplicate to their BND mates, so both breakpoints can be represented and egally annotated in following steps.

### Small variants calling
Small variants calling can either be done with `Clair3`, the caller recommended by nanopore, or with PEPPER-Margin-DeepVariant, or both.

#### [`Clair3`](https://github.com/HKU-BAL/Clair3)
``` 
run_clair3.sh \
	    --bam_fn=${BAM} \
	    --ref_fn=${REF} \
	    --threads=$FLAGS_threads \
	    --platform="ont" \
	    --model_path=${CONDA_PREFIX}/bin/models/${MODEL_CLAIR} \
	    --output=${CLAIR_OUT} \
	    --remove_intermediate_dir
```
The following module "takes a Clair3 VCF and a Sniffle2 VCF as inputs. It switches the zygosity from homozygous to heterozygous of a Clair3 called SNP that matches the following two criteria: 1) AF<=0.7, and 2) the flanking 16bp of the SNP is inside one or more SV deletions given in the Sniffle2 VCF." Not working.

```
pypy3 /home/euphrasie/miniconda3/envs/promline/bin/clair3.py SwitchZygosityBasedOnSVCalls \
      --bam_fn ${BAM} \
      --clair3_vcf_input ${CLAIR_OUT}/merge_output.vcf.gz \
      --sv_vcf_input $VCF_SNF \
      --vcf_output ${CLAIR_OUT}_merge_output_switch.vcf \
      --threads ${FLAGS_threads}
```

#### [`PMDV`](https://github.com/kishwarshafin/pepper)
  ```
    docker run --ipc=host \
	    --gpus all \
	    -v "${BAMDIR}":"${BAMDIR}" \
	    -v "${OUTPUT_DIR_PMDV}":"${OUTPUT_DIR_PMDV}" \
	    -v "${REF}":"${REF}" \
	    kishwars/pepper_deepvariant:r0.8-gpu \
	    run_pepper_margin_deepvariant call_variant \
	    -o "${OUTPUT_DIR_PMDV}" \
	    -b "${BAM}" \
	    -f "${REF}" \
	    -p "${PMDV_PREFIX}" \
	    -t "{THREADS}" \
	    -g \
	    ${MODEL_PMDV}
  ```