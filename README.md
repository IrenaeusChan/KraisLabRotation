# KraisLabRotation
Collection of code snippets used during the rotation in Dr. John Krais' Lab

Get the BRCA Germline Files
```
cat TCGA_Metadata/carrier_to_mutect_files.txt | while read line; do
  location=$(echo $line | cut -d' ' -f2);
  echo $location;
  cp /storage1/fs1/krais/Active/IrenaeusChan/Mutect2_Annotation/$location/*.vcf.gz /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/BRCA1_Germline/;
done
```

Move the BRCA2 Germline Files Out of BRCA1
```
grep 'BRCA2' TCGA_Metadata/carrier.tcga.new.brca.filtered30VAF.txt | awk '{print $2}' | while read line; do
  grep $line mutectFileMetadata.csv | cut -d',' -f1;
done | while read line; do
  brca_file=$(grep $line mutectFileMetadata.csv | cut -d',' -f2);
  mv PolThetaAnalysis/BRCA1_Germline/$brca_file PolThetaAnalysis/BRCA2_Germline/;
done
```

Get the TCGA_BRCA Files WITHOUT the Germline Mutations
```
grep 'TCGA_BRCA' mutectFileMetadata.csv | cut -d',' -f2 | while read line; do
  file_name=$(anti_grep $line TCGA_Metadata/BRCA_Germline_VCFs.txt);
  location=$(grep $line mutectFileMetadata.csv | cut -d',' -f1);
  cp Mutect2_Annotation/$location/$file_name PolThetaAnalysis/TCGA_BRCA/ ;
done
```

Do the same for TCGA_OV
```
grep 'TCGA_OV' mutectFileMetadata.csv | cut -d',' -f2 | while read line; do
  file_name=$(anti_grep $line TCGA_Metadata/BRCA_Germline_VCFs.txt);
  location=$(grep $line mutectFileMetadata.csv | cut -d',' -f1);
  cp Mutect2_Annotation/$location/$file_name PolThetaAnalysis/TCGA_OV/;
done
```

Index the VCF Files
```
for vcf in PolThetaAnalysis/BRCA1_Germline/*.vcf.gz; do bsub1 kboltonlab/bst:1.0 tabix $vcf; done
for vcf in PolThetaAnalysis/BRCA2_Germline/*.vcf.gz; do bsub1 kboltonlab/bst:1.0 tabix $vcf; done
for vcf in PolThetaAnalysis/TCGA_BRCA/*.vcf.gz; do bsub1 kboltonlab/bst:1.0 tabix $vcf; done
for vcf in TCGA_OV/*.vcf.gz; do bsub1 kboltonlab/bst:1.0 tabix $vcf; done
```

Run the Analysis Script (MainAPP_IC.py) on all of the files
```
# BRCA1_Germline
for vcf in /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/BRCA1_Germline/*.vcf.gz; do
  CompareGroup="BRCA1_Germline"
  Sample_Name=$(basename $vcf .vcf.gz).analysis
  VCF_File=$vcf
  sed "s|<Sample_Name>|$Sample_Name|g" /storage1/fs1/krais/Active/IrenaeusChan/VCF_Analysis/input_params_template.txt > /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
  sed -i "s|<VCF_File>|$VCF_File|g" /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
  sed -i "s|<Group>|$CompareGroup|g" /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
  bsub4 kboltonlab/krais_vcf_analysis:latest python /app/MainAPP_IC.py --options_file /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
done

#BRCA2_Germline
for vcf in /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/BRCA2_Germline/*.vcf.gz; do
  CompareGroup="BRCA2_Germline"
  Sample_Name=$(basename $vcf .vcf.gz).analysis
  VCF_File=$vcf
  sed "s|<Sample_Name>|$Sample_Name|g" /storage1/fs1/krais/Active/IrenaeusChan/VCF_Analysis/input_params_template.txt > /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
  sed -i "s|<VCF_File>|$VCF_File|g" /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
  sed -i "s|<Group>|$CompareGroup|g" /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
  bsub4 kboltonlab/krais_vcf_analysis:latest python /app/MainAPP_IC.py --options_file /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
done

#TCGA_OV
for vcf in /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/TCGA_OV/*.vcf.gz; do
  CompareGroup="TCGA_OV"
  Sample_Name=$(basename $vcf .vcf.gz).analysis
  VCF_File=$vcf
  sed "s|<Sample_Name>|$Sample_Name|g" /storage1/fs1/krais/Active/IrenaeusChan/VCF_Analysis/input_params_template.txt > /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
  sed -i "s|<VCF_File>|$VCF_File|g" /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
  sed -i "s|<Group>|$CompareGroup|g" /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
  bsub4 kboltonlab/krais_vcf_analysis:latest python /app/MainAPP_IC.py --options_file /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
done

#TCGA_BRCA
for vcf in /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/TCGA_BRCA/*.vcf.gz; do
  CompareGroup="TCGA_BRCA"
  Sample_Name=$(basename $vcf .vcf.gz).analysis
  VCF_File=$vcf
  sed "s|<Sample_Name>|$Sample_Name|g" /storage1/fs1/krais/Active/IrenaeusChan/VCF_Analysis/input_params_template.txt > /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
  sed -i "s|<VCF_File>|$VCF_File|g" /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
  sed -i "s|<Group>|$CompareGroup|g" /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
  bsub4 kboltonlab/krais_vcf_analysis:latest python /app/MainAPP_IC.py --options_file /storage1/fs1/krais/Active/IrenaeusChan/PolThetaAnalysis/InputParams/input_params_${Sample_Name}.txt
done
```

From the analysis files, only grab the INDELs to create a summary file
```
cat header.txt > BRCA1_Germline_InsDel.tsv
for vcf in BRCA1_Germline/*.vcf; do 
  sn=$(basename $vcf .wxs.MuTect2.somatic_annotation.analysis.vcf)
  echo $sn
  grep -v '^#' $vcf | awk -F'\t' -v sn="$sn" '$12=="Ins" || $12=="Del" {print $0"\t"sn}' >> BRCA1_Germline_InsDel.tsv
done

cat header.txt > BRCA2_Germline_InsDel.tsv
for vcf in BRCA2_Germline/*.vcf; do 
  sn=$(basename $vcf .wxs.MuTect2.somatic_annotation.analysis.vcf)
  echo $sn
  grep -v '^#' $vcf | awk -F'\t' -v sn="$sn" '$12=="Ins" || $12=="Del" {print $0"\t"sn}' >> BRCA2_Germline_InsDel.tsv
done

cat header.txt > TCGA_OV_InsDel.tsv
for vcf in TCGA_OV/*.vcf; do 
  sn=$(basename $vcf .wxs.MuTect2.somatic_annotation.analysis.vcf)
  echo $sn
  grep -v '^#' $vcf | awk -F'\t' -v sn="$sn" '$12=="Ins" || $12=="Del" {print $0"\t"sn}' >> TCGA_OV_InsDel.tsv
done

cat header.txt > TCGA_BRCA_InsDel.tsv
for vcf in TCGA_BRCA/*.vcf; do 
  sn=$(basename $vcf .wxs.MuTect2.somatic_annotation.analysis.vcf)
  echo $sn
  grep -v '^#' $vcf | awk -F'\t' -v sn="$sn" '$12=="Ins" || $12=="Del" {print $0"\t"sn}' >> TCGA_BRCA_InsDel.tsv
done
```

# Analysis for Dana-Farber Cancer Institute Ovariant Data
Realign hg19 to hg38
```
for bam in /storage1/fs1/krais/Active/DFCI_OV/BAM_Files/*.bam; do 
  sn=$(basename $bam .bam); 
  bsub8 fredhutch/bwa:0.7.17-samtools-1.10 $scripts/hengli_realign.sh $bam /storage1/fs1/krais/Active/DFCI_OV/HG38_BAMs/${sn}.hg38.bam $HG38_REF_DH; 
done
```

Running Mutect and HaplotypeCaller on CGOV Data
```
for normal_bam in /storage1/fs1/krais/Active/DFCI_OV/HG38_BAMs/*N*.hg38.bam; do
  sn=$(basename $normal_bam N.hg38.bam);
  echo $normal_bam
  bsub4 broadinstitute/gatk:4.2.0.0 /gatk/gatk HaplotypeCaller --java-options "-Xmx12g" --native-pair-hmm-threads 3 -O /storage1/fs1/krais/Active/IrenaeusChan/CGOV_Analysis/HaplotypeCaller/${sn}.haplotypecaller.g.vcf.gz -R $HG38_REF_DH -I ${normal_bam} --read-index ${normal_bam}.bai -ERC GVCF --max-reads-per-alignment-start 0
  ls ${sn}T*.bam | while read tumor_bam; do
    echo $tumor_bam;
    tumor_sn=$(basename $tumor_bam .hg38.bam);
    bsub4 broadinstitute/gatk:4.2.0.0 /gatk/gatk Mutect2 --java-options "-Xmx12g" --native-pair-hmm-threads 3 -O /storage1/fs1/krais/Active/IrenaeusChan/CGOV_Analysis/Mutect2/${tumor_sn}.mutect.vcf.gz -R $HG38_REF_DH -I ${tumor_bam} --read-index ${tumor_bam}.bai -tumor "$tumor_sn" -I ${normal_bam} --read-index ${normal_bam}.bai -normal "$sn" --max-reads-per-alignment-start 0
    bsub4 broadinstitute/gatk:4.2.0.0 /gatk/gatk HaplotypeCaller --java-options "-Xmx12g" --native-pair-hmm-threads 3 -O /storage1/fs1/krais/Active/IrenaeusChan/CGOV_Analysis/HaplotypeCaller/${tumor_sn}.haplotypecaller.g.vcf.gz -R $HG38_REF_DH -I ${tumor_bam} --read-index ${tumor_bam}.bai -ERC GVCF --max-reads-per-alignment-start 0
    done;
done
```

# Analysis for the MG126 Assay
Running Mutect on the 126
```
for dir in $(ls -d /storage1/fs1/krais/Active/SR007150_Krais_MG126_Capture/watchmaker/SA*); do
  sn=$(echo $dir | cut -d'/' -f8);
  echo $sn;
  bsub4 broadinstitute/gatk:4.2.0.0 /gatk/gatk Mutect2 --java-options "-Xmx12g" --native-pair-hmm-threads 3 -O ${sn}.mutect.vcf.gz -R $HG38_REF_DH -I ${dir}/${sn}.cram --read-index ${dir}/${sn}.cram.crai --max-reads-per-alignment-start 0
done
for mutect_vcf in /storage1/fs1/krais/Active/IrenaeusChan/MG126_Analysis/Mutect2/*.vcf.gz; do
  sn=$(basename $mutect_vcf .vcf.gz);
  echo $sn;
  bsub4 broadinstitute/gatk:4.2.0.0 /gatk/gatk FilterMutectCalls -R $HG38_REF_DH -V ${mutect_vcf} -O ${sn}.filtered.vcf.gz
done
```
