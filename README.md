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
