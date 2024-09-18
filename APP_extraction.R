
APP_extraction = function(VEP_annotation = ""){
  APP_SNP_summary = readxl::read_xlsx(APP_xlsx_address)
  pages = c("proteasome","TAPs","ERAPs","CALR","CANX","TAPBP","ERp57","Beta_2m","JAK2")
  for(i in seq_along(target)){
    index_page = which(pages %in% target[i])
    gene_summary = readxl::read_xlsx(APP_xlsx_address,sheet = index_page)
    if(length(VEP_annotation) != 0){
      gene_summary = gene_summary[which(grepl(paste(VEP_annotation, collapse = "|"), gene_summary$`VEP Annotation`)),]
    }
    genes = names(table(gene_summary$GENE))
    for(gene in seq_along(genes)){
      chr = unique(gene_summary$Chromosome[which(gene_summary$GENE == genes[gene])])
      if(chr != 6){
        next
      }
      SNPs = gene_summary$rsIDs[which(gene_summary$GENE == genes[gene])]
      # change the name here
      vcf = vcfR::read.vcfR(paste0(chr_vcf_address,"chr",chr,".dose.vcf.gz"))
      SNPS = vcf[vcf@fix[, "ID"] %in% SNPs, ]
      af_filteres_index = which(as.numeric(extract.info(SNPS,"AF")) >= AF)
      not_af_filteres_index = which(as.numeric(extract.info(SNPS,"AF")) < AF)
      SNPS_filtered_names = vcf@fix[, "ID"][vcf@fix[, "ID"] %in% SNPs][af_filteres_index]
      not_SNPS_filtered_names = vcf@fix[, "ID"][vcf@fix[, "ID"] %in% SNPs][not_af_filteres_index]
      SNPS_filtered = vcf[vcf@fix[, "ID"] %in% SNPS_filtered_names, ]
      
      SNPS_found_res = data.frame(SNPS_filtered_names,af = as.numeric(extract.info(SNPS,"AF"))[af_filteres_index])
      SNPS_found_but_filtered = data.frame(not_SNPS_filtered_names,af = as.numeric(extract.info(SNPS,"AF"))[-af_filteres_index])
      SNPS_not_found = SNPs[which(!SNPs %in% vcf@fix[, "ID"] )]
      
      write.vcf(SNPS_filtered, file = paste0(res_address,genes[gene],VEP_annotation,".vcf"))
      write.table(c(paste0("SNPS filtered:",paste0(SNPS_found_but_filtered$not_SNPS_filtered_names,collapse = ","),"\n",
                           "AF:",paste0(SNPS_found_but_filtered$af,collapse = ",")),
                    paste0("SNPS not found:",paste0(SNPS_not_found,collapse = ","))),
                  file = paste0(res_address,genes[gene],VEP_annotation,".txt"))
    }
  }
}