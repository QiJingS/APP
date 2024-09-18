source("packages.R")
source("APP_extraction.R")
APP_extraction("missense")
APP_SNP_summary = readxl::read_xlsx(APP_xlsx_address)
VEP_annotation = c("missense")
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
    cat(i,gene)
    vcf_file_add = paste0(res_address,genes[gene],VEP_annotation,".vcf")
    extract = snp_extraction_combination$new()
    extraction_VCF = extract$data_reorgan(vcf_file_add,vcfR::read.vcfR(vcf_file_add)@fix[,"ID"])
    genotypes_numerical_res = extract$genotypes_numerical()
    genotypes = extract$genotypes_combination()
    
    write.table(genotypes_numerical_res$additive_model,file = paste0(allotype_res_address,genes[gene],"_additive.txt"))
    write.table(genotypes_numerical_res$dominant_model,file = paste0(allotype_res_address,genes[gene],"_dominate.txt"))
    write.table(genotypes_numerical_res$recessive_model,file = paste0(allotype_res_address,genes[gene],"_recessive.txt"))
    write.table(genotypes_numerical_res$overdominant_model,file = paste0(allotype_res_address,genes[gene],"_overdominate.txt"))
    write.table(genotypes[[2]],file = paste0(allotype_res_address,genes[gene],"_allotypes.txt"))
  }
}

