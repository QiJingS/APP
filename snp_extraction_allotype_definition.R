
snp_extraction_combination = R6::R6Class("snp_extraction_combination",
                                         public = list(
                                           name = NULL,
                                           data = NULL,
                                           data_reorgan = function(data_vcf_address,snp_list){
                                             data_vcf_address = data_vcf_address
                                             data = read.table(data_vcf_address)
                                             self$name = strsplit(readLines(data_vcf_address)[grep("^#CHROM", readLines(data_vcf_address))],"\t")[[1]]
                                             colnames(data) = self$name
                                             data[,-1:-9] = sapply(
                                               data[,-1:-9],
                                               function(x) sapply(strsplit(x, ":"), `[`, 1)
                                             )
                                             self$data = data[match(snp_list,data$ID),]
                                             return(self$data)
                                           },
                                           genotypes_numerical = function(){
                                             data = self$data
                                             if(nrow(data) == 1){
                                               cat("PLEASE ENTER MORE THAN ONE SNP TO COMBINE AS AN ALLOTYPE.")
                                             }else{
                                             vcf = vcfR::read.vcfR(data_vcf_address)
                                             gt_matrix <- extract.gt(vcf, element = "GT")
                                             
                                             gt_numeric <- apply(gt_matrix, 2, function(x) {
                                               x[x == "./."] <- NA  # Treat missing data as NA
                                               sapply(x, function(gt) {
                                                 if (is.na(gt)) return(NA)
                                                 sum(as.numeric(unlist(strsplit(gt, "[|/]"))))  # Sum alleles
                                               })
                                             })
                                             # Dominant Model
                                             # Group 0/0 as 0 (unaffected), and 0/1 or 1/1 as 1 (affected)
                                             dominant_model <- gt_numeric
                                             dominant_model[gt_numeric == 0] <- 0
                                             dominant_model[gt_numeric == 1 | gt_numeric == 2] <- 1
                                             
                                             # Recessive Model
                                             # Group 0/0 and 0/1 as 0 (unaffected), and 1/1 as 1 (affected)
                                             recessive_model <- gt_numeric
                                             recessive_model[gt_numeric == 0 | gt_numeric == 1] <- 0
                                             recessive_model[gt_numeric == 2] <- 1
                                             
                                             # Overdominant Model
                                             # Group 0/0 and 1/1 as 0 (unaffected), and 0/1 as 1 (affected)
                                             overdominant_model <- gt_numeric
                                             overdominant_model[gt_numeric == 0 | gt_numeric == 2] <- 0
                                             overdominant_model[gt_numeric == 1] <- 1
                                             return(list(additive_model = t(gt_numeric),
                                                         dominant_model = t(dominant_model),
                                                         recessive_model = t(recessive_model),
                                                         overdominant_model = t(overdominant_model)))
                                             }
                                           },
                                           
                                           
                                           genotypes_combination = function(){
                                             data = self$data
                                             if(nrow(data) == 1){
                                               cat("PLEASE ENTER MORE THAN ONE SNP TO COMBINE AS AN ALLOTYPE.")
                                             }else{
                                          
                                             combination_snps = data.frame()
                                             for(l in 10:ncol(data)){
                                               for(left_right in 1:2){
                                                 individual = sapply(strsplit(data[,l],split = "\\|"),'[',left_right)
                                                 combination = rep(NA,length(individual))
                                                 for(i in 1:length(individual)){
                                                   if(individual[i] == '0'){
                                                     combination[i] = data$REF[i]
                                                   }else{
                                                     combination[i] = data$ALT[i]
                                                   }
                                                 }
                                                 combination_snps[l-9,left_right] = paste(combination,collapse = "")
                                               }
                                               for(p in 1:nrow(data)){
                                                 combination_snps[l-9,p+2] = paste0(substr(combination_snps[l-9,1],p,p),substr(combination_snps[l-9,2],p,p),collapse = " ")
                                               }
                                             }
                                             rownames(combination_snps) = self$name[-1:-9]
                                             colnames(combination_snps) = c("comb1","comb2",data$ID)
                                             if(nchar(combination_snps[1,1]) == 1){
                                               combination_snps = combination_snps[-1:-2]
                                             }
                                             
                                             data = combination_snps
                                             header = unique(c(names(table(data$comb1)),names(table(data$comb2))))
                                             allotype_individual_number = data.frame()
                                             if(dim(data)[2] >= 2){
                                               for(i in 1:dim(data)[1]){
                                                 allotype_individual_number[i,1:length(header)] = 0
                                                 for(j in 1:2){
                                                   find_index = which(header == data[i,j])
                                                   if(allotype_individual_number[i,find_index] == 1){
                                                     allotype_individual_number[i,find_index] = 2
                                                   }else{
                                                     allotype_individual_number[i,find_index] = 1
                                                   }
                                                 }
                                               }
                                               names(allotype_individual_number) = header
                                               rownames(allotype_individual_number) = rownames(data)
                                             }
                                             #
                                             
                                             #combination_snps[,-1:-2]
                                             for(w in 3:dim(combination_snps)[2]){
                                             get_rep_name = names(which(sapply(names(table(combination_snps[,w])), function(x) length(unique(strsplit(x, "")[[1]]))) == 2))
                                             combination_snps[which(combination_snps[,w] == get_rep_name[2]),w] = get_rep_name[1]
                                             }
                                             return(list(combination_snps = combination_snps,
                                                    allotype_individual_number = allotype_individual_number))
                                            }
                                           },
                                           
                                           allotype_combination = function(amino_acid_address){
                                             amino_acid_changes = read.table(amino_acid_address)
                                             geno = self$genotypes_combination()
                                             data = geno$combination_snps
                                             for(i in 1:nrow(data)){
                                               for(j in 1:2){
                                                 for(p in 1:nrow(amino_acid_changes)){
                                                   if(substr(data[i,j],p,p) == amino_acid_changes[p,1]){
                                                     substr(data[i,j],p,p) = amino_acid_changes[p,3]
                                                   }else if(substr(data[i,j],p,p) == amino_acid_changes[p,2]){
                                                     substr(data[i,j],p,p) = amino_acid_changes[p,4]
                                                   }
                                                   if(substr(data[i,p+2],j,j) == amino_acid_changes[p,1]){
                                                     substr(data[i,p+2],j,j) = amino_acid_changes[p,3]
                                                   }else{
                                                     substr(data[i,p+2],j,j) = amino_acid_changes[p,4]
                                                   }
                                                 }
                                               }
                                             }
                                             
                                             header = unique(c(names(table(data$comb1)),names(table(data$comb2))))
                                             allotype_individual_number = data.frame()
                                             for(i in 1:dim(data)[1]){
                                               allotype_individual_number[i,1:length(header)] = 0
                                               for(j in 1:2){
                                                 find_index = which(header == data[i,j])
                                                 if(allotype_individual_number[i,find_index] == 1){
                                                   allotype_individual_number[i,find_index] = 2
                                                 }else{
                                                   allotype_individual_number[i,find_index] = 1
                                                 }
                                               }
                                             }
                                             names(allotype_individual_number) = header
                                             rownames(allotype_individual_number) = rownames(data)
                                             return(list(allotype_combination = data,
                                                         allotype_individual_number = allotype_individual_number))
                                           }
                                         )
)
