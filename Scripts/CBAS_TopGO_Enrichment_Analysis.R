library("topGO")

library("BiocParallel")
register(MulticoreParam(4))



###
#
#using the background and degs obtained with script CBAS_DeSeq2_Analysis and the GO annotations for the transcripts in the universe use topGO to assess enrichment
#
####

#function GO analysis
alldegs2Function_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/all_degs_p001_2LFC.functions.gos")
allback2Function_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/Back_degs_p001_2LFC.functions.gos")

universe_for_function_topgo<-c(alldegs2Function_GO, allback2Function_GO)

#get names of genes in different groups
genes_of_interest_for_function<-names(alldegs2Function_GO)
genes_in_background_for_function<-names(allback2Function_GO)
genes_in_universe_for_function<-names(universe_for_function_topgo)

#get the list of genes of interest in the universe
genes_of_interest_list_for_function<-factor(as.integer(genes_in_universe_for_function %in% genes_of_interest_for_function))
names(genes_of_interest_list_for_function)<-genes_in_universe_for_function

#create a topGO object for Molecular Function
functions_go_data<-new("topGOdata", ontology="MF", allGenes=genes_of_interest_list_for_function, annot=annFUN.gene2GO, gene2GO=universe_for_function_topgo)

#create a test object
resultsFisher_functions<-runTest(functions_go_data, algorithm = "classic", statistic="fisher")

#generate table of results
GenTable(functions_go_data, classic=resultsFisher_functions, topNodes=50)
#generate graph of results
showSigOfNodes(functions_go_data, score(resultsFisher_functions), useInfo = "all")

#process GO analysis
alldegs2Process_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/all_degs_p001_2LFC.processes.gos")
allback2Process_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/Back_degs_p001_2LFC.processes.gos")

universe_for_process_topgo<-c(alldegs2Process_GO, allback2Process_GO)

#get names of genes in different groups
genes_of_interest_for_process<-names(alldegs2Process_GO)
genes_in_background_for_process<-names(allback2Process_GO)
genes_in_universe_for_process<-names(universe_for_process_topgo)

#get the list of genes of interest in the universe
genes_of_interest_list_for_process<-factor(as.integer(genes_in_universe_for_process %in% genes_of_interest_for_process))
names(genes_of_interest_list_for_process)<-genes_in_universe_for_process

#create a topGO object for Molecular Function
process_go_data<-new("topGOdata", ontology="BP", allGenes=genes_of_interest_list_for_process, annot=annFUN.gene2GO, gene2GO=universe_for_process_topgo)

#create a test object
resultsFisher_process<-runTest(process_go_data, algorithm = "classic", statistic="fisher")

#generate table of results
GenTable(process_go_data, classic=resultsFisher_process, topNodes=50)
#generate graph of results
showSigOfNodes(process_go_data, score(resultsFisher_process), useInfo = "all")


#component GO analysis
alldegs2Component_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/all_degs_p001_2LFC.components.gos")
allback2Component_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/Back_degs_p001_2LFC.components.gos")

universe_for_component_topgo<-c(alldegs2Component_GO, allback2Component_GO)

#get names of genes in different groups
genes_of_interest_for_component<-names(alldegs2Component_GO)
genes_in_background_for_component<-names(allback2Component_GO)
genes_in_universe_for_component<-names(universe_for_component_topgo)

#get the list of genes of interest in the universe
genes_of_interest_list_for_component<-factor(as.integer(genes_in_universe_for_component %in% genes_of_interest_for_component))
names(genes_of_interest_list_for_component)<-genes_in_universe_for_component

#create a topGO object for Molecular Function
component_go_data<-new("topGOdata", ontology="CC", allGenes=genes_of_interest_list_for_component, annot=annFUN.gene2GO, gene2GO=universe_for_component_topgo)

#create a test object
resultsFisher_component<-runTest(component_go_data, algorithm = "classic", statistic="fisher")

#generate table of results
GenTable(component_go_data, classic=resultsFisher_component, topNodes=50)
#generate graph of results
showSigOfNodes(component_go_data, score(resultsFisher_component), useInfo = "all")


##################################
#overexpressed genes only


#function GO analysis
overdegs2Function_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/over_degs_p001_2LFC.functions.gos")
overback2Function_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/Back_over_degs_p001_2LFC.functions.gos")

over_universe_for_function_topgo<-c(overdegs2Function_GO, overback2Function_GO)

#get names of genes in different groups
over_genes_of_interest_for_function<-names(overdegs2Function_GO)
over_genes_in_background_for_function<-names(overback2Function_GO)
over_genes_in_universe_for_function<-names(over_universe_for_function_topgo)

#get the list of genes of interest in the universe
over_genes_of_interest_list_for_function<-factor(as.integer(over_genes_in_universe_for_function %in% over_genes_of_interest_for_function))
names(over_genes_of_interest_list_for_function)<-over_genes_in_universe_for_function

#create a topGO object for Molecular Function
over_functions_go_data<-new("topGOdata", ontology="MF", allGenes=over_genes_of_interest_list_for_function, annot=annFUN.gene2GO, gene2GO=over_universe_for_function_topgo)

#create a test object
over_resultsFisher_functions<-runTest(over_functions_go_data, algorithm = "classic", statistic="fisher")

#generate table of results
GenTable(over_functions_go_data, classic=over_resultsFisher_functions, topNodes=50)
#generate graph of results
showSigOfNodes(over_functions_go_data, score(over_resultsFisher_functions), useInfo = "all")


#process GO analysis
overdegs2Process_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/over_degs_p001_2LFC.processes.gos")
overback2Process_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/Back_over_degs_p001_2LFC.processes.gos")

over_universe_for_process_topgo<-c(overdegs2Process_GO, overback2Process_GO)

#get names of genes in different groups
over_genes_of_interest_for_process<-names(overdegs2Process_GO)
over_genes_in_background_for_process<-names(overback2Process_GO)
over_genes_in_universe_for_process<-names(over_universe_for_process_topgo)

#get the list of genes of interest in the universe
over_genes_of_interest_list_for_process<-factor(as.integer(over_genes_in_universe_for_process %in% over_genes_of_interest_for_process))
names(over_genes_of_interest_list_for_process)<-over_genes_in_universe_for_process

#create a topGO object for Molecular Function
over_process_go_data<-new("topGOdata", ontology="BP", allGenes=over_genes_of_interest_list_for_process, annot=annFUN.gene2GO, gene2GO=over_universe_for_process_topgo)

#create a test object
over_resultsFisher_process<-runTest(over_process_go_data, algorithm = "classic", statistic="fisher")

#generate table of results
GenTable(over_process_go_data, classic=over_resultsFisher_process, topNodes=50)
#generate graph of results
showSigOfNodes(over_process_go_data, score(over_resultsFisher_process), useInfo = "all")


#component GO analysis
overdegs2Component_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/over_degs_p001_2LFC.components.gos")
overback2Component_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/Back_over_degs_p001_2LFC.components.gos")

over_universe_for_component_topgo<-c(overdegs2Component_GO, overback2Component_GO)

#get names of genes in different groups
over_genes_of_interest_for_component<-names(overdegs2Component_GO)
over_genes_in_background_for_component<-names(overback2Component_GO)
over_genes_in_universe_for_component<-names(over_universe_for_component_topgo)

#get the list of genes of interest in the universe
over_genes_of_interest_list_for_component<-factor(as.integer(over_genes_in_universe_for_component %in% over_genes_of_interest_for_component))
names(over_genes_of_interest_list_for_component)<-over_genes_in_universe_for_component

#create a topGO object for Molecular Function
over_component_go_data<-new("topGOdata", ontology="CC", allGenes=over_genes_of_interest_list_for_component, annot=annFUN.gene2GO, gene2GO=over_universe_for_component_topgo)

#create a test object
over_resultsFisher_component<-runTest(over_component_go_data, algorithm = "classic", statistic="fisher")

#generate table of results
GenTable(over_component_go_data, classic=over_resultsFisher_component, topNodes=50)
#generate graph of results
showSigOfNodes(over_component_go_data, score(over_resultsFisher_component), useInfo = "all")

######
#underexpressed genes only

#function GO analysis
underdegs2Function_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/under_degs_p001_2LFC.functions.gos")
underback2Function_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/Back_under_degs_p001_2LFC.functions.gos")

under_universe_for_function_topgo<-c(underdegs2Function_GO, underback2Function_GO)

#get names of genes in different groups
under_genes_of_interest_for_function<-names(underdegs2Function_GO)
under_genes_in_background_for_function<-names(underback2Function_GO)
under_genes_in_universe_for_function<-names(under_universe_for_function_topgo)

#get the list of genes of interest in the universe
under_genes_of_interest_list_for_function<-factor(as.integer(under_genes_in_universe_for_function %in% under_genes_of_interest_for_function))
names(under_genes_of_interest_list_for_function)<-under_genes_in_universe_for_function

#create a topGO object for Molecular Function
under_functions_go_data<-new("topGOdata", ontology="MF", allGenes=under_genes_of_interest_list_for_function, annot=annFUN.gene2GO, gene2GO=under_universe_for_function_topgo)

#create a test object
under_resultsFisher_functions<-runTest(under_functions_go_data, algorithm = "classic", statistic="fisher")

#generate table of results
GenTable(under_functions_go_data, classic=under_resultsFisher_functions, topNodes=50)
g#enerate graph of results
showSigOfNodes(under_functions_go_data, score(under_resultsFisher_functions), useInfo = "all")


#process GO analysis
underdegs2Process_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/under_degs_p001_2LFC.processes.gos")
underback2Process_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/Back_under_degs_p001_2LFC.processes.gos")

under_universe_for_process_topgo<-c(underdegs2Process_GO, underback2Process_GO)

#get names of genes in different groups
under_genes_of_interest_for_process<-names(underdegs2Process_GO)
under_genes_in_background_for_process<-names(underback2Process_GO)
under_genes_in_universe_for_process<-names(under_universe_for_process_topgo)

#get the list of genes of interest in the universe
under_genes_of_interest_list_for_process<-factor(as.integer(under_genes_in_universe_for_process %in% under_genes_of_interest_for_process))
names(under_genes_of_interest_list_for_process)<-under_genes_in_universe_for_process

#create a topGO object for Molecular Function
under_process_go_data<-new("topGOdata", ontology="BP", allGenes=under_genes_of_interest_list_for_process, annot=annFUN.gene2GO, gene2GO=under_universe_for_process_topgo)

#create a test object
under_resultsFisher_process<-runTest(under_process_go_data, algorithm = "classic", statistic="fisher")

#generate table of results
GenTable(under_process_go_data, classic=under_resultsFisher_process, topNodes=50)
#generate graph of results
showSigOfNodes(under_process_go_data, score(under_resultsFisher_process), useInfo = "all")


#component GO analysis
underdegs2Component_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/under_degs_p001_2LFC.components.gos")
underback2Component_GO<-readMappings("~/Desktop/CBAS/DeSeq2/Results/Back_under_degs_p001_2LFC.components.gos")

under_universe_for_component_topgo<-c(underdegs2Component_GO, underback2Component_GO)

#get names of genes in different groups
under_genes_of_interest_for_component<-names(underdegs2Component_GO)
under_genes_in_background_for_component<-names(underback2Component_GO)
under_genes_in_universe_for_component<-names(under_universe_for_component_topgo)

#get the list of genes of interest in the universe
under_genes_of_interest_list_for_component<-factor(as.integer(under_genes_in_universe_for_component %in% under_genes_of_interest_for_component))
names(under_genes_of_interest_list_for_component)<-under_genes_in_universe_for_component

#create a topGO object for Molecular Function
under_component_go_data<-new("topGOdata", ontology="CC", allGenes=under_genes_of_interest_list_for_component, annot=annFUN.gene2GO, gene2GO=under_universe_for_component_topgo)

#create a test object
under_resultsFisher_component<-runTest(under_component_go_data, algorithm = "classic", statistic="fisher")

#generate table of results
GenTable(under_component_go_data, classic=under_resultsFisher_component, topNodes=50)
#generate graph of results
showSigOfNodes(under_component_go_data, score(under_resultsFisher_component), useInfo = "all")



