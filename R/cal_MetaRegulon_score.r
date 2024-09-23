#' Get the date of the state, using for the MetaRegulon analysis

#' @param object A seurat object
#' @param feature Feature used to cluster cells, usually use the colnames of the meta.data of the seurat object 
#' @param state State of the interested cell types pr metabolic state
#' @param output_path Path to save the metacell result
#' @param file_name File name of the metacell result
#  @export

cal_MetaRegulon=function(object,feature,state,interested_MM,MM_list,markers,lisa_file,sample_tech,output_path,file_name){
    tmp_dir <- file.path(output_path,file_name)
    if (!dir.exists(tmp_dir)) {
        dir.create(tmp_dir)
    }
    if (sample_tech=='scRNA'){
        if (length(object@assays)>1) {
        # if 'integrated' exists，set it as the default
            DefaultAssay(object) <- "integrated"
            data=GetAssayData(object)
            data_use=data[,rownames(object@meta.data[object[[feature]]==state,])]
        } else {
            data=GetAssayData(object)
            data_use=data[,rownames(object@meta.data[object[[feature]]==state,])]
        }
    }
    if (sample_tech=='bulk'){
        data_use=object
    }
    
    file_res=paste0(tmp_dir,'/',file_name,".rds") ## get the state data to do the following analysis
    if (!file.exists(file_res)){
        saveRDS(data_use,file_res)
    }
     ## get the genome-wide correlation
    gg_cor_res=paste0(output_path,file_name,"/",file_name,":gg_activity_cor.rds")
    if (!file.exists(gg_cor_res)){
        data=readRDS(file_res)
        data=t(as.matrix(data))
        result=data %>%
            apply(2,function(x){10*(2**x - 1)}) %>%
            apply(2,function(x){log2(mean(x) + 1)}) 
        expressed_genes <- names(result[order(result, decreasing = TRUE)][1:3000]) ## screen the expressed genes
        data_exp=data[,expressed_genes]
        data_mg=data[rownames(data_exp),colnames(data) %in% unique(unlist(MM_list))]
        cor=cor(data_exp,data_mg)
        get_cor_res=lapply(as.list(interested_MM),get_cor,cor=cor,MM_list=MM_list,class='gene')
        names(get_cor_res)=interested_MM
        saveRDS(get_cor_res,gg_cor_res)    
    }

     ## get the ligand activity
    ligand_target_matrix = readRDS('/fs/home/tangke/metabolism/MetroTiME/MetroTiME_1214/fibroblast/regulation/Ref/ligand_target_matrix.rds')
    lr_network = readRDS('/fs/home/tangke/metabolism/MetroTiME/MetroTiME_1214/fibroblast/regulation/Ref/lr_network.rds')
    
    lr_activity_res=paste0(output_path,file_name,"/",file_name,":lr_activity.rds")
    if (!file.exists(lr_activity_res)){
        marker_state=markers[markers$cluster==state,'gene']
        data=t(as.matrix(readRDS(file_res)))
        data_use=data[,intersect(colnames(data),marker_state)]
        expressed_genes <-colnames(data_use)
        
        ligands = lr_network$from %>% unique()
        expressed_ligands = intersect(ligands,unique(markers$gene))
        receptors = lr_network$to %>% unique()
        expressed_receptors = intersect(receptors,marker_state) 
        
        potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% .$from %>% unique()
        background_expressed_genes = expressed_genes %>% .[. %in% rownames(ligand_target_matrix)]
        
        expression_scaled = data_use %>% .[,background_expressed_genes] %>% scale_quantile()
        
        ligand_activities = predict_single_cell_ligand_activities(cell_ids = rownames(data_use), 
                                                            expression_scaled = expression_scaled, 
                                                            ligand_target_matrix = ligand_target_matrix, 
                                                            potential_ligands = potential_ligands)
        normalized_ligand_activities = as.data.frame(normalize_single_cell_ligand_activities(ligand_activities))
        rownames(normalized_ligand_activities)=normalized_ligand_activities[,1]
        normalized_ligand_activities=normalized_ligand_activities[,-1]
        saveRDS(normalized_ligand_activities,lr_activity_res)
    }

    ## get the ligand-receptor interaction
    lr_cor_res=paste0(output_path,file_name,"/",file_name,":lr_activity_cor.rds")
    if (!file.exists(lr_cor_res)){
        data=readRDS(file_res)
        data_lr=readRDS(lr_activity_res)
        data_mg=t(as.matrix(data[rownames(data) %in% unique(unlist(MM_list)),rownames(data_lr)]))
        cor=cor(data_lr,data_mg)
        
        get_cor_res=lapply(as.list(interested_MM),get_cor,cor=cor,MM_list=MM_list,class='ligand')
        names(get_cor_res)=interested_MM
        saveRDS(get_cor_res,lr_cor_res)
    }   

    ## get the tr activity
    tr_activity_res=paste0(output_path,file_name,"/",file_name,":tr_activity.rds")
    if (!file.exists(tr_activity_res)){
        data=readRDS(file_res)
        lisa_res=read.csv(paste0(lisa_file),sep='\t')
        tf_gene<-unique(lisa_res$factor)

        M_TR=cal_tf_score(tf_gene,data) ## tf activity
        cal_target_score_res=lapply(as.list(colnames(M_TR)),
                                                cal_target_score,
                                                lisa_res=lisa_res,
                                                data=data) ## tf target score
        M_Target=Reduce(function(x,y){cbind(x,y)},cal_target_score_res)
        lisa_res_use=lisa_res[!duplicated(lisa_res$factor),c('factor','summary_p_value')]
        lisa_res_use$score=(-log(lisa_res_use[,c('summary_p_value')]))
        rownames(lisa_res_use)=lisa_res_use[,'factor']
        lisa_res_use=lisa_res_use[colnames(M_TR),] ## 

        M_actifvity=(M_TR+M_Target) ## sum the tf target score and tf activity score
        M_actifvity_scale=as.data.frame(lapply(M_actifvity, function(x) (x - min(x)) / (max(x) - min(x)))) ##scale
        rownames(M_actifvity_scale)=rownames(M_actifvity)
        S=M_actifvity_scale*lisa_res_use$score  ## scripro
        saveRDS(S,tr_activity_res)
    }

    ## get the tr-target interaction
    tr_cor_res=paste0(output_path,file_name,"/",file_name,":tr_activity_cor.rds")
    if (!file.exists(tr_cor_res)){
        data=readRDS(file_res)
        data_tf=readRDS(tr_activity_res)
        data_t=t(as.matrix(data[,rownames(data_tf)]))
        result=data_t %>%
            apply(2,function(x){10*(2**x - 1)}) %>%
            apply(2,function(x){log2(mean(x) + 1)}) 
        expressed_genes <- names(result[order(result, decreasing = TRUE)][1:3000]) 
        data_tf=data_tf[,colnames(data_tf) %in% expressed_genes]
        data_mg=t(as.matrix(data[rownames(data) %in% unique(unlist(MM_list)),rownames(data_tf)]))
        cor=cor(data_tf,data_mg)
        get_cor_res=lapply(as.list(interested_MM),get_cor,cor=cor,MM_list=MM_list,class='tr')

        names(get_cor_res)=interested_MM
        saveRDS(get_cor_res,tr_cor_res)
    }

    MetaRegulon_res=paste0(output_path,file_name,"/",file_name,":recom.rds")
    if (!file.exists(MetaRegulon_res)){
        tr_cor=readRDS(paste0(output_path,file_name,"/",file_name,":tr_activity_cor.rds"))
        ligand_cor=readRDS(paste0(output_path,file_name,"/",file_name,":lr_activity_cor.rds"))
        gene_cor=readRDS(paste0(output_path,file_name,"/",file_name,":gg_activity_cor.rds"))
        tr_activity=readRDS(paste0(output_path,file_name,"/",file_name,":tr_activity.rds"))
        lr_activity=readRDS(paste0(output_path,file_name,"/",file_name,":lr_activity.rds"))

        data=t(as.matrix(readRDS(file_res)))

        repref_each_res=lapply(as.list(interested_MM),repref_each,tr_cor=tr_cor,tf_activity=tr_activity,lr_activity=lr_activity,ligand_cor=ligand_cor,gene_cor=gene_cor,MM_list=MM_list,data=data)
        names(repref_each_res)=interested_MM
        saveRDS(repref_each_res,MetaRegulon_res)
    }
}

repref_each=function(tr_cor,tf_activity,lr_activity,ligand_cor,gene_cor,pathway_index,MM_list,data){
    if (pathway_index %in% names(tr_cor) & pathway_index %in% names(ligand_cor) & pathway_index %in% names(gene_cor)){
        
        tr_cor_use=tr_cor[[pathway_index]]
        ligand_cor_use=ligand_cor[[pathway_index]]
        gene_cor_use=gene_cor[[pathway_index]]
        
        if (nrow(tr_cor_use)>0 & nrow(ligand_cor_use)>0 &nrow(gene_cor_use)>0){
            li=list(TR_Target_interaction=tr_cor_use,Ligand_Receptor_interaction=ligand_cor_use,Gene_Gene_interaction=gene_cor_use) #
            li_merge=Reduce(function(x,y){merge(x,y,by='MR',all=TRUE)},li)
            rownames(li_merge)=li_merge$MR

            a=li_merge[!is.na(li_merge$TR_Target_interaction>0),]
            a=rownames(a[order(a$TR_Target_interaction,decreasing=TRUE),])
            b=li_merge[!is.na(li_merge$Ligand_Receptor_interaction>0),]
            b=rownames(b[order(b$Ligand_Receptor_interaction,decreasing=TRUE),])
            c=li_merge[!is.na(li_merge$Gene_Gene_interaction>0),]
            c=rownames(c[order(c$Gene_Gene_interaction,decreasing=TRUE),])

            glist=list(a=a,b=b,c=c)
            ag=aggregateRanks(glist,exact=TRUE)
            ag=ag[rownames(li_merge),]
            li_merge$ag_score=ag$Score
            li_merge[is.na(li_merge)]<-0
            li_merge<-li_merge[,-1]
            res_pr=psel(li_merge, high(TR_Target_interaction) *high(Ligand_Receptor_interaction)*high(Gene_Gene_interaction), top_level=30)
            res_pr$Final_score=res_pr$ag_score*res_pr$.level
            res_pr$gene=rownames(res_pr)
            res_pr <- res_pr %>%
                group_by(.level) %>%
                arrange(ag_score,.by_group = TRUE) %>%
                as.data.frame()
            res_pr$rank=1:nrow(res_pr)
            rownames(res_pr)=res_pr$gene

            ## the source of the regulation
            ## identify the regulation source
            res_pr$resource <- ifelse(res_pr$TR_Target_interaction > 0 | res_pr$Gene_Gene_interaction > 0,
                          'intrinsic', 'extrinsic')

            if (length(MM_list[[pathway_index]])<30){
                res_pr$direction='' ## setting the direction of these factor and downstream metabolic reaction
                times=ifelse(nrow(data)>300,5,1)
                factor_tf=res_pr[res_pr$TR_Target_interaction>0,'gene']
                factor_gg=res_pr[res_pr$TR_Target_interaction==0 & res_pr$Ligand_Receptor_interaction==0,'gene']
                factor_lr_paracrine=res_pr[res_pr$TR_Target_interaction==0 & res_pr$Gene_Gene_interaction==0,'gene']
                factor_lr_autocrine=res_pr[res_pr$TR_Target_interaction==0 & !res_pr$Gene_Gene_interaction==0 & !res_pr$Ligand_Receptor_interaction==0,'gene']

                ## for ligand, if the ligand is paracrine, then it is regulator, if the ligand is autocrine, MetroSCREEN will identify it's casauation
                res_pr$direction=ifelse(rownames(res_pr) %in% factor_lr_paracrine,'regulator',res_pr$direction)

                metabolic_gene=unique(MM_list[[pathway_index]])

                if (length(factor_tf)>0){
                    factor_tf_res=do.call(rbind,
                        lapply(as.list(factor_tf),get_causal,data=data,tr_activity=tr_activity,lr_activity=lr_activity,MM_list=MM_list,pathway_index=pathway_index,times=times,factor_class='TF'))
                    rownames(factor_tf_res)=factor_tf_res$factor
                    res_pr$direction=ifelse(rownames(res_pr) %in% factor_tf,factor_tf_res$direction,res_pr$direction)
                }

                if (length(factor_gg)>0){
                    factor_gg_res=do.call(rbind,
                        lapply(as.list(factor_gg),get_causal,data=data,tr_activity=tr_activity,lr_activity=lr_activity,MM_list=MM_list,pathway_index=pathway_index,times=times,factor_class='GG'))
                    rownames(factor_tf_res)=factor_tf_res$factor
                    res_pr$direction=ifelse(rownames(res_pr) %in% factor_tf,factor_tf_res$direction,res_pr$direction)
                }

                if (length(factor_lr_autocrine)>0){
                    factor_lr_autocrine_res=do.call(rbind,
                        lapply(as.list(factor_lr_autocrine),get_causal,data=data,tr_activity=tr_activity,lr_activity=lr_activity,MM_list=MM_list,pathway_index=pathway_index,times=times,factor_class='Ligand'))
                    rownames(factor_tf_res)=factor_tf_res$factor
                    res_pr$direction=ifelse(rownames(res_pr) %in% factor_tf,factor_tf_res$direction,res_pr$direction)
                }
                return(res_pr)
            }
        }
    }
}

get_causal <- function(data,tr_activity,lr_activity,factor,MM_list,pathway_index,times,factor_class) {
    if (factor_class=='TF'){ ## TF and gene use different matrix to get the max correlation
        df_tf=as.data.frame(tr_activity[,factor])
        rownames(df_tf)=rownames(tr_activity)
        colnames(df_tf)=factor

        df_gene=as.data.frame(data[,colnames(data) %in% MM_list[[pathway_index]]])
        unique_counts <- sapply(as.data.frame(df_gene), function(col) length(unique(col)))
        df_gene_filtered <- df_gene[, unique_counts >= 10]
        cor=cor(df_tf,df_gene_filtered)[1,]

        cor=cor[order(cor,decreasing = TRUE)]
        gene_max=names(cor[1])
    }
    if (factor_class=='GG'){
        df_gg=as.data.frame(data[,factor])
        rownames(df_gg)=rownames(data)
        colnames(df_gg)=factor

        df_gene=as.data.frame(data[,colnames(data) %in% MM_list[[pathway_index]]])
        unique_counts <- sapply(as.data.frame(df_gene), function(col) length(unique(col)))
        df_gene_filtered <- df_gene[, unique_counts >= 10]
        cor=cor(df_gg,df_gene_filtered)[1,]

        cor=cor[order(cor,decreasing = TRUE)]
        gene_max=names(cor[1]) ## get the highst correlation of the MM
    }
    if (factor_class=='Ligand'){
        df_lr=as.data.frame(lr_activity[,factor])
        rownames(df_lr)=rownames(lr_activity)
        colnames(df_lr)=factor

        df_gene=as.data.frame(data[,colnames(data) %in% MM_list[[pathway_index]]])
        unique_counts <- sapply(as.data.frame(df_gene), function(col) length(unique(col)))
        df_gene_filtered <- df_gene[, unique_counts >= 10]
        cor=cor(df_lr,df_gene_filtered)[1,]

        cor=cor[order(cor,decreasing = TRUE)]
        gene_max=names(cor[1])
    }
    
    df_use=data[,colnames(data) %in% c(factor,unique(MM_list[[pathway_index]]))]
    results <- list()
    for (i in 1:times) { ## identify the causation
            rows <- if (times == 1) {
                1:nrow(df_use)
            } else {
                sample(1:nrow(df_use), 300)
            }

            df_final=df_use[rows,]
            unique_counts <- sapply(as.data.frame(df_final), function(col) length(unique(col)))
            df_final_filtered <- df_final[, unique_counts >= 10]

            # build causal network
            kpc <- kpc(suffStat = list(data = df_final_filtered, ic.method = "dcc.gamma"),
                       indepTest = kernelCItest,
                       alpha = 0.1,
                       labels = colnames(df_final_filtered),
                       u2pd = "relaxed",
                       skel.method = "stable.fast",
                       verbose = TRUE)

            g <- graph_from_graphnel(kpc@graph, name = TRUE, weight = TRUE, unlist.attrs = TRUE)
            edges_df <- as.data.frame(get.edgelist(g), stringsAsFactors=FALSE)
            edges_df_use=edges_df[edges_df$V1==factor | edges_df$V2==factor,]
            colnames(edges_df_use)=c('from','to')

            if (nrow(edges_df_use)>0){
                res <- if (any(edges_df_use$from == factor & edges_df_use$to == gene_max)) {
                    1
                } else {
                    0 
                }
                results[[i]] <- res
            }
            if (nrow(edges_df_use)==0){
                results[[i]] <- 0
            }
    }
    ## when there are at least one in three or two in five times repeat the results, we think it is rubost
    result=sum(unlist(results))/length(results) 
    df=data.frame(factor=factor,direction=ifelse(result>0.3,'regulator','effector'))
    return(df)
}


cal_tf_score=function(tf,data){
    tf_orig=data[rownames(data) %in% tf,] ## lisa res 
    tf_orig_zscore=apply(tf_orig,1,function(x){(x-mean(x))/sd(x)}) ## normalize
    cols_with_nan <- apply(tf_orig_zscore, 2, function(x) any(is.nan(x))) ## remove na
    tf_orig_zscore_filtered <- tf_orig_zscore[, !cols_with_nan]
    tf_orig_zscore_filtered[tf_orig_zscore_filtered>4]=4 
    tf_orig_zscore_filtered[tf_orig_zscore_filtered<(-4)]=(-4)            
    return(tf_orig_zscore_filtered)
}
                           
cal_target_score=function(lisa_res,factor,data){
    sample_id=lisa_res[lisa_res$factor==factor,'sample_id'][1] 
    if (nchar(sample_id)>2){
        dir_id=substring(sample_id,1,3)
    }else{
        dir_id='000'
    }
    tf_target=read.table(paste0('/fs/data/cistrome/RP/',dir_id,'/',sample_id,'_gene_score.txt'))
    tf_target_use=tf_target[!duplicated(tf_target$V7),]
    target=tf_target_use[tf_target_use$V5>2,'V7']
    if (length(target)<500){
        tf_target_use_order=tf_target_use[order(tf_target_use$V5,decreasing = TRUE) & tf_target_use$V5>0,]
        if (nrow(tf_target_use_order>=500)){
            target=tf_target_use_order[1:500,'V7']
        }else{
            target=tf_target_use_order[,'V7']
        }
    }
    
    target_orig=data[rownames(data) %in% unique(target),] 
    target_orig_zscore=apply(target_orig,1,function(x){(x-mean(x))/sd(x)})
    cols_with_nan <- apply(target_orig_zscore, 2, function(x) any(is.nan(x)))
    target_orig_zscore_filtered <- target_orig_zscore[, !cols_with_nan] 
    target_orig_zscore_filtered_mean=apply(target_orig_zscore_filtered,1,mean)                       
    target_orig_mean=as.data.frame(target_orig_zscore_filtered_mean)
                           
    target_orig_mean[target_orig_mean>4]=4 
    target_orig_mean[target_orig_mean<(-4)]=(-4)      
                           
    colnames(target_orig_mean)=factor
    return(target_orig_mean)
}

get_cor=function(cor,MM_list,MM_index,class){
    cor=cor[,colnames(cor) %in% MM_list[[MM_index]]]
    if(is.matrix(cor)){
        cor <- cor[, colSums(is.na(cor)) != nrow(cor)]
        if(is.matrix(cor)){
            cor_max=as.data.frame(apply(na.omit(cor),1,max))
            cor_max$MR=rownames(cor_max)
        }
        else{
            cor_max=as.data.frame(na.omit(cor))
            cor_max$MR=rownames(cor_max)
        }
    } else{
        cor_max=as.data.frame(na.omit(cor))
        cor_max$MR=rownames(cor_max)
    }
    if (class=='gene'){
        colnames(cor_max)[1]='Gene_Gene_interaction'
        cor_max=cor_max[!rownames(cor_max) %in% MM_list[[MM_index]],]
        cor_max=cor_max[cor_max$Gene_Gene_interaction>0.3,]
        return(cor_max)
    }
    if (class=='ligand'){
        colnames(cor_max)[1]='Ligand_Receptor_interaction'
        cor_max=cor_max[!rownames(cor_max) %in% MM_list[[MM_index]],]
        cor_max=cor_max[cor_max$Ligand_Receptor_interaction>0,]
        return(cor_max)
    }
    if (class=='tr'){
        colnames(cor_max)[1]='TR_Target_interaction'
        cor_max=cor_max[!rownames(cor_max) %in% MM_list[[MM_index]],]
        cor_max=cor_max[cor_max$TR_Target_interaction>0,]
        return(cor_max)
    }
    
}


