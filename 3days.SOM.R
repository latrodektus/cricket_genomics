library(edgeR)
library(genefilter)

data = read.table("TF.count.csv", header=T, row.names=1, com='')

rnaseqMatrix = round(data)
conditions=factor(c(rep("LV_F_3",3),rep("LV_F_13",3),rep("LV_M_3",3),rep("LV_M_13",3),rep("No_F_3",3),rep("No_F_13",3),rep("No_M_3",3),rep("No_M_13",3)))
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)

Ave_LV_F_3=rowMeans(exp_study$counts[,c(1,2,3)])
Ave_LV_F_13=rowMeans(exp_study$counts[,c(4,5,6)])
Ave_LV_M_3=rowMeans(exp_study$counts[,c(7,8,9)])
Ave_LV_M_13=rowMeans(exp_study$counts[,c(10,11,12)])
Ave_No_F_3=rowMeans(exp_study$counts[,c(13,14,15)])
Ave_No_F_13=rowMeans(exp_study$counts[,c(16,17,18)])
Ave_No_M_3=rowMeans(exp_study$counts[,c(19,20,21)])
Ave_No_M_13=rowMeans(exp_study$counts[,c(22,23,24)])
AveCounts <- data.frame(Ave_LV_F_3,Ave_LV_M_3,Ave_No_F_3,Ave_No_M_3)

d_Ave <- DGEList(counts=AveCounts)
d_Ave <- calcNormFactors(d_Ave)

d_Ave <- estimateGLMCommonDisp(d_Ave, verbose=TRUE)
Ave.LogCPM <- predFC(d_Ave)
logCPM.scaled <- genescale(Ave.LogCPM, axis=1, method="Z")
genes.id=rownames(data)

#run SOMs
library(kohonen)
library(lattice)
convert_to_format_for_lattice <- function (  som, input_matrix, mic_id, ypd_id, ncol, nrow) 
{

        data.dimension <- dim( input_matrix)
        data.total_length <- data.dimension[1] * data.dimension[2]
        data.num_mic_id <- data.dimension[1]

        som_data_pivot <- data.frame(  id_a= rep(0, data.total_length),  id_b= rep(0, data.total_length), 
                                       x = rep(0, data.total_length), y= rep(0, data.total_length), cell=rep(0, data.total_length) )

        row_count <- 1
        for ( i in 1:(ncol*nrow) )
        {    
                temp_input_matrix_data <- input_matrix[ som$unit.classif == i, ]
                temp_mic_id     <- mic_id[ som$unit.classif == i] 
                temp_ypd_id     <- ypd_id[ som$unit.classif == i] 
                temp_num_genes  <- length(which(som$unit.classif == i ) )

                if ( temp_num_genes  > 0 )
                {     
                        for ( j in 1:temp_num_genes )
                        {
                                my_k_length <- 0
                                if ( temp_num_genes > 1 )
                                {
                                      my_k_length <- length(temp_input_matrix_data[1,]  )
                                } else 
                                {
                                      my_k_length <- length(temp_input_matrix_data )
                                }       

                                for ( k in 1:my_k_length )
                                {
                                      som_data_pivot[row_count,"id_a"] <-  as.character(temp_mic_id[j])
                                      som_data_pivot[row_count,"id_b"] <-  as.character(temp_ypd_id[j])
                                      som_data_pivot[row_count,"x"] <-  k 

                                      som_data_pivot[row_count,"cell"] <-  i 
                                
                                      if ( temp_num_genes > 1 )
                                      {
                                         som_data_pivot[row_count,"y"] <-  temp_input_matrix_data[j, k] 
                                      } else if( temp_num_genes == 1 )
                                      {
                                         som_data_pivot[row_count,"y"] <-  temp_input_matrix_data[k] 
                                      }  
                                                                
  
                                      row_count <- row_count + 1
                                }
                        }
                } else
                {

                    # treat empty cells                  
                    som_data_pivot <-  som_data_pivot 
    
                   for ( num_time_points in 1:data.dimension[2] )
                   {     
                     som_data_pivot <- rbind( som_data_pivot, 
                                                c(NA, NA, num_time_points, NA, i ))

                   }
                }

        }

        return ( som_data_pivot) 
}

rectangular_5_by_5 <- somgrid(xdim = 5, ydim = 5, topo = c("rectangular"))
### Need to set the seed so that every time you run the code, the results are the same
set.seed(7)

 all.som <- som(logCPM.scaled, grid=rectangular_5_by_5 , rlen = 100, toroidal = FALSE,   keep.data = TRUE)  

som_data_pivot <- convert_to_format_for_lattice(all.som, logCPM.scaled, genes.id, genes.id, ncol=5, nrow=5) 

### This variable gives a colour to each SOM cell. The order of the colours must match the 'data = som_data_pivot' variable.
colours <- factor (som_data_pivot[,'cell'] , labels=rainbow(5*5) )
scales_list=list(x=list(labels=c("","F","M","F","M"),at=seq(0,4,1),cex=0.3))
### This provides the grouping to the 'panel' function. We can then use 'subscripts' to select the groups specific to each cell.
grouping <- factor(som_data_pivot[,'id_a'])
xyplot (y ~ x | cell, 
	data=som_data_pivot, 
	groups= grouping, 
	layout=c(5,5),
	strip=FALSE,
    scales=scales_list,
	xlab="Sample", 
	ylab="Expression Level",
	user.defined.color=as.character(colours) ,user.defined.groups= grouping,scales,
	panel = function(x, y,user.defined.color,user.defined.groups,...,subscripts) {
	my_fill <- user.defined.color[subscripts]
	my_group <- user.defined.groups[subscripts]
	panel.xyplot(x, y, type=c("l", "p"), col=my_fill, groups=my_group, subscripts=TRUE) 
	panel.text(1,1.5,subscripts)
	} )
cell_for_each_gene_id <- unique(som_data_pivot[,c("id_a", "id_b", "cell")])
write.table(cell_for_each_gene_id, "3days_cell_for_each_geneID_5x5.txt", sep="\t")