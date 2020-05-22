################################################################################
# Offspring simulator; see README
################################################################################
# LICENCE

# MIT License
# 
# Copyright (c) University of Helsinki
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

################################################################################
# LIBRARIES
library(data.table) # Make sure this is newest version
set.seed(10)

################################################################################

################################################################################
# FUNCTIONS
################################################################################

splitAt <- function(x, pos){
  unname(split(x, cumsum(seq_along(x) %in% pos) ) )
} 

haplo_combine <- function( recomb_indeces, haplo1, haplo2 ){
  
  # Split haplotypes into blocks based on recombination spots (recomb_indeces)
  blocks1 <- splitAt( haplo1, recomb_indeces)
  blocks2 <- splitAt( haplo2, recomb_indeces)
  # Merge haplotypes to list so that it can be iteratively accessed
  blocks_list <- list( blocks1, blocks2 )
  
  # Create block_order list to keep track of which haplotype should be accessed 
  # next in block_list 
  block_order <- rep( c(1,2), 100 )[ 1:(length(recomb_indeces)+1) ]
  
  # Create new haplotype by choosing blocks from block_list
  haplo_new <- c()
  for( block in 1:length(block_order) ){
    # We select from block list first correct haplotype 1/2, second correct block
    haplo_new <- append( haplo_new, blocks_list[[ block_order[ block ] ]][[ block ]] )
  }
  
  return( haplo_new )  
}


sample_crossover <- function( recombination ){
  rec_events <- rbinom( nrow(recombination), size = 1, prob = recombination$recom.rate )
  rec_events <- as.logical(rec_events)
  rec_indeces <- which( rec_events )
  print( paste( "Number of Cross-overs:", sum(rec_events) ) )
  print( paste( "Events at ", paste( rec_indeces, collapse = " ") ) )
  return( rec_indeces )
  
}


makeOffspring <- function( parent1, parent2, phase_sample, phase_file ){
  
  # FIND PARENTAL HAPLOTYPES FROM PHASE FILES
  phase_id_order <- c( "X", "X", "X", rep( phase_sample$V1, each = 2 ) )
  parent1_indeces <- which( phase_id_order == parent1 ) 
  parent2_indeces <- which( phase_id_order == parent2 ) 
  
  print( parent1_indeces )
  
  # Read haplotypes from file
  parent1_hap1 <- fread( phase_file, skip = parent1_indeces[1]-1 , nrows = 1, header = F )
  parent1_hap2 <- fread( phase_file, skip = parent1_indeces[2]-1 , nrows = 1, header = F )
  
  parent2_hap1 <- fread( phase_file, skip = parent2_indeces[1]-1 , nrows = 1, header = F )
  parent2_hap2 <- fread( phase_file, skip = parent2_indeces[2]-1 , nrows = 1, header = F )
  
  # Convert haplotypes to numeric vector
  parent1_hap1 <- as.numeric( unlist(strsplit( as.character( parent1_hap1[1,1]), "" ) ) )
  parent1_hap2 <- as.numeric( unlist(strsplit( as.character( parent1_hap2[1,1]), "" ) ) )
  
  parent2_hap1 <- as.numeric( unlist(strsplit( as.character( parent2_hap1[1,1]), "" ) ) )
  parent2_hap2 <- as.numeric( unlist(strsplit( as.character( parent2_hap2[1,1]), "" ) ) )


  # MODEL CROSSING OVER N TIMES
  # Sample cross-over events
  rec_indeces1 <- sample_crossover(recombination)
  rec_indeces2 <- sample_crossover(recombination)
  
  # Perform cross-over by combining haplotypes
  if( length(rec_indeces1) == 0 ){
    # No cross-over happens
    parent1_hap1_new <- parent1_hap1
    parent1_hap2_new <- parent1_hap2 
  } else {
    # Cross-over between haplotypes
    parent1_hap1_new <- haplo_combine( rec_indeces1, parent1_hap1, parent1_hap2 )
    parent1_hap2_new <- haplo_combine( rec_indeces1, parent1_hap2, parent1_hap1 )
  }
  
  if( length(rec_indeces2) == 0 ){
    # No cross-over happens
    parent2_hap1_new <- parent2_hap1
    parent2_hap2_new <- parent2_hap2
  } else {
    # Cross-over between haplotypes
    parent2_hap1_new <- haplo_combine( rec_indeces2, parent2_hap1, parent2_hap2 )
    parent2_hap2_new <- haplo_combine( rec_indeces2, parent2_hap2, parent2_hap1 )
  }
  
  # SAMPLE OFFSPRING 
  off_hap1 <- ifelse( sample(1:2, 1) == 1 , list(parent1_hap1_new), list(parent1_hap2_new) )
  off_hap2 <- ifelse( sample(1:2, 1) == 1 , list(parent2_hap1_new), list(parent2_hap2_new) )
  
  return( list( off_hap1, off_hap2 ) )
}


makeGeneration <- function( parents, input_phase_file, output_file_phase, 
                            input_label_file, output_label_file,
                            recombination ){
 
  phase_sample <- read.table(input_label_file, as.is = T)
  
  # Make offspring haplotypes and saves the to phase file
  for( p in 1:nrow(parents) ){
    
    print(p)
    print(parents)
    off_haps <- makeOffspring( parents[ p, 1] , parents[ p, 2 ], 
                               phase_sample, input_phase_file )
  
    # Write offspring into file
    fwrite( list( paste(off_haps[[1]][[1]], collapse = "") ), 
            output_file_phase, append = T )
    fwrite( list( paste(off_haps[[2]][[1]], collapse = "") ), 
            output_file_phase, append = T )
  }

  
  # Make label_infile for offsprings
  sample_list <- paste0( "F", generation+1, "_", 1:nrow(parents) )
  sample_df <- as.data.frame(matrix( NA, nrow = nrow(parents), ncol = 3 ) )
  sample_df$V1 <- sample_list
  sample_df$V2 <- sample_list
  sample_df$V3 <- 1
  write.table( sample_df, output_label_file, append = T,
               quote = F, row.names = F, col.names = F)
  
  return( sample_df )
}


makeParentalFile <- function( sample_file, output_file ){
  
  if( nrow(sample_file) %% 2 == 1 ){
    sample_file <- sample_file[ 1:(nrow(sample_file)-1), ]
  }
  
  new_parents <- as.data.frame( matrix( NA, nrow = floor(nrow(sample_file) / 2), ncol = 2 ) )
  new_parents[ , 1 ] <- sample_file[ c(T,F) , 1 ]
  new_parents[ , 2 ] <- sample_file[ c(F,T) , 1 ]
  
  write.table( new_parents, output_file ,
               quote = F, row.names = F, col.names = F )
  
  print( paste( "Parents ready." ) )
}


################################################################################
################################################################################
# INPUT TEST FILES
file_name <- "New"
chr <- 1
recomb_file <- "Testi_recombination.txt"
input_phase_file <- "Testi.phase"
input_label_file <- "Test_label_infile.txt"
input_parent_file <- "Testi.F0.txt"

output_file_phase <- paste0( file_name, "_chr", chr, ".phase" )
output_label_file <- paste0( file_name, "_chr", chr, "_ladel_infile.txt" )
recombination <- read.table(recomb_file , h = T, as.is = T) # Note: recombination 
# is defined as global variable (data.frame) and not given as an argument for makeOffspring function.


generation <- 0

while( TRUE ){
  
  parents <- read.table(input_parent_file, as.is = T)
  
  if( generation == 0){
  
    print( paste( "Generation: 0") )    
  
    # Initialize output phase file, this needs to be updated at the end
    write( 9999, output_file_phase )
    write( nrow(recombination), output_file_phase, append = T )
    write( paste( c("P", recombination$position) , collapse = " " ), 
           output_file_phase, append = T )
    
    # Read first generation haplotypes from input file, others from output file
    # MakeGeneration generates new haplotypes, saves them to file and 
    # returns label_infile of new samples
    new_samples <- makeGeneration( parents, input_phase_file, output_file_phase, 
                                   input_label_file, output_label_file,
                                   recombination )
    
  } else {
    
    print( paste( "Generation:", generation) )  
    
    # Read next gen haplotypes from pre-existing output file
    # MakeGeneration generates new haplotypes, saves them to file and 
    # returns label_infile of new samples
    new_samples <- makeGeneration( parents, output_file_phase, output_file_phase, 
                                   output_label_file, output_label_file,
                                   recombination )
  }
  
    
    if( nrow(new_samples) >= 2 ){
      
      generation <- generation + 1
      
      output_parent_file <- paste0( file_name, "_chr", chr, "_parents_F", generation,".txt" )
      makeParentalFile( new_samples, output_parent_file )
      
      input_parent_file <- output_parent_file 
      
    } else {
      break
    }

}

################################################################################
################################################################################