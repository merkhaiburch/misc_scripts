# ---------------------------------------------------------------
# Author.... Merritt Burch
# Contact... mbb262@cornell.edu
# Date...... 2020-10-28 
# Updated... 2021-10-11
#
# Description 
#   - Script creates paired, or triplet ranges out of phg haplotypes
#   - to use when filtering GWAS results
# ---------------------------------------------------------------

# ---------------
## Functions ----
# ---------------

# Function to create new ranges by pairing old PHG ref ranges
pair_ref_ranges <- function(ranges) {
  paired_ranges <- c()
  for (i in seq_len(10)){
    
    # sort by chrom and start position
    ranges <- ranges %>% arrange(seqid, start) 
    
    # filer to single chrom
    single_chrom <- ranges %>% filter(seqid == i)
    print(dim(single_chrom))
    
    # if chromosome has an even number of ranges:
    # Pair sorted ranges together by taking start of first range and end of second range to create a new range
    if (nrow(single_chrom) %% 2 == 0){
      # Get starting points of new ranges, should be every even column
      tmp_start <- single_chrom$start[seq(1, nrow(single_chrom), 2)]
      
      # get end points for new ranges, should be every odd column
      tmp_end <- single_chrom$end[seq(2, nrow(single_chrom), 2)]
      
      # Combine
      tmp <- cbind(tmp_start, tmp_end) %>% data.table::as.data.table()
      
    } else {
      
      # If there is an odd number of columns, remove the last row (i.e. make the df have an even # of rows)
      # use the same methods as above, make the last/removed row it's own range (don't pair it with another range)
      # Get starting points of new ranges, should be every even column
      tmp_start <- single_chrom$start[seq(1, nrow(single_chrom)-1, 2)]
      
      # get end points for new ranges, should be every odd column
      tmp_end <- single_chrom$end[seq(2, nrow(single_chrom)-1, 2)]
      
      # Combine
      tmp <- cbind(tmp_start, tmp_end) %>% data.table::as.data.table()
      
      # rbind on last unpaired range
      lala <- single_chrom[nrow(single_chrom),2:3]
      tmp <- rbind(tmp, lala, use.names = FALSE)
    }
    # Add chromosome info
    tmp$seqid <- rep(i, nrow(tmp))
    
    # Add new ref range id
    tmp$rr_id <- paste0("paired_rr_", seq(1, nrow(tmp)))
    
    # Rearrage, rename columns
    colnames(tmp) <- c("start", "end", "seqid", "rr_id")
    tmp <- tmp %>% select("seqid", "start", "end", "rr_id")
    
    # Add to outside df
    paired_ranges <- rbind(paired_ranges, tmp)
  }
  return(paired_ranges)
}

# Use function
paired_ranges <- pair_ref_ranges(ranges)

# Count paired haplotype length
# temp <- paired_ranges$end-paired_ranges$start
# summary(temp)
# length(temp)
# 
# lala <- ranges$end - ranges$start
# summary(lala)
# length(lala)


# ----------------------------------
# Function to make triplet sliding 
#   window haplotyoes
# ----------------------------------

# Take all three window haplotypes, gather top 1% of hits
first_triplet_window <- function(ranges) {
  triplet_ranges <- c()
  for (i in seq_len(10)){
    
    # sort by chrom and start position
    ranges <- ranges %>% arrange(seqid, start)
    
    # filer to single chrom
    single_chrom <- ranges %>% filter(seqid == i)
    print(dim(single_chrom))
    
    # if chromosome has an odd number of ranges:
    # Pair sorted ranges together by taking start of first range and end of second range to create a new range
    if (nrow(single_chrom) %% 3 == 0){
      
      print(message("There are an odd number of rows"))
      
      # Get starting points of new ranges, should be every even column
      tmp_start <- single_chrom$start[seq(1, nrow(single_chrom), 3)]
      
      # get end points for new ranges, should be every odd column
      tmp_end <- single_chrom$end[seq(3, nrow(single_chrom), 3)]
      
      # Combine
      tmp <- cbind(tmp_start, tmp_end) %>% data.table::as.data.table()
      
    } else if (nrow(single_chrom) %% 3 == 1) {
      
      print(message("There are an even number of rows and \n and there is one additional row"))
      
      # If there is an even number of columns, remove the last 1 or 2 row (i.e. make the df have an odd # of rows)
      # use the same methods as above, make the last/removed row it's own range (don't pair it with another range)
      # Get starting points of new ranges, should be every even column
      tmp_start <- single_chrom$start[seq(1, nrow(single_chrom)-1, 3)]
      
      # get end points for new ranges, should be every odd column
      tmp_end <- single_chrom$end[seq(3, nrow(single_chrom)-1, 3)]
      
      # Combine
      tmp <- cbind(tmp_start, tmp_end) %>% data.table::as.data.table()
      
      # rbind on last unpaired range
      lala <- single_chrom[nrow(single_chrom),2:3]
      
      # Combine with triplet ranges
      tmp <- rbind(tmp, lala, use.names = FALSE)
      
    } else if (nrow(single_chrom) %% 3 == 2){
      
      print(message("There are an even number of rows and \n and there are two additional rows"))
      
      # If there is an even number of columns, remove the last 1 or 2 row (i.e. make the df have an odd # of rows)
      # use the same methods as above, make the last/removed row it's own range (don't pair it with another range)
      # Get starting points of new ranges, should be every even column
      tmp_start <- single_chrom$start[seq(1, nrow(single_chrom)-2, 3)]
      
      # get end points for new ranges, should be every odd column
      tmp_end <- single_chrom$end[seq(3, nrow(single_chrom)-2, 3)]
      
      # Combine
      tmp <- cbind(tmp_start, tmp_end) %>% data.table::as.data.table()
      
      # If there are two ranges remaining, turn into a single range
      lala <- data.frame(single_chrom$start[nrow(single_chrom)-1], single_chrom$end[nrow(single_chrom)]) %>% data.table::as.data.table()
      
      # Change colnames
      colnames(lala) <- c("start", "end")
      
      # Combine with triplet ranges
      tmp <- rbind(tmp, lala, use.names = FALSE)
    }
    
    # Add back chromosome info
    tmp$seqid <- rep(i, nrow(tmp))
    
    # Rearrage, rename columns
    colnames(tmp) <- c("start", "end", "seqid")
    tmp <- tmp %>% select("seqid", "start", "end")
    
    # Add to outside df
    triplet_ranges <- rbind(triplet_ranges, tmp)
  }
  
  return(triplet_ranges)
}
first_triplet_haps <- first_triplet_window(ranges)
first_triplet_haps$rr_id <- paste0("first_triplet_rr_", seq(1, nrow(first_triplet_haps)))


# Move one haplotype over, get three haplotypes, turn into a single haplotype
second_triplet_window <- function(ranges) {
  triplet_ranges <- c()
  for (i in seq_len(10)){
    
    # sort by chrom and start position
    ranges <- ranges %>% arrange(seqid, start) 
    
    # filer to single chrom
    single_chrom <- ranges %>% filter(seqid == i)
    print(dim(single_chrom))
    
    rows_working_with <- nrow(single_chrom)-1
    
    # if chromosome has an odd number of ranges:
    # Pair sorted ranges together by taking start of first range and end of second range to create a new range
    if (rows_working_with %% 3 == 0){
      
      print(message("There are an odd number of rows"))
      
      # Get starting points of new ranges, should be every even column
      tmp_start <- single_chrom$start[seq(2, nrow(single_chrom), 3)]
      
      # get end points for new ranges, should be every odd column
      tmp_end <- single_chrom$end[seq(4, nrow(single_chrom), 3)]
      
      # Combine
      tmp <- cbind(tmp_start, tmp_end) %>% data.table::as.data.table()
      
    } else if (rows_working_with %% 3 == 1) {
      
      print(message("There are an even number of rows and \n and there is one additional row"))
      
      # If there is an even number of columns, remove the last 1 or 2 row (i.e. make the df have an odd # of rows)
      # use the same methods as above, make the last/removed row it's own range (don't pair it with another range)
      # Get starting points of new ranges, should be every even column
      tmp_start <- single_chrom$start[seq(2, nrow(single_chrom)-1, 3)]
      
      # get end points for new ranges, should be every odd column
      tmp_end <- single_chrom$end[seq(4, nrow(single_chrom)-1, 3)]
      
      # Combine
      tmp <- cbind(tmp_start, tmp_end) %>% data.table::as.data.table()
      
      # rbind on last unpaired range
      lala <- single_chrom[nrow(single_chrom),2:3]
      
      # Combine with triplet ranges
      tmp <- rbind(tmp, lala, use.names = FALSE)
      
    } else if (rows_working_with %% 3 == 2){
      
      print(message("There are an even number of rows and \n and there are two additional rows"))
      
      # If there is an even number of columns, remove the last 1 or 2 row (i.e. make the df have an odd # of rows)
      # use the same methods as above, make the last/removed row it's own range (don't pair it with another range)
      # Get starting points of new ranges, should be every even column
      tmp_start <- single_chrom$start[seq(2, nrow(single_chrom)-2, 3)]
      
      # get end points for new ranges, should be every odd column
      tmp_end <- single_chrom$end[seq(4, nrow(single_chrom)-2, 3)]
      
      # Combine
      tmp <- cbind(tmp_start, tmp_end) %>% data.table::as.data.table()
      
      # If there are two ranges remaining, turn into a single range
      lala <- data.frame(single_chrom$start[nrow(single_chrom)-1], single_chrom$end[nrow(single_chrom)]) %>% data.table::as.data.table()
      
      # Change colnames
      colnames(lala) <- c("start", "end")
      
      # Combine with triplet ranges
      tmp <- rbind(tmp, lala, use.names = FALSE)
    }
    
    # Add back chromosome info
    tmp$seqid <- rep(i, nrow(tmp))
    
    # Rearrage, rename columns
    colnames(tmp) <- c("start", "end", "seqid")
    tmp <- tmp %>% select("seqid", "start", "end")
    
    # Add to outside df
    triplet_ranges <- rbind(triplet_ranges, tmp)
  }
  
  return(triplet_ranges)
}
second_triplet_haps <- second_triplet_window(ranges) 
second_triplet_haps$rr_id <- paste0("second_triplet_rr_", seq(1, nrow(second_triplet_haps)))

# join all ranges
all_triplets <- rbind(first_triplet_haps, second_triplet_haps)
