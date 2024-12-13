# for a given BT output (in file identified by filename variable) get
# the node assigned by BT to the target classification string
# input and relevant results are in the commented lines after the script

get_node_id <- function(clstring, filename) {
  # read the file clstrings_ids.txt where Glottolog classification strings
  # are listed along with all member languages,
  # find the relevant string, and get the member languages as a vector
  # sl stands for classification strings and languages
  # the matching will not work in a few cases where a language is missing
  # from the BT output but present in clstrings_ids.txt; those are 
  # languages for which Glottolog does not provide coordinates
  # they can be found in glottolog_classification_strings.txt, where 
  # coordinates will be given as NA
  # for this reason such languages are retried from 
  # glottolog_classification_strings.txt and removed when matching 
  # vectors containing languages
  gcs <- read.table(file="glottolog_classification_strings.txt", header=TRUE,
                    sep="\t", quote="")
  no_coor <- gcs$GlottologID[is.na(gcs$lat)]
  sl <- read.table(file="clstrings_ids.txt", sep="\t", header=TRUE, quote="")
  w_s <- which(sl$clades==clstring)
  if (length(w_s)==0) {
    stop("Could not find the classification string")
  }
  id_string_file <- sl$ids[w_s]
  id_vector_file_full <- strsplit(id_string_file, "_")[[1]]
  remove <- which(id_vector_file_full %in% no_coor)
  if (length(remove) > 0) {
    id_vector_file <- id_vector_file_full[-remove]
  } else {
    id_vector_file <- id_vector_file_full
  }
  # now read the BT output file and find the node number
  BT_all <- readLines(paste0("BT_outfiles/", filename))
  start_discard <- grep("Itter", BT_all)
  BT <- BT_all[1:(start_discard - 1)]
  parse_line <- function(line) {
    l1 <- strsplit(line, "\\\t")[[1]]
    languages <- l1[3:length(l1)]
    l2 <- l1[1]
    l3 <- strsplit(l2, "Node-")[[1]][2]
    node <- as.numeric(l3)
    return(list(node,languages))
  }
  found <- FALSE
  j <- 0
  while(found==FALSE) {
    j <- j + 1
    if(j==(length(BT)+1)) {
      cat("A node could not be identified containing\n")
      print(id_vector_file)
      return(NULL)
    }
    parse_line_out <- parse_line(BT[j])
    foundnode <- parse_line_out[[1]]
    lgs <- parse_line_out[[2]]
    if (setequal(lgs,id_vector_file)==TRUE) {
      found==TRUE
      return(foundnode)
    }
  }
  cat("Node could not be identified\n")
  return(NULL)
}

# get_node_id(clstring, filename)
# filename <- "indo1319_AncStates.txt"
# clstring <- "Indo-European,ClassicalIndo-European"  # 19
# clstring <- "Indo-European,ClassicalIndo-European,Indo-Iranian"  # 269
# clstring <- "Indo-European,ClassicalIndo-European,Balto-Slavic"  # 30
# clstring <- "Indo-European,ClassicalIndo-European,Albanian"  # 20

# filename <- "atla1278_AncStates.txt"
# clstring <- "Atlantic-Congo"  # 0
# clstring <- "Atlantic-Congo,Volta-Congo"  # 96
# clstring <- "Atlantic-Congo,Mel"  # 5
# clstring <- "Atlantic-Congo,North-CentralAtlantic"  # 20
# clstring <- "Atlantic-Congo,Volta-Congo,NorthVolta-Congo" # 1731
# clstring <- "Atlantic-Congo,Volta-Congo,Benue-Congo"  # 97

# filename <- "afro1255_AncStates.txt"
# clstring <- "Afro-Asiatic"  # 0
# clstring <- "Afro-Asiatic,Chadic"  # 36
# clstring <- "Afro-Asiatic,Cushitic"  # 367

# filename <- "araw1281_AncStates.txt"
# clstring <- "Arawakan"  # 0
# clstring <- "Arawakan,SouthernMaipuran"  # 83
# clstring <- "Arawakan,SouthernMaipuran,Kampa-Amuesha"  # 93
# clstring <- "Arawakan,Central-EasternMaipuran"  # 22
# clstring <- "Arawakan,CaribbeanArawakan,AntilleanArawakan"  # 7
# clstring <- "Arawakan,CaribbeanArawakan,Guajiro-Paraujano"  # 16
# clstring <- "Arawakan,SouthernMaipuran,BolivianArawakan"  # 84
# clstring <- "Arawakan,SouthernMaipuran,Kampa-Amuesha,Pre-AndineMaipuran"  # 94

# filename <- "chib1249_AncStates.txt"
# clstring <- "Chibchan,CoreChibchan"
# clstring <- "Chibchan,CoreChibchan,Magdalenic"  # 17
# clstring <- "Chibchan,CoreChibchan,IsthmicChibchan,EasternIsthmicChibchan"
# clstring <- "Chibchan,CoreChibchan,IsthmicChibchan"  # 1
# clstring <- "Chibchan,CoreChibchan,Magdalenic,SouthernMagdalenic"  # 27
# clstring <- "Chibchan,CoreChibchan,Magdalenic,NorthernMagdalenic"  # 18
# clstring <- "Chibchan,CoreChibchan,IsthmicChibchan,EasternIsthmicChibchan,Guaymiic"  # 5
# clstring <- "Chibchan,CoreChibchan,IsthmicChibchan,EasternIsthmicChibchan,Kuna"  # 8



