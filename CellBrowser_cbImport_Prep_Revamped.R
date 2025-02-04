# Tyler Therron, MS

#LIBRARY PATH
libs <- .libPaths("/path/")

library(Seurat)
library(dplyr)
library(tidyr)
library(readr)

# INPUT for script
args <- commandArgs(TRUE) # allows script to run in CmdLine

seurat_object_final <- as.character(args[1]) # INPUT Seurat S4 .rds object locale
desired_output_name <- as.character(args[2]) # NAME of Seurat object that is output
output_directory <- as.character(args[3]) # OUTPUT locale - include the '/' at the end because of how code writes the output 

#  ================================= ================================= ================================= ================================= ================================= ================================= ================================= ================================= =================================
# seurat_object_final <- as.character("CD64pos_CITEseq_edits2.Rds") # INPUT Seurat S4 .rds object locale
# desired_output_name <- as.character("TEST_cbPrep") # NAME of Seurat object that is output
# default_ident <- as.character("seurat_clusters") # HAS to be a FACTOR, will be loaded identity CB
# output_directory <- as.character("/Users/ttm3567/Documents/August2024/")
# ===================================================================================================

#. -------------- FXNs --------------
# Function to check if the data contains non-integer values
is_normalized <- function(seurat_obj) {
  rna_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  return(any(rna_data != round(rna_data))) # Returns TRUE if any non-integers are found
}

check_seurat_version <- function(seurat_obj) {
  # Check if 'layers' exists in the RNA assay
  if ("layers" %in% slotNames(seurat_obj@assays[["RNA"]])) {
    return("v5")
  } else {
    return("v4")
  }
}

# Function to check for the presence of ADT assay
process_ADT_assay <- function(seurat_obj) {
  if ("ADT" %in% names(seurat_obj@assays)) {
    adt_assay <- GetAssayData(seurat_obj, assay = "ADT")
    print("grabbed ADT assay")
    
    # Pull name of ADT features
    adt_features <- rownames(adt_assay)
    print("grabbed ADT names")
    
    # Append "ADT_" to the beginning of all the ADT feature names
    adt_features <- paste0("ADT_", adt_features)
    print("appended _ADT to features")
    
    # Replace ADT assay feature names with edited feature names
    rownames(adt_assay) <- adt_features
    print("replaced feature ADT names")
    
    # Check if Seurat object is v5 or v4
    if ("layers" %in% slotNames(seurat_obj@assays$RNA)) {
      # Seurat v5 logic
      combined_assay <- rbind(GetAssayData(seurat_obj, assay = "RNA", slot = "counts"), adt_assay)
      seurat_obj@assays[["RNA"]]@layers[["counts"]] <- combined_assay
      print("Seurat v5: Replaced RNA assay counts in layers")
    } else {
      # Seurat v4 logic
      combined_assay <- rbind(GetAssayData(seurat_obj, assay = "RNA", slot = "counts"), adt_assay)
      seurat_obj@assays$RNA@counts <- combined_assay
      print("Seurat v4: Replaced RNA assay counts")
    }
  } else {
    print("only RNA assay so not combining any ADT counts")
  }
  
  return(seurat_obj)
}

align_cells <- function(seurat_obj) {
  if ("ADT" %in% names(seurat_obj@assays)) {
    print("if 1")
    # Check if Seurat object is v5 or v4
    if ("layers" %in% slotNames(seurat_obj@assays$RNA)) {
      # Seurat v5 logic
      print("if 2")
      rna_cells <- colnames(seurat_obj@assays[["RNA"]]@layers[["counts"]])
      adt_cells <- colnames(seurat_obj@assays[["ADT"]]@counts)
    } else {
      print("else 2")
      # Seurat v4 logic
      rna_cells <- colnames(seurat_obj@assays$RNA@counts)
      adt_cells <- colnames(seurat_obj@assays$ADT@counts)
    }
    
    # Identify cells present in RNA but missing in ADT
    missing_cells <- setdiff(rna_cells, adt_cells)
    print("missing_cells")
    if (length(missing_cells) > 0) {
      if ("layers" %in% slotNames(seurat_obj@assays$RNA)) {
        print("if 3")
        # Seurat v5 logic
        zeros_matrix <- matrix(0, nrow = nrow(seurat_obj@assays[["ADT"]]@counts), ncol = length(missing_cells))
        rownames(zeros_matrix) <- rownames(seurat_obj@assays[["ADT"]]@counts)
        colnames(zeros_matrix) <- missing_cells
        
        new_adt_counts <- cbind(seurat_obj@assays[["ADT"]]@counts, zeros_matrix)
        new_adt_counts <- new_adt_counts[, rna_cells]
        
        seurat_obj@assays[["ADT"]]@counts <- new_adt_counts
      } else {
        print("else 3")
        # Seurat v4 logic
        zeros_matrix <- matrix(0, nrow = nrow(seurat_obj@assays$ADT@counts), ncol = length(missing_cells))
        rownames(zeros_matrix) <- rownames(seurat_obj@assays$ADT@counts)
        colnames(zeros_matrix) <- missing_cells
        
        new_adt_counts <- cbind(seurat_obj@assays$ADT@counts, zeros_matrix)
        new_adt_counts <- new_adt_counts[, rna_cells]
        
        seurat_obj@assays$ADT@counts <- new_adt_counts
	seurat_obj@assays$ADT@data <- new_adt_counts
      }
      
      print("Aligned ADT cells with RNA cells")
    } else {
      print("No missing cells, ADT and RNA assays are aligned")
    }
  } else {
    print("ADT assay not present, skipping alignment")
  }
  
  return(seurat_obj)
}

choose_valid_ident <- function(seurat_obj) {
  # Get the metadata columns that are factors
  factor_columns <- names(Filter(is.factor, seurat_obj@meta.data))
  
  # Check each factor column for the number of levels
  for (column in factor_columns) {
    levels_count <- length(levels(seurat_obj@meta.data[[column]]))
    
    # If a factor column has more than one level, use it as the default ident
    if (levels_count > 1) {
      Idents(seurat_obj) <- seurat_obj@meta.data[[column]]
      print(paste("Default Ident set to:", column))
      return(column)
    }
  }
  
  stop("No valid metadata column with more than one factor level found.")
}
#. -------------- FXNs --------------
# ===================================================================================================

seurat_object_readIN <- readRDS(seurat_object_final)
print("read in data")

seurat_object_readIN@active.assay <- 'RNA'
print("active assay is forsure RNA now")

default_ident <- seurat_object_readIN@active.ident

# Proceed with attaching markers if the ident was successfully set
seurat_object_readIN@misc$markers <- FindAllMarkers(seurat_object_readIN)
print("attached markers")


Seurat_Obj <- seurat_object_readIN
# --- adding in ADT protein features ---
# need to normalize the data first

# Apply normalization only if the data is not already normalized
if (!is_normalized(Seurat_Obj)) {
  Seurat_Obj <- NormalizeData(Seurat_Obj, normalization.method = "LogNormalize", scale.factor = 10000)
  print("Data normalized")
} else {
  print("Data is already normalized, skipping normalization step")
}

seurat_version <- check_seurat_version(Seurat_Obj)

if (seurat_version == "v5") {
  # For Seurat v5
  Seurat_Obj@assays[["RNA"]]@layers[["counts"]] <- Seurat_Obj@assays[["RNA"]]@layers[["data"]]
  print("Seurat v5: layers slot used for counts")
} else {
  # For Seurat v4
  Seurat_Obj@assays$RNA@counts <- Seurat_Obj@assays[["RNA"]]@data
  print("Seurat v4: counts slot used for data")
}

# Align ADT assay with RNA cells if ADT assay is present
Seurat_Obj <- align_cells(Seurat_Obj)

# combine ADT and RNA assays:
Seurat_Obj <- process_ADT_assay(Seurat_Obj)

# 11. Inspect results
if (seurat_version == "v5") {
  # For Seurat v5
  tail(rownames(data.frame(Seurat_Obj@assays[["RNA"]]@layers[["counts"]])), n = 200)
} else {
  # For Seurat v4
  tail(rownames(data.frame(Seurat_Obj@assays$RNA@counts)), n = 200)
}


# 4. specify the default Cell Identities -> Idents(object) <- "metadata_column_name" - ** HAS to be last step and HAS to be a FACTOR **
Idents(Seurat_Obj) <- default_ident 
#Idents(Seurat_Obj)
print("default Ident set")

# 12. Save the seurat object as RDS
saveRDS(Seurat_Obj, paste0(output_directory,desired_output_name,".rds"))
print("successfully saved")


