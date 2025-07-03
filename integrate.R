library(Seurat)


DCIS_SCT_path <- "/share/crsp/lab/dalawson/share/DCIS_integration/Merged_DCIS_all_SCT.rds"
DCIS_LT_path <- "/share/crsp/lab/dalawson/share/1_Lab_sequencing/Sharmila/DCIS/Merged/Merged_DCIS_all.rds"

KCIS_SCT_path <- "/share/crsp/lab/dalawson/share/1_Lab_sequencing/Sharmila/DCIS/fariba.rds"

DCIS_SCT <- readRDS(DCIS_SCT_path)
DCIS_LogTransform <- readRDS(DCIS_LT_path)
KCIS_SCT <- readRDS(KCIS_SCT_path)


#head(DCIS_LogTransform)
#DCIS_LogTransform@assays$RNA@counts
#DCIS_SCT_fixed <- SCTransform(DCIS_LogTransform, verbose = TRUE)
#saveRDS(DCIS_SCT_fixed, "/share/crsp/lab/dalawson/share/DCIS_integration/DCIS_SCT_fixed.rds")
#print(KCIS_SCT@assays$RNA@counts[1:5, 1:5])



#KCIS_SCT_fixed <- SCTransform(KCIS_SCT, verbose = TRUE)
#saveRDS(KCIS_SCT_fixed, "/share/crsp/lab/dalawson/share/DCIS_integration/KCIS_SCT_fixed.rds")


fix_metadata <- function(obj) {
  for (col in colnames(obj@meta.data)) {
    if (is.factor(obj@meta.data[[col]])) {
      obj@meta.data[[col]] <- as.character(obj@meta.data[[col]])
    }
  }
  return(obj)
}


get_safe_features <- function(obj) {
  modeled <- unique(unlist(lapply(obj[["SCT"]]@SCTModel.list, function(m) {
    rownames(m@feature.attributes)
  })))
  scaled <- rownames(obj[["SCT"]]@scale.data)
  intersect(modeled, scaled)
}



KCIS_SCT_fixed <- readRDS("/share/crsp/lab/dalawson/share/DCIS_integration/KCIS_SCT_fixed.rds")
DCIS_SCT_fixed <- readRDS("/share/crsp/lab/dalawson/share/DCIS_integration/DCIS_SCT_fixed.rds")

# check if anything is broken
#DCIS_SCT_test <- subset(DCIS_SCT_fixed, cells = sample(colnames(DCIS_SCT), 2000))
#DCIS_SCT_test <- fix_metadata(DCIS_SCT_test)

#safe_features <- get_safe_features(DCIS_SCT_test)
#features <- head(safe_features, 2000)

#obj_list <- list(DCIS_SCT_test, DCIS_SCT_test)

#obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)

#anchors <- FindIntegrationAnchors(
#  object.list = obj_list,
#  normalization.method = "SCT",
#  anchor.features = features
#)

#integrated <- IntegrateData(
#  anchorset = anchors,
#  normalization.method = "SCT"
#)

#print("DCIS is good")

#KCIS_SCT_test <- subset(KCIS_SCT_fixed, cells = sample(colnames(KCIS_SCT_fixed), 2000))
#KCIS_SCT_test <- fix_metadata(KCIS_SCT_test)

#safe_features <- get_safe_features(KCIS_SCT_test)
#features <- head(safe_features, 2000)

#obj_list <- list(KCIS_SCT_test, KCIS_SCT_test)

#obj_list <- PrepSCTIntegration(obj_list, anchor.features = features)

#anchors <- FindIntegrationAnchors(
#  object.list = obj_list,
#  normalization.method = "SCT",
#  anchor.features = features
#)

#integrated <- IntegrateData(
#  anchorset = anchors,
#  normalization.method = "SCT"
#)

#print("KCIS is good")

#DCIS_SCT_down <- subset(DCIS_SCT_fixed, cells = sample(colnames(DCIS_SCT), 10000))
#KCIS_SCT_down <- subset(KCIS_SCT_fixed, cells = sample(colnames(KCIS_SCT), 10000))



DCIS_SCT_down <- fix_metadata(DCIS_SCT_fixed)
KCIS_SCT_down <- fix_metadata(KCIS_SCT_fixed)

# check quality of data
head(DCIS_SCT_down@meta.data)
str(DCIS_SCT_down@meta.data)

head(KCIS_SCT_down@meta.data)
str(KCIS_SCT_down@meta.data)



SCT_objs <- list(DCIS_SCT_down, KCIS_SCT_down)

# integrate SCT data
scale_DCIS <- rownames(DCIS_SCT_down[["SCT"]]@scale.data)
scale_KCIS <- rownames(KCIS_SCT_down[["SCT"]]@scale.data)

usable_features <- intersect(scale_DCIS, scale_KCIS)
length(usable_features)


# Use them directly (skip SelectIntegrationFeatures which may reintroduce bad ones)
features <- usable_features[1:min(3000, length(usable_features))]
sct_genes_1 <- rownames(SCT_objs[[1]][["SCT"]])
sct_genes_2 <- rownames(SCT_objs[[2]][["SCT"]])

# Restrict features to only valid genes
features <- intersect(features, sct_genes_1)
features <- intersect(features, sct_genes_2)

# Confirm it's clean
length(features)
all(features %in% rownames(SCT_objs[[1]][["SCT"]]))
all(features %in% rownames(SCT_objs[[2]][["SCT"]]))


features <- as.character(features)
features <- features[!is.na(features)]
features <- unique(features)
str(features)
# Proceed with integration
SCT_objs <- list(DCIS_SCT_down, KCIS_SCT_down)
SCT_objs <- PrepSCTIntegration(object.list = SCT_objs, anchor.features = features)

#anchors <- FindIntegrationAnchors(
#  object.list = SCT_objs,
#  normalization.method = "SCT",
#  anchor.features = features
#)

#print("anchors found")
#integrated <- IntegrateData(
#  anchorset = anchors,
#  normalization.method = "SCT"
#)

#saveRDS(integrated, "/share/crsp/lab/dalawson/share/DCIS_integration/DCIS_KCIS_SCT_integration.rds")

# LogTransformed integration
KCIS_LogTransform <- NormalizeData(KCIS_SCT, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

DCIS_Log_down <- fix_metadata(DCIS_LogTransform)
KCIS_Log_down <- fix_metadata(KCIS_LogTransform)

# check quality of data
head(DCIS_Log_down@meta.data)
str(DCIS_Log_down@meta.data)

head(KCIS_Log_down@meta.data)
str(KCIS_Log_down@meta.data)



Log_objs <- list(DCIS_Log_down, KCIS_Log_down)

features <- SelectIntegrationFeatures(object.list = Log_objs, nfeatures = 3000)
Log_objs <- lapply(Log_objs, function(obj) {
  obj <- ScaleData(obj, features = features, verbose = FALSE)
  obj <- RunPCA(obj, features = features, verbose = FALSE)
  return(obj)
})


anchors <- FindIntegrationAnchors(
  object.list = Log_objs,
  anchor.features = features,
  normalization.method = "LogNormalize"
)
integrated <- IntegrateData(
  anchorset = anchors,
  normalization.method = "LogNormalize"
)

saveRDS(integrated, "/share/crsp/lab/dalawson/share/DCIS_integration/DCIS_KCIS_Log_integration.rds")
