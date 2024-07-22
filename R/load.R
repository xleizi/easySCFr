read_matrice_matrix <- function(h5data) {
  if ("data" %in% names(h5data) && "indices" %in% names(h5data) && "indptr" %in% names(h5data)) {
    i <- h5data[["indices"]][]
    p <- h5data[["indptr"]][]
    x <- h5data[["data"]][]

    X <- if (getEncodingType(h5data) == "csr_matrix") {
      Matrix::sparseMatrix(i = i, p = p, x = x, dims = rev(getShape(h5data)), index1 = FALSE)
    } else {
      Matrix::t(Matrix::sparseMatrix(j = i, p = p, x = x, dims = getShape(h5data), index1 = FALSE))
    }

    return(X)
  } else {
    return(h5data[])
  }
}

read_matrix <- function(h5layer) {
  data_list <- list()
  for (name in names(h5layer)) {
    group_h5layer <- h5layer[[name]]
    if (getEncodingType(h5layer) == "array") {
      data_list[[name]] <- group_h5layer[, ]
    } else {
      data_list[[name]] <- read_matrice_matrix(group_h5layer)
    }
  }
  all_data <- do.call(cbind, data_list)
  return(all_data)
}

h5_to_X <- function(h5, assay = "RNA", layer = "rawdata") {
  h5layer <- h5[["assay"]][[assay]][["layers"]][[layer]]
  varType <- if (layer == "rawdata") "rawvar" else "var"
  X <- read_matrix(h5layer)
  colnames(X) <- h5[["obs"]][["_index"]][]
  rownames(X) <- h5[["var"]][[varType]][["_index"]][]
  return(X)
}

h5_to_DF <- function(h5data) {
  rownamesStr <- NULL
  for (name in names(h5data)) {
    if (name == "_index") {
      rownamesStr <- h5data[[name]][]
      next
    }
  }
  df <- data.frame()
  for (name in names(h5data)) {
    if (name == "_index") {
      next
    }
    if (getEncodingType(h5data[[name]]) == "categorical") {
      values_attr <- h5data[[name]]
      labelName <- values_attr[["categories"]][]
      values <- values_attr[["codes"]][]
      values <- factor(as.integer(values), labels = labelName)
    } else if (getEncodingType(h5data[[name]]) %in% c("array", "string-array")) {
      values <- h5data[[name]][]
    } else {
      stop("Unknown encoding type")
    }
    df <- addDF(df, setNames(data.frame(values), name), "col")
  }
  if (!is.null(rownamesStr)) {
    rownames(df) <- rownamesStr
  }
  if ("column-order" %in% hdf5r::h5attr_names(h5data)) {
    colnamesOrder <- hdf5r::h5attr(h5data, "column-order")
    df <- df[, colnamesOrder]
  }

  return(df)
}

openH5 <- function(FileName) {
  if (!hdf5r::is_hdf5(FileName)) {
    stop("File is not a hdf5 file.")
  }
  h5 <- hdf5r::H5File$new(FileName, "r")
  if (!("var" %in% names(h5))) {
    stop("var not found.")
  }
  if (!("obs" %in% names(h5))) {
    stop("obs not found.")
  }
  if (!("assay" %in% names(h5))) {
    stop("assay not found.")
  }
  return(h5)
}

# 读取 uns 数据
h5_to_uns <- function(h5obj) {
  data <- list()

  # 检查 h5obj 是否为组
  if (is(h5obj, "H5Group")) {
    keys <- names(h5obj)
    for (key in keys) {
      # 递归读取
      data[[key]] <- h5_to_uns(h5obj[[key]])
    }
  } else { # h5obj 是数据集
    item <- h5obj[]
    if (is.raw(item[1])) {
      item_text <- rawToChar(item)
      # 尝试将内容作为 JSON 解析，否则保持原样
      data <- tryCatch(
        {
          jsonlite::fromJSON(item_text)
        },
        error = function(e) {
          item_text
        }
      )
    } else {
      data <- item
    }
  }

  return(data)
}

h5_to_obs <- function(h5, obsName = "obs") {
  return(h5_to_DF(h5[["obs"]]))
}

sce_add_h5_to_graphs <- function(sce, h5, cellNames, graphsName = "graphs") {
  for (name in names(h5[[graphsName]])) {
    tryCatch({
      Graphmt <- read_matrice_matrix(h5[[graphsName]][[name]])
      rownames(Graphmt) <- cellNames
      colnames(Graphmt) <- cellNames
      sce@graphs[[name]] <- Seurat::as.Graph(Graphmt)
    }, error = function(e) {
      print(paste0("Error in reading graph ", name, ": ", e$message))
    })
  }
  return(sce)
}

sce_add_h5_to_var <- function(sce, h5, assay, varName = "var", SeuratVersion = checkSeuratVersion()) {
  if (SeuratVersion == 5) {
    sce@assays[[assay]] <- Seurat::AddMetaData(sce@assays[[assay]], h5_to_DF(h5[[varName]][["rawvar"]]))
  } else if (SeuratVersion == 4) {
    sce@assays[[assay]]@meta.features <- h5_to_DF(h5[[varName]][["rawvar"]])
  } else {
    stop("Unsupported Seurat version")
  }
  return(sce)
}

sce_add_h5_to_reductions <- function(sce, h5, cellNames, assay = "RNA", reductionsName = "reductions") {
  for (name in names(h5[[reductionsName]])) {
    tryCatch({
      matrix <- h5[[reductionsName]][[name]][, ]
      colnames(matrix) <- cellNames
      sce@reductions[[name]] <- Seurat::CreateDimReducObject(
        embeddings = t(matrix),
        key = paste0(name, "_"),
        assay = assay
      )
    }, error = function(e) {
      print(paste0("Error in reading reductions ", name, ": ", e$message))
    })
  }
  return(sce)
}

#' readH5
#'
#' @param FileName
#' @param assay
#' @param SeuratVersion
#' @param calData
#' @param calScale
#' @param calFeatures
#' @param readType
#'
#' @return
#' @export
#'
#' @examples
readH5 <- function(FileName,
                   assay = "RNA",
                   SeuratVersion = checkSeuratVersion(),
                   calData = TRUE,
                   calScale = FALSE,
                   calFeatures = FALSE,
                   datatype = "Seurat") {
  options(warn = -1)
  h5 <- openH5(FileName)
  tryCatch(
    {
      sce <- h5_to_seurat(h5, assay, SeuratVersion, calData, calScale, calFeatures)
    },
    error = function(e) {
      print(e)
    },
    finally = {
      h5$close_all()
    }
  )

  if (datatype == "Seurat") {
    return(sce)
  }
}

h5_to_seurat <- function(h5,
                         assay = "RNA",
                         SeuratVersion = checkSeuratVersion(),
                         calData = TRUE,
                         calScale = FALSE,
                         calFeatures = FALSE) {
  cellNames <- h5[["names_obs"]][]
  geneNames <- h5[["names_var"]][]

  # 1 如果存在data, data -> sce@data, rawdata -> sce@counts,var的行名 -> sce@assays[["RNA"]]@var.feature
  # 2 如果不存在data, rawdata -> sce@counts
  # assay
  layersNames <- names(h5[["assay"]][[assay]][["layers"]])
  rawX <- h5_to_X(h5, assay, "rawdata")
  sce <- Seurat::CreateSeuratObject(counts = rawX, assay = assay, project = "scanpy")
  if (calData) {
    sce <- Seurat::NormalizeData(sce)
  }

  # Features
  if (calFeatures) {
    if (is.null(sce@assays[[assay]]@layers$scale.data)) {
      if ("data" %in% layersNames) {
        print("there is no scale data in h5 files and calScale is FALSE, skip FindVariableFeatures")
        Seurat::VariableFeatures(sce) <- rownames(h5_to_DF(h5[["var"]][["var"]]))
      } else {
        print("there is no scale data in h5 files and calScale is FALSE, please set both 'calData = TRUE' and 'calScale = TRUE' to calculate scale data")
      }
    } else {
      sce <- Seurat::FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
    }
  } else {
    if ("data" %in% layersNames) {
      print("load features from h5 file")
      Seurat::VariableFeatures(sce) <- rownames(h5_to_DF(h5[["var"]][["var"]]))
    } else {
      print("there is no features in h5 files and calFeatures is FALSE, please set 'calFeatures = TRUE' to calculate features")
    }
  }

  # var
  sce <- sce_add_h5_to_var(sce, h5, assay, "var")

  # 如果calScale为TRUE，则计算scale数据，不读取data
  # 如果calScale为FALSE，且data存在，则读取data，不计算scale数据
  # 如果calScale为FALSE，且data存在，则读取data，判断dim是否一致，不一致则抛出提醒，是否calScale
  # 如果calScale为FALSE，且data不存在，则抛出提醒，是否calScale
  # Scale
  if (calScale) {
    if (calData) {
      sce <- Seurat::ScaleData(sce)
    } else {
      if ("data" %in% layersNames) {
        print("please set 'calData = TRUE' if you want to calculate scale data, load scale data from h5 file")
        scaleX <- h5_to_X(h5, assay, "data")
        sce@assays[[assay]]@layers$scale.data <- NULL
        sce@assays[[assay]]@layers$scale.data <- scaleX
      } else {
        stop("there is no scale data in h5 files and calScale is FALSE,
          please set both 'calData = TRUE' and 'calScale = TRUE'
          to calculate scale data")
      }
    }
  } else {
    if ("data" %in% layersNames) {
      scaleX <- h5_to_X(h5, assay, "data")
      sce@assays[[assay]]@layers$scale.data <- NULL
      sce@assays[[assay]]@layers$scale.data <- scaleX
      if (all(dim(scaleX) != dim(sce))) {
        print("the dim of scale data is not equal to the dim of raw data, please set 'calScale = TRUE' to calculate new scale data")
      }
    } else {
      print("there is no scale data in h5 files and calScale is FALSE, please set both 'calData = TRUE' and 'calScale = TRUE' to calculate scale data")
    }
  }

  # obs
  sce <- Seurat::AddMetaData(sce, h5_to_obs(h5, "obs"))

  # graphs
  if ("graphs" %in% names(h5)) {
    sce <- sce_add_h5_to_graphs(sce, h5, cellNames, "graphs")
  }

  # reductions
  if ("reductions" %in% names(h5)) {
    sce <- sce_add_h5_to_reductions(sce, h5, cellNames, assay, reductionsName = "reductions")
  }

  # JoinLayers
  if (SeuratVersion == 5) {
    if(exists("scaleX") && all(dim(scaleX) == dim(sce))){
      sce <- SeuratObject::JoinLayers(sce, assay = assay)
    }
  }

  # 添加uns
  uns <- h5_to_uns(h5[["uns"]])
  for (unsName in names(uns)) {
    # class(uns[[unsName]]) <- "SeuratCommand"
    # attr(uns[[unsName]], "package") <- "SeuratObject"
    sce@commands[[unsName]] <- uns[[unsName]]
  }

  return(sce)
}
