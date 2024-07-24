df_to_h5 <- function(h5data, df) {
  h5AddAttribute(h5data, "encoding-type", "dataframe")
  h5AddAttribute(h5data, "shape", dim(df))
  h5AddAttribute(h5data, "column-order", colnames(df))
  h5AddAttribute(h5data, "_index", "_index")
  h5AddAttribute(h5data, "encoding-version", "0.2.0")
  h5dfAddIndex(h5data, df)
  for (name in colnames(df)) {
    if (name == "_index") {
      next
    }
    if (is.character(df[[name]]) & length(unique(df[[name]])) < 50) {
      df[[name]] <- as.factor(df[[name]])
    }
    if (is.factor(df[[name]])) {
      h5data2 <- h5data$create_group(name)
      h5data2[["categories"]] <- levels(df[[name]])
      h5AddAttribute(h5data2[["categories"]], "encoding-type", "string-array")
      h5AddAttribute(h5data2[["categories"]], "encoding-version", "0.2.0")

      h5data2[["codes"]] <- as.integer(df[[name]]) - 1L
      h5AddAttribute(h5data2[["codes"]], "encoding-type", "array")
      h5AddAttribute(h5data2[["codes"]], "encoding-version", "0.2.0")

      h5AddAttribute(h5data2, "encoding-type", "categorical")
      h5AddAttribute(h5data2, "encoding-version", "0.2.0")
      h5AddAttribute(h5data2, "ordered", is.ordered(df[[name]]))
    } else {
      h5data[[name]] <- df[[name]]
      h5AddAttribute(h5data[[name]], "encoding-type", "array")
      h5AddAttribute(h5data[[name]], "encoding-version", "0.2.0")
    }
  }
}


X_to_h5 <- function(
    sce_X, h5data, data_name = "data",
    split_save = TRUE,
    max_cells_per_subset = 5000) {
  if (nrow(sce_X) == 0 | ncol(sce_X) == 0) {
    return(NULL)
  }
  if ("dgCMatrix" %in% class(sce_X)) {
    sce_X <- Matrix::t(sce_X)
    sce_X <- as(sce_X, "RsparseMatrix")
  }
  if ("Graph" %in% class(sce_X)) {
    sce_X <- as(sce_X, "RsparseMatrix")
  }
  if ("matrix" %in% class(sce_X) | attr(class(sce_X), "package") == "BPCells") {
    sce_X <- t(sce_X)
  }
  handle_data_splitting(sce_X, h5data, data_name, split_save, max_cells_per_subset)
}

reductions_to_h5 <- function(sce, h5file, group_path) {
  reductionsh5 <- h5file$create_group(group_path)
  h5AddAttribute(reductionsh5, "encoding-type", "dict")
  h5AddAttribute(reductionsh5, "encoding-version", "0.1.0")
  for (name in names(sce@reductions)) {
    reductionsh5[[name]] <- t(sce@reductions[[name]]@cell.embeddings)
    h5AddAttribute(reductionsh5[[name]], "encoding-type", "array")
    h5AddAttribute(reductionsh5[[name]], "encoding-version", "0.2.0")
  }
}

graphs_to_h5 <- function(sce, h5file, group_path) {
  graphsh5 <- h5file$create_group(group_path)
  for (name in names(sce@graphs)) {
    tryCatch(
      {
        graphs_X <- as(sce@graphs[[name]], "RsparseMatrix")
        write_matrix(graphsh5, matrix = graphs_X, name = name)
      },
      error = function(e) {
        print(paste0("Error in writing graph: ", name, e))
      }
    )
  }
}

handle_data_splitting <- function(
    sce_X, h5data, data_name,
    split_save = TRUE,
    max_cells_per_subset = 5000) {
  h5data_Name <- h5data$create_group(data_name)
  h5AddAttribute(h5data_Name, "encoding-type", "csr_matrix")
  h5AddAttribute(h5data_Name, "encoding-version", "0.1.0")
  h5AddAttribute(h5data_Name, "shape", dim(sce_X))


  num_cells <- nrow(sce_X)
  if (num_cells > 200000 & !split_save) {
    print("Large dataset detected. Splitting the save process.")
    split_save <- TRUE
  }
  if (split_save) {
    num_subsets <- ceiling(num_cells / max_cells_per_subset)
    # i <- 1
    for (i in 1:num_subsets) {
      start_idx <- ((i - 1) * max_cells_per_subset) + 1
      end_idx <- min(i * max_cells_per_subset, num_cells)
      subset <- sce_X[start_idx:end_idx, ]
      save_path <- sprintf("%s/%s_%05d.npz", data_name, data_name, i - 1)
      write_matrix(h5data, matrix = subset, name = save_path)
    }
  } else {
    write_matrix(h5data, matrix = sce_X, name = sprintf("%s/%s_00000.npz", data_name, data_name))
  }
}

# 插入稀疏矩阵
write_matrix <- function(h5data, matrix, name) {
  if (attr(class(matrix), "package") == "BPCells") {
    # matrix <- as(matrix, "dgRMatrix")
    # matrix <- as(matrix, "RsparseMatrix")
    matrix <- as(matrix, "sparseMatrix")
    matrix <- as(matrix, "RsparseMatrix")
  }
  if ("matrix" %in% class(matrix)) {
    h5data[[name]] <- matrix
    h5AddAttribute(h5data[[name]], "encoding-type", "array")
    h5AddAttribute(h5data[[name]], "encoding-version", "0.2.0")
    h5AddAttribute(h5data[[name]], "shape", dim(matrix))
  } else if ("dgRMatrix" %in% class(matrix)) {
    h5data_sub <- h5data$create_group(name)
    h5data_sub[["data"]] <- slot(object = matrix, name = "x")
    h5data_sub[["indices"]] <- slot(object = matrix, name = "j")
    h5data_sub[["indptr"]] <- slot(object = matrix, name = "p")
    h5AddAttribute(h5data_sub, "encoding-type", "csr_matrix")
    h5AddAttribute(h5data_sub, "encoding-version", "0.1.0")
    h5AddAttribute(h5data_sub, "shape", dim(matrix))
  } else {
    print(paste0("Error in writing matrix: ", class(matrix)))
  }
}



# 将数据写入 HDF5 文件
commands_to_h5 <- function(data, h5file, group_path) {
  for (i in 1:length(data)) {
    if ("SeuratCommand" %in% class(data[[i]])) {
      tryCatch(
        {
          data[[i]] <- as.list(data[[i]])
        },
        error = function(e) {
          print(paste0("Error in converting SeuratCommand to list: ", e))
        }
      )
    }
  }
  # 创建或获取组
  if (!h5file$exists(group_path)) {
    group <- h5file$create_group(group_path)
  } else {
    group <- h5file$open_group(group_path)
  }

  # key <- names(data)[1]
  # 遍历列表数据
  for (key in names(data)) {
    item <- data[[key]]
    # 如果元素是列表，递归调用函数
    if (is.list(item)) {
      new_group_path <- paste0(group_path, "/", key)
      commands_to_h5(item, h5file, new_group_path)
    } else if (is.function(item)) {
      item <- capture.output(item)
      item <- paste(item, collapse = "\n")
      group[[key]] <- item
    } else if (is.function(item)) {
      data_type <- h5file$create_type("string", size = nchar(item))
      dataset <- group$create_dataset(key, dtype = data_type, dims = 1)
      group[[key]] <- item
    } else {
      group[[key]] <- item
    }
  }
}

#' saveH5
#'
#' @param FileName
#' @param data
#' @param assay
#' @param save_graph
#' @param SeuratVersion
#' @param split_save
#' @param max_cells_per_subset
#'
#' @return
#' @export
#'
#' @examples
saveH5 <- function(FileName, data, assay = "RNA",
                   save_graph = TRUE,
                   SeuratVersion = checkSeuratVersion(),
                   split_save = TRUE,
                   max_cells_per_subset = 5000) {
  if ("Seurat" %in% class(data)) {
    Seurat_to_H5(FileName, data, assay = assay, save_graph = save_graph, SeuratVersion = SeuratVersion, split_save = split_save, max_cells_per_subset = max_cells_per_subset)
  }
}

Seurat_to_H5 <- function(FileName, sce, assay = "RNA",
                         save_graph = TRUE,
                         SeuratVersion = checkSeuratVersion(),
                         split_save = TRUE,
                         max_cells_per_subset = 5000) {
  if (SeuratVersion == 5) {
    if ("Assay" %in% class(sce[[assay]])) {
      sce[[assay]] <- as(object = sce[[assay]], Class = "Assay5")
    }
    sce <- SeuratObject::JoinLayers(sce, assay = assay)
  }
  h5 <- hdf5r::H5File$new(FileName, "w")

  # X
  assayList <- h5$create_group("assay")
  assayListsCur <- assayList$create_group(assay)
  layersList <- assayListsCur$create_group("layers")


  if (SeuratVersion == 4) {
    rawData <- SeuratObject::GetAssayData(object = sce, assay = assay, slot = "counts")
    scaleData <- SeuratObject::GetAssayData(object = sce, assay = assay, slot = "scale.data")
    # var
    varObjName <- "meta.features"
  } else if (SeuratVersion == 5) {
    rawData <- SeuratObject::GetAssayData(object = sce, assay = assay, layer = "counts")
    scaleData <- SeuratObject::GetAssayData(object = sce, assay = assay, layer = "scale.data")
    # var
    varObjName <- "meta.data"
  }

  # sce_X = rawData
  # h5data = layersList
  # data_name = "rawdata"
  X_to_h5(
    sce_X = rawData, h5data = layersList, data_name = "rawdata",
    split_save = split_save, max_cells_per_subset = max_cells_per_subset
  )

  X_to_h5(
    sce_X = scaleData, h5data = layersList, data_name = "data",
    split_save = split_save, max_cells_per_subset = max_cells_per_subset
  )

  # var
  varList <- h5$create_group("var")

  rawvarh5 <- varList$create_group("rawvar")
  rawvar <- slot(object = sce@assays[[assay]], name = varObjName)
  df_to_h5(rawvarh5, rawvar)

  if (!is.null(Seurat::VariableFeatures(sce))) {
    var <- rawvar[rownames(rawvar) %in% Seurat::VariableFeatures(sce), ]
    varh5 <- varList$create_group("var")
    df_to_h5(varh5, var)
  }

  # obs
  obs <- h5$create_group("obs")
  df_to_h5(obs, sce@meta.data)

  # reductions
  if (length(sce@reductions) > 0) {
    reductions_to_h5(sce, h5, "reductions")
    # tryCatch(
    #     {
    #         reductions_to_h5(sce, h5, "reductions")
    #     },
    #     error = function(e) {
    #         print(paste0("Error in writing reductions: ", e))
    #     }
    # )
  }


  # graphs
  if (length(sce@graphs) > 0 & save_graph) {
    graphs_to_h5(sce, h5, "graphs")
  }

  # names
  h5[["names_var"]] <- rownames(sce)
  h5[["names_obs"]] <- colnames(sce)

  commands_to_h5(sce@commands, h5, "uns")
  h5$close_all()
}
