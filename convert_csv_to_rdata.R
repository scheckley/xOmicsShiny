# Complete function to convert CSV files to RData format with biomaRt annotation
# Includes all logic from xOmicsShiny_ConvertCSV2Rdata.R and biomaRt from inputdata.R

convert_csv_to_rdata <- function(sample_file_path, exp_file_path, comp_file_path = NULL,
                                 protein_id_file_path = NULL, project_name, species = "human",
                                 auto_annotate = FALSE, id_type = NULL, include_description = TRUE,
                                 fill_missing_names = TRUE, create_network = FALSE,
                                 temp_dir = "unlisted/") {
  require(tidyverse)
  require(reshape2)
  if (create_network) require(Hmisc)
  if (auto_annotate) require(biomaRt)

  # Ensure temp directory exists
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }

  # Helper function to clean up empty rows/columns
  cleanup_empty <- function(df) {
    df.empty <- (is.na(df) | df == "")
    selCol <- !(colSums(df.empty) == nrow(df))
    selRow <- !(rowSums(df.empty) == ncol(df))
    return(df[selRow, selCol])
  }

  # Read and process sample metadata
  MetaData <- read.csv(sample_file_path, header = TRUE, check.names = FALSE)
  MetaData <- cleanup_empty(MetaData)

  # Read and process expression data
  exp_file <- exp_file_path
  if (str_detect(exp_file, "gz$")) {
    exp_file <- gzfile(exp_file, "rt")
  }
  if (str_detect(exp_file, "zip$")) {
    fnames <- as.character(unzip(exp_file, list = TRUE)$Name)
    exp_file <- unz(exp_file, fnames[1])
  }

  data_wide <- read.csv(exp_file, row.names = 1, header = TRUE, check.names = FALSE)
  data_wide <- cleanup_empty(data_wide)

  # Create data_long from data_wide
  data_long <- reshape2::melt(as.matrix(data_wide))
  colnames(data_long) <- c("UniqueID", "sampleid", "expr")
  data_long <- data_long %>%
    mutate(sampleid = as.character(sampleid)) %>%
    left_join(MetaData %>% dplyr::select(sampleid, group), by = "sampleid")

  # Process comparison data if provided
  results_long <- NULL
  tests <- c()

  if (!is.null(comp_file_path) && file.exists(comp_file_path)) {
    comp_file <- comp_file_path
    if (str_detect(comp_file, "gz$")) {
      comp_file <- gzfile(comp_file, "rt")
    }
    if (str_detect(comp_file, "zip$")) {
      fnames <- as.character(unzip(comp_file, list = TRUE)$Name)
      comp_file <- unz(comp_file, fnames[1])
    }

    results_long <- read.csv(comp_file, header = TRUE, check.names = FALSE)
    results_long <- cleanup_empty(results_long)
    tests <- sort(unique(results_long$test))
  }

  # Get all unique IDs from data
  if (!is.null(results_long) && nrow(results_long) > 0) {
    IDs <- rownames(data_wide)
    IDs2 <- results_long$UniqueID
    IDall <- unique(c(IDs, IDs2))
    IDall <- IDall[!is.na(IDall)]
    IDall <- IDall[!(IDall == "")]
  } else {
    IDall <- rownames(data_wide)
  }

  # Process protein/gene annotation
  ProteinGeneName <- NULL

  if (!auto_annotate && !is.null(protein_id_file_path) && file.exists(protein_id_file_path)) {
    # Use provided annotation file
    ProteinGeneName <- read.csv(protein_id_file_path, header = TRUE, check.names = FALSE)
    ProteinGeneName <- cleanup_empty(ProteinGeneName)
    ProteinGeneName <- ProteinGeneName %>% dplyr::filter(UniqueID %in% IDall)

    if (!"Protein.ID" %in% names(ProteinGeneName)) {
      ProteinGeneName$Protein.ID <- NA
    }
  } else if (auto_annotate && !is.null(id_type)) {
    # Auto-annotation using biomaRt

    if (str_detect(id_type, "UniProt")) {
      # Protein name matching using UniProt
      if (file.exists("db/ProteinInfo.rds")) {
        ProteinInfo <- readRDS("db/ProteinInfo.rds") %>%
          dplyr::mutate(Gene_Name = str_replace(Gene_Name, " .+", ""))

        if (id_type == "UniProtKB Protein ID") {
          ProteinGeneName <- data.frame(id = 1:length(IDall), UniqueID = IDall) %>%
            dplyr::left_join(ProteinInfo %>%
              dplyr::transmute(
                UniqueID = UniProtKB.AC, Gene.Name = Gene_Name,
                Protein.ID = UniProtKB.AC, Description = Protein_Name
              ) %>%
              dplyr::filter(!duplicated(UniqueID)), by = "UniqueID")
        } else {
          ProteinGeneName <- data.frame(id = 1:length(IDall), UniqueID = IDall) %>%
            dplyr::left_join(ProteinInfo %>%
              dplyr::transmute(
                UniqueID = UniProtKB.ID, Gene.Name = Gene_Name,
                Protein.ID = UniProtKB.AC, Description = Protein_Name
              ) %>%
              dplyr::filter(!duplicated(UniqueID)), by = "UniqueID")
        }

        if (!include_description) {
          ProteinGeneName <- ProteinGeneName %>% dplyr::select(-Description)
        }

        if (fill_missing_names) {
          ProteinGeneName <- ProteinGeneName %>%
            dplyr::mutate(Gene.Name = ifelse(is.na(Gene.Name), UniqueID, Gene.Name))
        }
      } else {
        stop("UniProt annotation requested but db/ProteinInfo.rds file not found")
      }
    } else {
      # Gene annotation using biomaRt
      if (species == "rat") {
        ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl")
      } else if (species == "mouse") {
        ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
      } else {
        ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
      }

      if (id_type == "Ensembl Gene ID") {
        filter_type <- "ensembl_gene_id"
        IDall_old <- IDall
        IDall <- str_replace(IDall, "\\.\\d+$", "")
        EID <- data.frame(IDall_old, IDall)
      } else if (id_type == "NCBI GeneID") {
        filter_type <- "entrezgene_id"
      } else if (id_type == "Gene Symbol") {
        filter_type <- "external_gene_name"
      } else {
        stop("Unsupported ID type for biomaRt annotation: ", id_type)
      }

      E_attributes <- c("ensembl_gene_id", "external_gene_name", "gene_biotype", "entrezgene_id")
      if (include_description) {
        E_attributes <- c(E_attributes, "description")
      }

      output <- biomaRt::getBM(
        attributes = E_attributes, filters = filter_type,
        values = IDall, mart = ensembl, useCache = FALSE
      )

      if (nrow(output) == 0) {
        stop("No gene annotation extracted from Biomart. Check Species and ID Type.")
      }

      output <- output %>% dplyr::arrange(entrezgene_id, ensembl_gene_id)

      if (!include_description) {
        output$description <- NA
      }

      F_TYPE <- sym(filter_type)

      ProteinGeneName <- data.frame(UniqueID = IDall) %>%
        dplyr::left_join(output %>%
          dplyr::transmute(
            UniqueID = !!F_TYPE, Gene.Name = external_gene_name,
            GeneType = gene_biotype, Description = description
          ) %>%
          dplyr::filter(!duplicated(UniqueID)), by = "UniqueID")

      if (id_type == "Ensembl Gene ID") {
        ProteinGeneName <- EID %>%
          dplyr::transmute(UniqueID = IDall_old, Unique1 = IDall) %>%
          left_join(ProteinGeneName %>% mutate(Unique1 = UniqueID) %>% dplyr::select(-UniqueID), by = "Unique1") %>%
          dplyr::select(-Unique1)
      }

      ProteinGeneName$Protein.ID <- NA
      ProteinGeneName <- ProteinGeneName %>%
        dplyr::select(UniqueID, Gene.Name, Protein.ID, GeneType, Description)

      if (!include_description) {
        ProteinGeneName <- ProteinGeneName %>% dplyr::select(-Description)
      }

      if (fill_missing_names) {
        ProteinGeneName <- ProteinGeneName %>%
          dplyr::mutate(Gene.Name = ifelse(is.na(Gene.Name), UniqueID, Gene.Name))
      }
    }
  } else {
    # Create minimal ProteinGeneName if no annotation
    ProteinGeneName <- data.frame(
      UniqueID = IDall,
      Gene.Name = IDall,
      Protein.ID = NA,
      stringsAsFactors = FALSE
    )

    if (fill_missing_names) {
      ProteinGeneName <- ProteinGeneName %>%
        dplyr::mutate(Gene.Name = ifelse(is.na(Gene.Name), UniqueID, Gene.Name))
    }
  }

  # Add id column if it doesn't exist
  if (!"id" %in% names(ProteinGeneName)) {
    ProteinGeneName$id <- 1:nrow(ProteinGeneName)
  }

  # Process groups and metadata
  if ("Order" %in% names(MetaData)) {
    groups <- group_order <- as.character(MetaData$Order[MetaData$Order != ""])
  } else {
    groups <- group_order <- MetaData %>%
      dplyr::pull(group) %>%
      unique()
  }

  samples <- sample_order <- as.character(MetaData$sampleid[order(match(MetaData$group, groups))])

  # Add comparison pairs and order to MetaData (required by the app)
  if (!"ComparePairs" %in% names(MetaData)) {
    MetaData$ComparePairs <- ""
    if (length(tests) > 0) {
      MetaData$ComparePairs[1:length(tests)] <- tests
    }
  }

  if (!"Order" %in% names(MetaData)) {
    MetaData$Order <- ""
    if (length(groups) > 0) {
      MetaData$Order[1:length(groups)] <- groups
    }
  }

  # Create MetaData_long
  MetaData_long <- MetaData %>%
    dplyr::select(-any_of(c("Order", "ComparePairs"))) %>%
    dplyr::mutate_all(as.character) %>%
    tidyr::pivot_longer(cols = -sampleid, names_to = "type", values_to = "group")

  # Create data_results
  data_results <- ProteinGeneName %>%
    dplyr::select(any_of(c("id", "UniqueID", "Gene.Name", "Protein.ID"))) %>%
    dplyr::left_join(
      data.frame(
        UniqueID = rownames(data_wide),
        Intensity = apply(data_wide, 1, mean, na.rm = TRUE)
      ) %>%
        dplyr::filter(!duplicated(UniqueID)),
      by = "UniqueID"
    )

  # Add group means and standard deviations
  sinfo1 <- data.frame(sampleid = names(data_wide)) %>%
    dplyr::left_join(MetaData %>% dplyr::select(sampleid, group), by = "sampleid")

  for (grp in unique(sinfo1$group)) {
    if (!is.na(grp)) {
      subdata <- data.frame(
        UniqueID = rownames(data_wide),
        t(apply(data_wide[, sinfo1$group == grp, drop = FALSE], 1, function(x) {
          return(setNames(
            c(mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)),
            paste(grp, c("Mean", "sd"), sep = "_")
          ))
        })),
        check.names = FALSE
      )
      data_results <- data_results %>%
        dplyr::left_join(subdata %>% dplyr::filter(!duplicated(UniqueID)), by = "UniqueID")
    }
  }

  # Add comparison results if available
  if (!is.null(results_long) && nrow(results_long) > 0) {
    for (ctr in tests) {
      subdata <- results_long %>%
        dplyr::filter(test == ctr) %>%
        dplyr::select(UniqueID, logFC, P.Value, Adj.P.Value)
      names(subdata)[2:4] <- stringr::str_c(ctr, "_", names(subdata)[2:4])
      data_results <- data_results %>%
        dplyr::left_join(subdata %>% dplyr::filter(!duplicated(UniqueID)), by = "UniqueID")
    }

    # Join ProteinGeneName with results_long
    results_long <- results_long %>%
      dplyr::mutate_if(is.factor, as.character) %>%
      dplyr::inner_join(ProteinGeneName, ., by = "UniqueID")
  }

  # Create network if requested (optional, can be time-consuming)
  network <- NULL
  file2 <- NULL
  if (create_network && nrow(data_wide) < 5000) {
    tryCatch(
      {
        cor_res <- Hmisc::rcorr(as.matrix(t(data_wide)))
        cormat <- cor_res$r
        pmat <- cor_res$P
        ut <- upper.tri(cormat)

        network <- tibble(
          from = rownames(cormat)[row(cormat)[ut]],
          to = rownames(cormat)[col(cormat)[ut]],
          cor = signif(cormat[ut], 2),
          p = signif(pmat[ut], 2),
          direction = as.integer(sign(cormat[ut]))
        ) %>%
          dplyr::filter(!is.na(cor) & abs(cor) > 0.7 & p < 0.05)

        # Limit network size
        if (nrow(network) > 2e6) {
          network <- network %>% dplyr::filter(abs(cor) > 0.8 & p < 0.005)
        }
        if (nrow(network) > 2e6) {
          network <- network %>% dplyr::filter(abs(cor) > 0.85 & p < 0.005)
        }

        cat("Final network size:", nrow(network), "\n")
      },
      error = function(e) {
        warning("Network creation failed: ", e$message)
      }
    )
  }

  # Generate unique ProjectID
  ProjectID <- make.names(project_name)
  if (nchar(ProjectID) > 45) {
    ProjectID <- substr(ProjectID, 1, 45)
  }
  ProjectID <- stringr::str_c(ProjectID, "_", stringi::stri_rand_strings(1, 6))

  # Save RData files
  file1 <- file.path(temp_dir, paste0(ProjectID, ".RData"))
  save(data_long, data_results, data_wide, MetaData, ProteinGeneName, results_long,
    file = file1
  )

  if (!is.null(network)) {
    file2 <- file.path(temp_dir, paste0(ProjectID, "_network.RData"))
    save(network, file = file2)
  }

  # Prepare return data structure (matching the expected format from inputdata.R)
  returnlist <- list(
    ProjectID = ProjectID,
    Name = project_name,
    Species = species,
    ShortName = project_name,
    Path = temp_dir,
    file1 = file1,
    file2 = file2,
    exp_unit = "Expression Level",
    MetaData = MetaData,
    MetaData_long = MetaData_long,
    ProteinGeneName = ProteinGeneName,
    ProteinGeneNameHeader = colnames(ProteinGeneName),
    data_long = data_long,
    data_wide = data_wide,
    results_long = results_long,
    data_results = data_results,
    groups = groups,
    group_order = group_order,
    samples = samples,
    sample_order = sample_order,
    tests = tests,
    tests_order = tests
  )

  return(returnlist)
}
