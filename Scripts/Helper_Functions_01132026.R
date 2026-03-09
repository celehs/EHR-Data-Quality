
# ----------------------- GLOBAL SETTINGS -----------------------

required_packages <- c(
  "dplyr", "stringr", "tidyr", "ggplot2", "gridExtra", 
  "ggtext", "knitr", "lubridate", "cowplot", "RColorBrewer", "purrr"
)

# Install any packages that are missing
installed <- rownames(installed.packages())
to_install <- setdiff(required_packages, installed)

if (length(to_install) > 0) {
  install.packages(to_install)
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Set global themes
theme_global <- theme_minimal() + 
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_textbox_simple(size = 16, lineheight = 1.3),  
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

update_geom_defaults("line", list(linewidth = 1))
update_geom_defaults("point", list(size = 3))

theme_bar <- theme_minimal() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_textbox_simple(size = 16, lineheight = 1.3),  
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

set_palette <- scale_color_brewer(palette = "Set2")
set_fill_palette <- scale_fill_brewer(palette = "Set2")

# ----------------------- MODULE 1 HELPERS -----------------------

clean_ONCE_data <- function(target_code, O2 = TRUE, path_code = NULL, path_nlp = NULL) {
  formatted_code <- gsub(":", "", target_code)
  formatted_code <- gsub("\\.", "_", formatted_code)
  
  if (O2) {
    path_code <- paste0("/n/data1/hsph/biostat/celehs/lab/sm731/ONCE_features/codified/", formatted_code, "_cod_features.csv")
    path_nlp <- paste0("/n/data1/hsph/biostat/celehs/lab/sm731/ONCE_features/NLP/", formatted_code, "_nlp_features.csv")
  } else {
    if (is.null(path_code) || is.null(path_nlp)) {
      stop("Must provide manual ONCE paths when O2 is FALSE.")
    }
  }
  
  code <- read.csv(path_code) %>%
    dplyr::select(Variable, Description, target_similarity) %>%
    rename(feature_id = Variable, description = Description) %>%
    mutate(description = tolower(description))
  
  nlp <- read.csv(path_nlp) %>%
    dplyr::select(cui, term, target_similarity) %>%
    rename(feature_id = cui, description = term) %>%
    mutate(description = tolower(description))
  
  return(list(code = code, nlp = nlp))
}

clean_raw_data <- function(df, date_col = "start_date", id_col = "feature_id") {
  df <- df %>%
    mutate(start_date_chr = as.character(.data[[date_col]]))
  
  bad_dates <- df %>%
    filter(!stringr::str_detect(start_date_chr, "^\\d{4}$|^\\d{4}-\\d{2}-\\d{2}$"))
  
  if (nrow(bad_dates) > 0) {
    warning(nrow(bad_dates), " rows have invalid date format in `", date_col, "` and were set to NA.")
  }
  
  df %>%
    mutate(
      patient_num = as.character(patient_num),
      year = dplyr::case_when(
        stringr::str_detect(start_date_chr, "^\\d{4}$") ~ start_date_chr,
        stringr::str_detect(start_date_chr, "^\\d{4}-\\d{2}-\\d{2}$") ~ as.character(
          lubridate::year(suppressWarnings(lubridate::ymd(start_date_chr)))
        ),
        TRUE ~ NA_character_
      ),
      feature_id_clean = stringr::str_replace(.data[[id_col]], "CCS-PCS", "CCS")
    ) %>%
    dplyr::select(patient_num, year, feature_id = feature_id_clean)
}

clean_labels <- function(df) {
  req <- c("patient_num", "label")
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop("labels file missing required columns: ", paste(miss, collapse = ", "))
  
  df %>%
    dplyr::transmute(
      patient_num = as.character(patient_num),
      label = as.integer(label)
    )
}

clean_data <- function(paths, date_col = "start_date", id_col = "feature_id") {
  if (is.null(names(paths)) || any(names(paths) == "")) {
    stop("`paths` must be a *named* list (e.g., paths$codified1, paths$nlp1, optional paths$labels).")
  }
  
  # Read everything except optional labels
  paths_main <- paths[names(paths) != "labels"]
  input_main <- lapply(paths_main, read.csv, stringsAsFactors = FALSE)
  cleaned_main <- lapply(input_main, clean_raw_data, date_col = date_col, id_col = id_col)
  names(cleaned_main) <- names(paths_main)
  
  # If labels path is provided AND exists, read + clean and attach
  if (!is.null(paths$labels) && is.character(paths$labels) && length(paths$labels) == 1 &&
      nzchar(paths$labels) && file.exists(paths$labels)) {
    labels_df <- read.csv(paths$labels, stringsAsFactors = FALSE)
    cleaned_main$labels <- clean_labels(labels_df)
  }
  
  return(cleaned_main)
}

clean_dictionary <- function(dictionary_path) {
  read.csv(dictionary_path) %>%
    select(feature_id, description) %>%
    distinct() %>%
    mutate(
      description = case_when(
        str_detect(description, "forms of") ~ str_replace(description, "(.*)forms of.*", "\\1forms"),
        TRUE ~ str_replace_all(description, "\\s*\\(.*?\\)|,.*", "")
      ),
      description = tolower(description)
    )
}

# ----------------------- MODULE 2 HELPERS -----------------------

generate_all_patient_data <- function(data_inputs_subset, target_feature, sample_labels) {
  counts_list <- list()
  summary_list <- list()
  
  for (name in names(data_inputs_subset)) {
    sample_key <- ifelse(grepl("2", name), "2", "1")
    label <- sample_labels[[sample_key]]
    df <- data_inputs_subset[[name]]
    
    counts_list[[label]] <- generate_patient_counts(df, target_feature)
    summary_list[[label]] <- generate_total_patient_summary(df, target_feature, label)
  }
  
  list(counts = counts_list, summary = summary_list)
}

# 2.a. Patient counts per year for one dataset
generate_patient_counts <- function(data, target_feature) {
  total_counts <- data %>%
    group_by(year) %>%
    summarise(total_patients = n_distinct(patient_num), .groups = "drop")
  
  target_counts <- data %>%
    filter(feature_id %in% target_feature) %>%      
    group_by(year) %>%
    summarise(target_patients = n_distinct(patient_num), .groups = "drop")
  
  total_counts %>%
    left_join(target_counts, by = "year") %>%
    mutate(target_patients = replace_na(target_patients, 0)) %>%  
    rename(Year = year, `Total` = total_patients, `Target` = target_patients) %>%
    pivot_longer(cols = c("Total", "Target"), names_to = "Type", values_to = "Patients")
}

# 2.b. Plot patients by year (single or dual samples)
generate_line_plot <- function(counts_list, title) {
  combined_data <- bind_rows(
    Map(function(df, label) mutate(df, Sample = label),
        counts_list,
        names(counts_list))
  )
  
  ggplot(combined_data, aes(
    x = Year,
    y = Patients,
    color = Type,
    group = interaction(Type, Sample),
    linetype = Sample
  )) +
    geom_line() +
    geom_point() +
    scale_linetype_manual(values = setNames(c("solid", "dashed")[seq_along(counts_list)], names(counts_list))) +
    labs(
      title = title,
      x = "",
      y = "Patients per Year",
      color = "Type",
      linetype = "Sample"
    ) +
    set_palette +
    theme_global
}

generate_line_plot <- function(counts_list, title) {
  combined_data <- bind_rows(
    Map(function(df, label) mutate(df, Sample = label), counts_list, names(counts_list))
  )
  
  ggplot(combined_data, aes(
    x = Year,
    y = Patients,
    color = Type,
    group = interaction(Type, Sample),
    linetype = Sample
  )) +
    geom_line() +
    geom_point() +
    scale_linetype_manual(values = setNames(c("solid", "dashed"), names(counts_list))) +  # Dynamic line types
    labs(
      title = title,
      x = "",
      y = "Patients per Year",
      color = "Type",
      linetype = "Sample"
    ) +
    set_palette +
    theme_global
}

# 2.c. Total patient summary
generate_total_patient_summary <- function(data, target_feature, label) {
  data.frame(
    Feature = c("Total Patients", "Patients with Target"),
    Total_Patients = c(
      n_distinct(data$patient_num),
      n_distinct(data$patient_num[data$feature_id %in% target_feature]) 
    ),
    Dataset = label
  )
}

# 2.d. Stacked bar plot for 1 or 2 datasets
generate_bar_plot <- function(summary_list, title) {
  combined_data <- bind_rows(summary_list)
  
  labels <- names(summary_list)
  first_label <- labels[1]
  second_label <- labels[2]
  
  ggplot(combined_data, aes(x = Feature, y = Total_Patients, fill = Feature)) +
    geom_bar(data = combined_data %>% filter(Dataset == second_label),
             stat = "identity", 
             position = position_identity(),
             width = 0.4,
             aes(x = as.numeric(as.factor(Feature)) + 0.1),
             color = "black",
             linetype = "dashed") +
    
    geom_bar(data = combined_data %>% filter(Dataset == first_label),
             stat = "identity", 
             position = position_identity(),
             width = 0.4,
             aes(x = as.numeric(as.factor(Feature)) - 0.1),
             color = "black") +
    
    scale_x_continuous(breaks = unique(as.numeric(as.factor(combined_data$Feature))), 
                       labels = unique(combined_data$Feature)) +
    labs(title = title, y = "Total Patients Across All Years", x = "Feature") +
    set_fill_palette +
    theme_bar
}

# Final wrapper function (module 2)
plot_target_prevalence <- function(
    data_inputs, target_code, target_cui, sample_labels,
    target_code_label = NULL, target_cui_label = NULL
) {
  
  # Labels for titles/subtitles
  if (is.null(target_code_label)) target_code_label <- paste(target_code, collapse = " + ")
  if (is.null(target_cui_label))  target_cui_label  <- paste(target_cui,  collapse = " + ")
  
  sample_sizes <- lapply(names(data_inputs), function(name) {
    data <- data_inputs[[name]]
    data_type <- ifelse(grepl("nlp", name), "NLP", "Codified")
    sample_key <- ifelse(grepl("2", name), "2", "1")
    sample_label <- sample_labels[[sample_key]]
    data.frame(
      Sample = sample_label,
      Dataset = data_type,
      `Number of Patients` = dplyr::n_distinct(data$patient_num)
    )
  }) %>% dplyr::bind_rows()
  
  # Helper: bar plot that works for 1 or 2 datasets
  generate_bar_plot_safe <- function(summary_list, title) {
    combined_data <- dplyr::bind_rows(summary_list)
    dataset_names <- unique(combined_data$Dataset)
    
    if (length(dataset_names) == 1) {
      ggplot2::ggplot(combined_data, ggplot2::aes(x = Feature, y = Total_Patients, fill = Feature)) +
        ggplot2::geom_bar(stat = "identity", width = 0.6, color = "black") +
        ggplot2::labs(title = title, y = "Total Patients Across All Years", x = "Feature") +
        set_fill_palette +
        theme_bar
    } else if (length(dataset_names) == 2) {
      
      first_label <- dataset_names[1]
      second_label <- dataset_names[2]
      
      ggplot2::ggplot(combined_data, ggplot2::aes(x = Feature, y = Total_Patients, fill = Feature)) +
        ggplot2::geom_bar(
          data = combined_data %>% dplyr::filter(Dataset == second_label),
          stat = "identity",
          position = ggplot2::position_identity(),
          width = 0.4,
          ggplot2::aes(x = as.numeric(as.factor(Feature)) + 0.1),
          color = "black",
          linetype = "dashed"
        ) +
        ggplot2::geom_bar(
          data = combined_data %>% dplyr::filter(Dataset == first_label),
          stat = "identity",
          position = ggplot2::position_identity(),
          width = 0.4,
          ggplot2::aes(x = as.numeric(as.factor(Feature)) - 0.1),
          color = "black"
        ) +
        ggplot2::scale_x_continuous(
          breaks = unique(as.numeric(as.factor(combined_data$Feature))),
          labels = unique(combined_data$Feature)
        ) +
        ggplot2::labs(title = title, y = "Total Patients Across All Years", x = "Feature") +
        set_fill_palette +
        theme_bar
    } else {
      stop("This pipeline supports only 1 or 2 samples.")
    }
  }
  
  # ---- NLP ----
  nlp_data_inputs <- data_inputs[grepl("nlp", names(data_inputs))]
  nlp_summaries <- generate_all_patient_data(nlp_data_inputs, target_cui, sample_labels)
  
  plot_nlp <- generate_line_plot(
    nlp_summaries$counts,
    paste0("Sample Size by Year (NLP)\nTarget: ", target_cui_label)
  )
  plot_bar_nlp <- generate_bar_plot_safe(
    nlp_summaries$summary,
    paste0("Total Patients vs. Target CUI\nTarget: ", target_cui_label)
  )
  
  # ---- Codified ----
  cod_data_inputs <- data_inputs[grepl("codified", names(data_inputs))]
  cod_summaries <- generate_all_patient_data(cod_data_inputs, target_code, sample_labels)
  
  plot_codified <- generate_line_plot(
    cod_summaries$counts,
    paste0("Sample Size by Year (Codified)\nTarget: ", target_code_label)
  )
  plot_bar_codified <- generate_bar_plot_safe(
    cod_summaries$summary,
    paste0("Total Patients vs. Target Code\nTarget: ", target_code_label)
  )
  
  # Legends
  legend_nlp <- suppressWarnings(cowplot::get_legend(plot_nlp + ggplot2::theme(legend.position = "right")))
  legend_codified <- suppressWarnings(cowplot::get_legend(plot_codified + ggplot2::theme(legend.position = "right")))
  
  combined_plot_nlp <- cowplot::plot_grid(
    plot_nlp + ggplot2::theme(legend.position = "none"),
    plot_bar_nlp,
    legend_nlp,
    ncol = 3,
    rel_widths = c(3, 2, 1)
  )
  
  combined_plot_codified <- cowplot::plot_grid(
    plot_codified + ggplot2::theme(legend.position = "none"),
    plot_bar_codified,
    legend_codified,
    ncol = 3,
    rel_widths = c(3, 2, 1)
  )
  
  return(list(
    sample_sizes = sample_sizes,
    nlp_plot = combined_plot_nlp,
    codified_plot = combined_plot_codified
  ))
}



plot_trends_v2 <- function(data_list, phecode_pattern, dictionary, y_var, y_label, title_suffix, sample_labels, custom_children = NULL) {
  phecode_description <- get_phecode_description(phecode_pattern, dictionary)
  
  if (!is.null(custom_children) && phecode_pattern %in% names(custom_children)) {
    included_codes <- c(phecode_pattern, custom_children[[phecode_pattern]])
  } else {
    included_codes <- dictionary %>%
      filter(
        str_detect(feature_id, phecode_pattern) &
          !str_detect(feature_id, "\\.\\d{2,}")
      ) %>%
      pull(feature_id)
  }
  
  code_info <- get_ordered_code_info(included_codes, dictionary)
  
  combined <- purrr::imap_dfr(data_list, function(df, key) {
    df %>% mutate(Sample = sample_labels[[key]])
  }) %>%
    filter(feature_id %in% code_info$ordered_codes) %>%
    mutate(
      is_parent = ifelse(feature_id == phecode_pattern, "Parent", "Child"),
      feature_id = factor(feature_id, levels = code_info$ordered_codes)
    )
  
  ggplot() +
    geom_line(data = combined,
              aes(x = year, y = !!sym(y_var), color = feature_id,
                  group = interaction(feature_id, Sample),
                  linetype = Sample,
                  size = is_parent)) +
    geom_point(data = combined,
               aes(x = year, y = !!sym(y_var), color = feature_id,
                   group = interaction(feature_id, Sample))) +
    scale_linetype_manual(values = setNames(c("solid", "dashed")[seq_along(data_list)], unname(sample_labels[seq_along(data_list)]))) +
    scale_color_manual(values = code_info$color_mapping, labels = code_info$legend_labels, name = "Code") +
    scale_size_manual(values = c("Parent" = 2, "Child" = 1), guide = "none") +
    labs(
      title = paste(title_suffix, "for Parent-Child Phecodes"),
      subtitle = paste0("<b><i>", phecode_pattern, "</i></b> | <i>", phecode_description, "</i>"),
      x = "",
      y = y_label
    ) +
    theme_global
}

plot_total_patients_v2 <- function(data_list, phecode_pattern, dictionary, sample_labels, custom_children = NULL) {
  # Step 1: Filter custom_children to only codes present in the data
  all_feature_ids <- unique(unlist(lapply(data_list, function(df) unique(df$feature_id))))
  if (!is.null(custom_children) && phecode_pattern %in% names(custom_children)) {
    valid_children <- intersect(custom_children[[phecode_pattern]], all_feature_ids)
    included_codes <- c(phecode_pattern, valid_children)
  } else {
    included_codes <- dictionary %>%
      filter(
        str_detect(feature_id, phecode_pattern) &
          !str_detect(feature_id, "\\.\\d{2,}")
      ) %>%
      pull(feature_id)
  }
  
  code_info <- get_ordered_code_info(included_codes, dictionary)
  
  # Step 2: Summarize total patients per code across all years
  all_data <- purrr::imap_dfr(data_list, function(data, key) {
    label <- sample_labels[[key]]
    data %>%
      filter(feature_id %in% code_info$ordered_codes) %>%
      group_by(feature_id) %>%
      summarise(Total_Patients = n_distinct(patient_num), .groups = "drop") %>%
      left_join(dictionary, by = "feature_id") %>%
      mutate(Dataset = label)
  })
  
  # Step 3: Clean up factor levels and axis setup
  codes_present <- unique(all_data$feature_id)
  ordered_codes <- code_info$ordered_codes[code_info$ordered_codes %in% codes_present]
  all_data$feature_id <- factor(all_data$feature_id, levels = ordered_codes)
  all_data <- all_data %>% mutate(x_pos = as.numeric(feature_id))
  x_levels <- levels(all_data$feature_id)
  x_breaks <- seq_along(x_levels)
  x_labels <- unname(code_info$legend_labels[x_levels])
  color_mapping <- code_info$color_mapping[x_levels]
  
  # Step 4: Get dataset names
  dataset_names <- unique(all_data$Dataset)
  if (length(dataset_names) != 2) {
    stop("plot_total_patients_v2 requires exactly 2 datasets.")
  }
  
  # Step 5: Plot
  ggplot(all_data, aes(x = x_pos, y = Total_Patients, fill = feature_id)) +
    geom_bar(data = all_data %>% filter(Dataset == dataset_names[2]),
             stat = "identity", position = position_identity(),
             width = 0.4, aes(x = x_pos + 0.1),
             linetype = "dashed", color = "black") +
    geom_bar(data = all_data %>% filter(Dataset == dataset_names[1]),
             stat = "identity", position = position_identity(),
             width = 0.4, aes(x = x_pos - 0.1),
             color = "black") +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels
    ) +
    scale_fill_manual(
      values = color_mapping,
      labels = x_labels,
      name = "Code"
    ) +
    labs(x = "", y = "Total Patients Across All Years") +
    theme_bar
}




plot_total_patients_v2 <- function(data_list, phecode_pattern, dictionary, sample_labels, custom_children = NULL) {
  # Step 1: Determine which codes to include
  all_feature_ids <- unique(unlist(lapply(data_list, \(df) unique(df$feature_id))))
  if (!is.null(custom_children) && phecode_pattern %in% names(custom_children)) {
    valid_children <- intersect(custom_children[[phecode_pattern]], all_feature_ids)
    included_codes <- c(phecode_pattern, valid_children)
  } else {
    included_codes <- dictionary %>%
      filter(
        str_detect(feature_id, phecode_pattern) &
          !str_detect(feature_id, "\\.\\d{2,}")
      ) %>%
      pull(feature_id)
  }
  
  code_info <- get_ordered_code_info(included_codes, dictionary)
  
  # Step 2: Summarize total patients per code across all years
  all_data <- purrr::imap_dfr(data_list, function(data, key) {
    label <- sample_labels[[key]]
    data %>%
      filter(feature_id %in% code_info$ordered_codes) %>%
      group_by(feature_id) %>%
      summarise(Total_Patients = n_distinct(patient_num), .groups = "drop") %>%
      left_join(dictionary, by = "feature_id") %>%
      mutate(Dataset = label)
  })
  
  # Step 3: Set up plot variables
  codes_present <- unique(all_data$feature_id)
  ordered_codes <- code_info$ordered_codes[code_info$ordered_codes %in% codes_present]
  all_data$feature_id <- factor(all_data$feature_id, levels = ordered_codes)
  all_data <- all_data %>% mutate(x_pos = as.numeric(feature_id))
  x_levels <- levels(all_data$feature_id)
  x_breaks <- seq_along(x_levels)
  x_labels <- unname(code_info$legend_labels[x_levels])
  color_mapping <- code_info$color_mapping[x_levels]
  
  # Step 4: Plot depending on number of datasets
  dataset_names <- unique(all_data$Dataset)
  
  if (length(dataset_names) == 1) {
    # One dataset — simple bar plot
    ggplot(all_data, aes(x = feature_id, y = Total_Patients, fill = feature_id)) +
      geom_bar(stat = "identity", position = position_dodge(), color = "black", width = 0.6) +
      scale_fill_manual(
        values = color_mapping,
        labels = x_labels,
        name = "Code"
      ) +
      labs(x = "", y = "Total Patients Across All Years") +
      theme_bar
  } else if (length(dataset_names) == 2) {
    # Two datasets — offset bars with dashed outline
    ggplot(all_data, aes(x = x_pos, y = Total_Patients, fill = feature_id)) +
      geom_bar(data = all_data %>% filter(Dataset == dataset_names[2]),
               stat = "identity", position = position_identity(),
               width = 0.4, aes(x = x_pos + 0.1),
               linetype = "dashed", color = "black") +
      geom_bar(data = all_data %>% filter(Dataset == dataset_names[1]),
               stat = "identity", position = position_identity(),
               width = 0.4, aes(x = x_pos - 0.1),
               color = "black") +
      scale_x_continuous(
        breaks = x_breaks,
        labels = x_labels
      ) +
      scale_fill_manual(
        values = color_mapping,
        labels = x_labels,
        name = "Code"
      ) +
      labs(x = "", y = "Total Patients Across All Years") +
      theme_bar
  } else {
    stop("plot_total_patients_v2 supports only 1 or 2 datasets.")
  }
}

generate_phecode_summary <- function(data) {
  total_patients <- data %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(Total_Patients = dplyr::n_distinct(patient_num), .groups = "drop")
  
  data %>%
    dplyr::group_by(year, feature_id) %>%
    dplyr::summarise(Count = dplyr::n(), Patients = dplyr::n_distinct(patient_num), .groups = "drop") %>%
    dplyr::left_join(total_patients, by = "year") %>%
    dplyr::mutate(Rate = Patients / Total_Patients)
}

get_phecode_description <- function(phecode_pattern, dictionary) {
  description <- dictionary %>%
    dplyr::filter(feature_id == phecode_pattern) %>%
    dplyr::pull(description)
  ifelse(length(description) == 0, "description not found", description)
}

get_ordered_code_info <- function(included_codes, dictionary) {
  ordered_codes <- included_codes[order(as.numeric(stringr::str_replace_all(included_codes, "[^0-9.]", "")))]
  
  code_info <- dictionary %>%
    dplyr::filter(feature_id %in% ordered_codes) %>%
    dplyr::mutate(feature_label = paste0(feature_id, " | ", description)) %>%
    dplyr::arrange(factor(feature_id, levels = ordered_codes))
  
  labels <- code_info$feature_label
  names(labels) <- code_info$feature_id
  
  colors <- setNames(RColorBrewer::brewer.pal(n = length(ordered_codes), name = "Set2"), ordered_codes)
  
  list(
    ordered_codes = ordered_codes,
    legend_labels = labels,
    color_mapping = colors
  )
}

# Final wrapper function (module 3)
analyze_code_hierarchy <- function(data_inputs, dictionary, sample_labels, phecodes, custom_children = NULL) {
  # Prepare codified data list using sample_labels
  codified_data_list <- list()
  if (!is.null(data_inputs$codified1)) codified_data_list[["1"]] <- data_inputs$codified1
  if (!is.null(data_inputs$codified2)) codified_data_list[["2"]] <- data_inputs$codified2
  
  # Generate Phecode summaries
  phecode_data <- lapply(codified_data_list, generate_phecode_summary)
  
  # Loop through Phecodes and generate plots
  plots <- lapply(phecodes, function(p) {
    rate_plot <- plot_trends_v2(phecode_data, p, dictionary, "Rate", "Rate", "Rates", sample_labels, custom_children)
    count_plot <- plot_trends_v2(phecode_data, p, dictionary, "Patients", "Patients per Year", "Patient Counts", sample_labels, custom_children)
    bar_plot <- plot_total_patients_v2(codified_data_list, p, dictionary, sample_labels, custom_children)
    
    combined_plot <- cowplot::plot_grid(
      count_plot + theme(legend.position = "none"),
      bar_plot + theme(legend.position = "none"),
      ncol = 2,
      rel_widths = c(3, 2)
    )
    
    list(phecode = p, rate_plot = rate_plot, combined_plot = combined_plot)
  })
  
  return(plots)
}

# ----------------------- MODULE 4 HELPERS -----------------------

collapse_one <- function(x, sep = " + ") {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("")
  paste(unique(as.character(x)), collapse = sep)
}

safe_desc <- function(dictionary, ids) {
  if (is.null(dictionary) || is.null(ids)) return("description not found")
  d <- dictionary %>% dplyr::filter(feature_id %in% ids) %>% dplyr::pull(description)
  d <- unique(d[!is.na(d)])
  if (length(d) == 0) "description not found" else collapse_one(d, sep = "; ")
}

summarize_target_code_trends <- function(codified_list, nlp_list, target_code, target_cui, sample_labels,
                                         target_code_label = "Target Code", target_cui_label = "Target CUI") {
  
  keys <- names(codified_list)
  
  all_summaries <- map2(seq_along(codified_list), seq_along(nlp_list), function(i, j) {
    codified <- codified_list[[i]]
    nlp <- nlp_list[[j]]
    key <- keys[[i]]
    label <- sample_labels[[key]]
    
    # --- patient-year “composite counts” for correlation ---
    wide_base <- bind_rows(
      codified %>% mutate(src = "codified"),
      nlp      %>% mutate(src = "nlp")
    ) %>%
      group_by(year, patient_num) %>%
      summarise(
        code_count = sum(feature_id %in% target_code),   # counts of target-code rows (can be 0+)
        cui_count  = sum(feature_id %in% target_cui),    # counts of target-cui rows (can be 0+)
        .groups = "drop"
      ) %>%
      mutate(Sample = label)
    
    # --- yearly patient counts & rates (ANY-of targets) ---
    target_phecode_summary <- codified %>%
      filter(feature_id %in% target_code) %>%
      group_by(year) %>%
      summarise(PheCode_Patients = n_distinct(patient_num), .groups = "drop")
    
    target_cui_summary <- nlp %>%
      filter(feature_id %in% target_cui) %>%
      group_by(year) %>%
      summarise(CUI_Patients = n_distinct(patient_num), .groups = "drop")
    
    target_total_patients <- bind_rows(codified, nlp) %>%
      group_by(year) %>%
      summarise(Total_Patients = n_distinct(patient_num), .groups = "drop")
    
    combined <- target_phecode_summary %>%
      full_join(target_cui_summary, by = "year") %>%
      left_join(target_total_patients, by = "year") %>%
      mutate(
        PheCode_Patients = replace_na(PheCode_Patients, 0),
        CUI_Patients     = replace_na(CUI_Patients, 0),
        PheCode_Rate = PheCode_Patients / Total_Patients,
        CUI_Rate     = CUI_Patients / Total_Patients,
        Sample = label
      )
    
    list(combined = combined, wide = wide_base)
  })
  
  list(
    combined = bind_rows(map(all_summaries, "combined")),
    wide     = bind_rows(map(all_summaries, "wide"))
  )
}

plot_target_code_trends <- function(summary_data, target_code, target_cui, dictionary_code, dictionary_nlp,
                                    target_code_label = NULL, target_cui_label = NULL) {
  
  # labels (single strings)
  if (is.null(target_code_label)) target_code_label <- collapse_one(target_code)
  if (is.null(target_cui_label))  target_cui_label  <- collapse_one(target_cui)
  
  # descriptions (single strings)
  phecode_description <- safe_desc(dictionary_code, target_code)
  cui_description     <- safe_desc(dictionary_nlp,  target_cui)
  
  subtitle_text <- paste0(
    "<i><b>Target Code(s):</b> ", target_code_label, "</i> | <i>", phecode_description, "</i>",
    "<br><b>Target CUI(s):</b> ", target_cui_label, "</i> | <i>", cui_description, "</i>"
    
  )
  # IMPORTANT: force scalar
  subtitle_text <- collapse_one(subtitle_text, sep = "")
  
  plot_rates <- ggplot(summary_data, aes(x = as.numeric(year))) +
    geom_line(aes(y = PheCode_Rate, color = "Target Code", linetype = Sample)) +
    geom_line(aes(y = CUI_Rate, color = "Target CUI", linetype = Sample)) +
    geom_point(aes(y = PheCode_Rate, color = "Target Code"), shape = 16) +
    geom_point(aes(y = CUI_Rate, color = "Target CUI"), shape = 16) +
    labs(
      title = "Rates for Target Code and CUI",
      subtitle = subtitle_text,
      x = "", y = "Rate", color = "Code", linetype = "Sample"
    ) +
    theme_global + set_palette +
    scale_x_continuous(breaks = sort(unique(as.numeric(summary_data$year))))
  
  plot_counts <- ggplot(summary_data, aes(x = as.numeric(year))) +
    geom_line(aes(y = PheCode_Patients, color = "Target Code", linetype = Sample)) +
    geom_line(aes(y = CUI_Patients, color = "Target CUI", linetype = Sample)) +
    geom_point(aes(y = PheCode_Patients, color = "Target Code"), shape = 16) +
    geom_point(aes(y = CUI_Patients, color = "Target CUI"), shape = 16) +
    labs(
      title = "Patient Counts for Target Code and CUI",
      subtitle = subtitle_text,
      x = "", y = "Patients per Year", color = "Code", linetype = "Sample"
    ) +
    theme_global + set_palette +
    scale_x_continuous(breaks = sort(unique(as.numeric(summary_data$year))))
  
  list(rates = plot_rates, counts = plot_counts)
}


plot_target_code_correlation <- function(wide_data, summary_data,
                                         target_code, target_cui,
                                         target_code_label = NULL, target_cui_label = NULL) {
  
  # Make scalar labels
  if (is.null(target_code_label)) target_code_label <- paste(unique(target_code), collapse = " + ")
  if (is.null(target_cui_label))  target_cui_label  <- paste(unique(target_cui),  collapse = " + ")
  
  subtitle_text <- paste0(
    "<i><b>Target Code(s):</b> ", target_code_label,
    "<br><b>Target CUI(s):</b> ", target_cui_label, "</i>"
  )
  subtitle_text <- paste0(subtitle_text, collapse = "")  # force scalar
  
  stopifnot(all(c("year", "Sample", "code_count", "cui_count") %in% names(wide_data)))
  
  correlation_by_year <- wide_data %>%
    dplyr::group_by(Sample, year) %>%
    dplyr::summarise(
      Correlation = {
        if (dplyr::n() < 5 || stats::sd(code_count, na.rm = TRUE) == 0 || stats::sd(cui_count, na.rm = TRUE) == 0) {
          NA_real_
        } else {
          stats::cor(code_count, cui_count, use = "pairwise.complete.obs", method = "spearman")
        }
      },
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(Correlation))
  
  ggplot2::ggplot(correlation_by_year, ggplot2::aes(x = as.numeric(year), y = Correlation, group = Sample)) +
    ggplot2::geom_line(ggplot2::aes(linetype = Sample), color = "black") +
    ggplot2::geom_point(color = "black", shape = 16) +
    ggplot2::labs(
      title = "Intra-Patient Correlation (Composite Target Code vs Composite Target CUI)",
      subtitle = subtitle_text,
      x = "", y = "Correlation", linetype = "Sample"
    ) +
    theme_global +
    ggplot2::scale_x_continuous(breaks = sort(unique(as.numeric(summary_data$year))))
}


# Final wrapper function (module 4)
code_cui_alignment <- function(data_inputs, target_code, target_cui, dictionary, sample_labels) {
  
  # Build codified and NLP lists from sample labels
  codified_list <- list()
  nlp_list <- list()
  if (!is.null(data_inputs$codified1)) codified_list[["1"]] <- data_inputs$codified1
  if (!is.null(data_inputs$codified2)) codified_list[["2"]] <- data_inputs$codified2
  if (!is.null(data_inputs$nlp1)) nlp_list[["1"]] <- data_inputs$nlp1
  if (!is.null(data_inputs$nlp2)) nlp_list[["2"]] <- data_inputs$nlp2
  
  # Summarize across codified and NLP
  summarized <- summarize_target_code_trends(
    codified_list = codified_list,
    nlp_list = nlp_list,
    target_code = target_code,
    target_cui = target_cui,
    sample_labels = sample_labels 
  )
  
  # Plot trends across samples/institutions
  trend_plots <- plot_target_code_trends(
    summary_data = summarized$combined,
    target_code = target_code,
    target_cui = target_cui,
    dictionary_code = dictionary,
    dictionary_nlp = dictionary
  )
  
  # Plot correlation across samples/institutions
  correlation_plot <- plot_target_code_correlation(
    wide_data   = summarized$wide,
    summary_data = summarized$combined,
    target_code = target_code,
    target_cui  = target_cui
  )
  
  return(list(
    rates_plot = trend_plots$rates,
    counts_plot = trend_plots$counts,
    correlation_plot = correlation_plot
  ))
}

# ----------------------- MODULE 5 HELPERS -----------------------

get_similarity_ordered_colors <- function(descriptions_df) {
  desc_ordered <- descriptions_df %>%
    arrange(desc(target_similarity)) %>%
    mutate(description_label = paste0(feature_id, " | ", description, " [", round(target_similarity, 3), "]"))
  
  labels <- desc_ordered$description_label
  names(labels) <- desc_ordered$feature_id
  
  colors <- setNames(RColorBrewer::brewer.pal(n = length(labels), name = "Set2"), labels)
  
  list(
    labels = labels,
    color_map = colors
  )
}

plot_related_code_trends_v2 <- function(Type, target_code, target_cui, codified_list, nlp_list,
                                        sample_labels, ONCE,
                                        type_dict = list("Diagnosis" = "PheCode",
                                                         "Medication" = "RXNORM",
                                                         "Lab" = c("LOINC","ShortName","Other Lab"),
                                                         "Procedure" = "CCS")) {
  
  if (!(Type %in% c(names(type_dict), "CUI"))) {
    return(list(NULL, paste("Invalid Type:", Type)))
  }
  
  is_cui <- (Type == "CUI")
  source_list <- if (is_cui) nlp_list else codified_list
  dict <- if (is_cui) ONCE$nlp else ONCE$code
  keyword <- if (!is_cui) type_dict[[Type]] else NULL
  
  # Vector targets are allowed
  target_features <- if (is_cui) target_cui else target_code
  target_features <- unique(as.character(target_features))
  target_label <- paste(target_features, collapse = " + ")
  
  # Keep only targets that actually appear in the data (prevents empty selections)
  present_features <- unique(unlist(lapply(source_list, function(df) df$feature_id)))
  target_present <- intersect(target_features, present_features)
  
  if (length(target_present) == 0) {
    return(list(NULL, paste("None of the target features are present in the", Type, "data.")))
  }
  
  # Related features are ranked from ONCE dict (dict is based on ONE anchor concept)
  # Exclude any target features
  related_features <- dict %>%
    dplyr::filter((if (is_cui) TRUE else stringr::str_detect(feature_id, paste0(keyword, collapse = "|")))) %>%
    dplyr::filter(!(feature_id %in% target_present)) %>%
    dplyr::filter(feature_id %in% present_features) %>%
    dplyr::arrange(dplyr::desc(target_similarity)) %>%
    dplyr::slice_head(n = 5)
  
  # Selected = targets present + related
  selected_features <- unique(c(target_present, related_features$feature_id))
  
  if (length(selected_features) == 0) {
    return(list(NULL, paste("No features found to plot for", Type)))
  }
  
  descriptions <- dict %>%
    dplyr::filter(feature_id %in% selected_features) %>%
    dplyr::distinct(feature_id, description, target_similarity) %>%
    dplyr::mutate(
      description = dplyr::if_else(is.na(description), "description not found", tolower(description)),
      target_similarity = dplyr::if_else(is.na(target_similarity), 0, target_similarity),
      description_label = paste0(feature_id, " | ", description, " [", round(target_similarity, 3), "]")
    )
  
  # Patient-year counts
  long_data <- purrr::imap_dfr(source_list, function(df, key) {
    df %>%
      dplyr::filter(feature_id %in% selected_features) %>%
      dplyr::group_by(year, feature_id) %>%
      dplyr::summarise(Patients = dplyr::n_distinct(patient_num), .groups = "drop") %>%
      dplyr::mutate(Sample = sample_labels[[key]])
  })
  
  # Total patients per year (denominator) — use codified_list (matches your current approach)
  total_data <- purrr::imap_dfr(codified_list, function(df, key) {
    df %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(Total_Patients = dplyr::n_distinct(patient_num), .groups = "drop") %>%
      dplyr::mutate(Sample = sample_labels[[key]])
  })
  
  combined_data <- long_data %>%
    dplyr::left_join(total_data, by = c("year", "Sample")) %>%
    dplyr::mutate(Rate = Patients / Total_Patients) %>%
    dplyr::left_join(descriptions, by = "feature_id") %>%
    dplyr::mutate(
      description_label = factor(description_label, levels = unique(descriptions$description_label)),
      is_target = ifelse(feature_id %in% target_present, "Target", "Other")
    )
  
  # Colors / legend order
  color_info <- get_similarity_ordered_colors(descriptions)
  
  # Subtitle must be scalar (prevents textbox_grob)
  target_desc <- descriptions %>%
    dplyr::filter(feature_id %in% target_present) %>%
    dplyr::pull(description) %>%
    unique()
  target_desc <- if (length(target_desc) == 0) "description not found" else paste(target_desc, collapse = "; ")
  
  subtitle_text <- paste0("<i><b>Target:</b> ", target_label, "</i> | <i>", target_desc, "</i>")
  subtitle_text <- paste0(subtitle_text, collapse = "")  
  
  plot_rates <- ggplot2::ggplot(combined_data, ggplot2::aes(x = as.numeric(year), y = Rate,
                                                            color = description_label, linetype = Sample,
                                                            group = interaction(feature_id, Sample), size = is_target)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_size_manual(values = c("Target" = 2, "Other" = 1), guide = "none") +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed")[seq_along(unique(combined_data$Sample))]) +
    ggplot2::scale_color_manual(values = color_info$color_map) +
    ggplot2::labs(title = paste("Rates for Related", Type, "Codes"),
                  subtitle = subtitle_text,
                  x = "", y = "Rate", color = Type, linetype = "Sample") +
    theme_global +
    ggplot2::scale_x_continuous(breaks = sort(unique(as.numeric(combined_data$year))))
  
  plot_counts <- ggplot2::ggplot(combined_data, ggplot2::aes(x = as.numeric(year), y = Patients,
                                                             color = description_label, linetype = Sample,
                                                             group = interaction(feature_id, Sample), size = is_target)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_size_manual(values = c("Target" = 2, "Other" = 1), guide = "none") +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed")[seq_along(unique(combined_data$Sample))]) +
    ggplot2::scale_color_manual(values = color_info$color_map) +
    ggplot2::labs(title = paste("Patient Counts for Related", Type, "Codes"),
                  subtitle = subtitle_text,
                  x = "", y = "Patients per Year", color = Type, linetype = "Sample") +
    theme_global +
    ggplot2::scale_x_continuous(breaks = sort(unique(as.numeric(combined_data$year))))
  
  # Total counts (bar) — same logic as your current, but safer for 1 vs 2 samples
  bar_counts <- purrr::imap_dfr(source_list, function(df, key) {
    df %>%
      dplyr::filter(feature_id %in% selected_features) %>%
      dplyr::distinct(patient_num, feature_id) %>%
      dplyr::group_by(feature_id) %>%
      dplyr::summarise(Total_Patients = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(Sample = sample_labels[[key]])
  }) %>%
    dplyr::left_join(descriptions, by = "feature_id") %>%
    dplyr::mutate(description_label = factor(description_label, levels = unique(descriptions$description_label))) %>%
    dplyr::filter(!is.na(description_label))
  
  sample_names <- unique(bar_counts$Sample)
  
  if (length(sample_names) == 1) {
    plot_bar <- ggplot2::ggplot(bar_counts, ggplot2::aes(x = description_label, y = Total_Patients, fill = description_label)) +
      ggplot2::geom_bar(stat = "identity", width = 0.6, color = "black") +
      ggplot2::scale_fill_manual(values = color_info$color_map) +
      ggplot2::labs(y = "Total Patients Across All Years", x = "") +
      theme_bar +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     legend.position = "none")
  } else if (length(sample_names) == 2) {
    plot_bar <- ggplot2::ggplot(bar_counts, ggplot2::aes(fill = description_label)) +
      ggplot2::geom_bar(data = bar_counts %>% dplyr::filter(Sample == sample_names[2]),
                        ggplot2::aes(x = as.numeric(description_label) + 0.1, y = Total_Patients),
                        stat = "identity", width = 0.4, color = "black", linetype = "dashed") +
      ggplot2::geom_bar(data = bar_counts %>% dplyr::filter(Sample == sample_names[1]),
                        ggplot2::aes(x = as.numeric(description_label) - 0.1, y = Total_Patients),
                        stat = "identity", width = 0.4, color = "black") +
      ggplot2::scale_x_continuous(breaks = unique(as.numeric(bar_counts$description_label)),
                                  labels = levels(bar_counts$description_label)) +
      ggplot2::scale_fill_manual(values = color_info$color_map) +
      ggplot2::labs(y = "Total Patients Across All Years", x = "") +
      theme_bar +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     legend.position = "none")
  } else {
    return(list(NULL, "plot_related_code_trends_v2 supports only 1 or 2 samples."))
  }
  
  combined_plot <- cowplot::plot_grid(
    plot_counts + ggplot2::theme(legend.position = "none"),
    plot_bar,
    ncol = 2,
    rel_widths = c(2, 1)
  )
  
  return(list(list(plot_rates, combined_plot), NULL))
}

# Final wrapper function (module 5)
plot_related_features <- function(data_inputs, target_code, target_cui, sample_labels,
                    O2 = FALSE,
                    manual_ONCE_path_code = NULL,
                    manual_ONCE_path_nlp = NULL,
                    types = c("Diagnosis", "Medication", "Lab", "Procedure", "CUI"),
                    type_dict = list("Diagnosis" = "PheCode", 
                                     "Medication" = "RXNORM", 
                                     "Lab" = c("LOINC","ShortName","Other Lab"), 
                                     "Procedure" = "CCS")) {
  
  # Load ONCE dictionary
  target_code_anchor <- target_code[1]
  ONCE <- clean_ONCE_data(target_code[1], O2, manual_ONCE_path_code, manual_ONCE_path_nlp)
  
  # Build codified and NLP lists
  codified_list <- list()
  nlp_list <- list()
  if (!is.null(data_inputs$codified1)) codified_list[["1"]] <- data_inputs$codified1
  if (!is.null(data_inputs$codified2)) codified_list[["2"]] <- data_inputs$codified2
  if (!is.null(data_inputs$nlp1)) nlp_list[["1"]] <- data_inputs$nlp1
  if (!is.null(data_inputs$nlp2)) nlp_list[["2"]] <- data_inputs$nlp2
  
  # Generate plots or messages for each selected type
  plot_results <- lapply(types, function(t) {
    result <- plot_related_code_trends_v2(
      Type = t,
      target_code = target_code,
      target_cui = target_cui,
      codified_list = codified_list,
      nlp_list = nlp_list,
      sample_labels = sample_labels,
      ONCE = ONCE,
      type_dict = type_dict
    )
    
    list(
      type = t,
      line_plot = if (!is.null(result[[1]])) result[[1]][[1]] else NULL,
      combined_plot = if (!is.null(result[[1]])) result[[1]][[2]] else NULL,
      message = if (is.null(result[[1]])) result[[2]] else NULL
    )
  })
  
  return(plot_results)
}

# ----------------------- MODULE 6 HELPERS -----------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(pROC)
})

# Confusion-matrix metrics
compute_conf_metrics <- function(truth, pred, positive = 1) {
  truth <- as.integer(truth)
  pred  <- as.integer(pred)
  
  TP <- sum(truth == positive & pred == positive, na.rm = TRUE)
  TN <- sum(truth != positive & pred != positive, na.rm = TRUE)
  FP <- sum(truth != positive & pred == positive, na.rm = TRUE)
  FN <- sum(truth == positive & pred != positive, na.rm = TRUE)
  
  sensitivity <- ifelse((TP + FN) == 0, NA_real_, TP / (TP + FN))
  specificity <- ifelse((TN + FP) == 0, NA_real_, TN / (TN + FP))
  ppv         <- ifelse((TP + FP) == 0, NA_real_, TP / (TP + FP))
  npv         <- ifelse((TN + FN) == 0, NA_real_, TN / (TN + FN))
  accuracy    <- ifelse((TP + TN + FP + FN) == 0, NA_real_, (TP + TN) / (TP + TN + FP + FN))
  
  tibble(
    TP = TP, TN = TN, FP = FP, FN = FN,
    sensitivity = sensitivity,
    specificity = specificity,
    ppv = ppv,
    npv = npv,
    accuracy = accuracy
  )
}

# ROC plot (theme + palette)
plot_roc_module6 <- function(roc_obj, title_text) {
  if (is.null(roc_obj)) return(NULL)
  
  coords_df <- pROC::coords(
    roc_obj, "all",
    ret = c("specificity", "sensitivity"),
    transpose = FALSE
  )
  
  roc_df <- tibble(
    FPR = 1 - coords_df[, "specificity"],
    TPR = coords_df[, "sensitivity"]
  )
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_line(aes(color = "ROC"), size = 1.2) +
    geom_abline(linetype = "dashed") +
    labs(
      title = title_text,
      subtitle = paste0("AUC = ", round(auc_val, 3)),
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = ""
    ) +
    theme_global + set_palette
}

# Metrics bar plot (Sens/Spec/PPV/NPV)
plot_metrics_module6 <- function(metrics_row, title_text) {
  metrics_long <- metrics_row %>%
    select(sensitivity, specificity, ppv, npv) %>%
    pivot_longer(everything(), names_to = "Metric", values_to = "Value") %>%
    mutate(
      Metric = factor(
        Metric,
        levels = c("sensitivity", "specificity", "ppv", "npv"),
        labels = c("Sensitivity", "Specificity", "PPV", "NPV")
      )
    )
  
  ggplot(metrics_long, aes(x = Metric, y = Value, fill = Metric)) +
    geom_col(width = 0.6) +
    geom_text(
      aes(label = ifelse(is.na(Value), "NA", scales::percent(Value, accuracy = 0.1))),
      vjust = -0.35,
      size = 4
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    labs(
      title = title_text,
      x = "",
      y = "Performance"
    ) +
    theme_global + set_fill_palette
}

# Confusion matrix tile plot
plot_confusion_module6 <- function(conf_row, title_text) {
  cm_df <- tibble(
    Truth = rep(c("Positive", "Negative"), each = 2),
    Predicted = rep(c("Positive", "Negative"), times = 2),
    Count = c(conf_row$TP, conf_row$FN, conf_row$FP, conf_row$TN),
    Cell = c("TP", "FN", "FP", "TN")
  )
  
  ggplot(cm_df, aes(x = Predicted, y = Truth, fill = Count)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Count), size = 6) +
    labs(
      title = title_text,
      x = "Predicted",
      y = "True"
    ) +
    theme_global + set_palette
}

# Main Module 6 function
#
# Evaluates labels against:
#   - ref_code: patient has >=1 target_code in codified{sample}
#   - ref_cui:  patient has >=1 target_cui  in nlp{sample}
#   - ref_both: patient has both
# Notes: Accepts binary labels (0/1). If you later use probabilities (0-1),
#   AUC/ROC remains valid; metrics use threshold (default 0.5).
# - Includes patients from labels OR codified OR nlp in the evaluation universe.
evaluate_labels_module6 <- function(data_inputs,
                                    target_code,
                                    target_cui,
                                    sample = "1",
                                    positive_label = 1,
                                    threshold = 0.5) {
  
  # ---- checks ----
  if (is.null(data_inputs$labels)) {
    stop("Module 6: data_inputs$labels not found. Provide labels via paths$labels or create them before Module 6.")
  }
  
  codified_name <- paste0("codified", sample)
  nlp_name      <- paste0("nlp", sample)
  
  if (is.null(data_inputs[[codified_name]])) stop("Module 6: Missing ", codified_name, " in data_inputs.")
  if (is.null(data_inputs[[nlp_name]]))      stop("Module 6: Missing ", nlp_name, " in data_inputs.")
  
  cod_df <- data_inputs[[codified_name]] %>% distinct(patient_num, feature_id)
  nlp_df <- data_inputs[[nlp_name]] %>% distinct(patient_num, feature_id)
  
  labels_df <- data_inputs$labels %>%
    transmute(
      patient_num = as.character(patient_num),
      label = as.numeric(label)
    ) %>%
    distinct(patient_num, .keep_all = TRUE)
  
  # ---- reference indicators ----
  ref_code <- cod_df %>%
    filter(feature_id %in% target_code) %>%
    distinct(patient_num) %>%
    mutate(ref_code = 1L)
  
  ref_cui <- nlp_df %>%
    filter(feature_id %in% target_cui) %>%
    distinct(patient_num) %>%
    mutate(ref_cui = 1L)
  
  all_patients <- tibble(
    patient_num = unique(c(
      as.character(cod_df$patient_num),
      as.character(nlp_df$patient_num),
      as.character(labels_df$patient_num)
    ))
  )
  
  eval_df <- all_patients %>%
    left_join(ref_code, by = "patient_num") %>%
    left_join(ref_cui,  by = "patient_num") %>%
    left_join(labels_df, by = "patient_num") %>%
    mutate(
      ref_code = ifelse(is.na(ref_code), 0L, ref_code),
      ref_cui  = ifelse(is.na(ref_cui),  0L, ref_cui),
      ref_both = as.integer(ref_code == 1L & ref_cui == 1L),
      label    = ifelse(is.na(label), 0, label),
      label_bin = as.integer(label >= threshold)
    )
  
  # ---- evaluator for a single reference ----
  eval_one <- function(truth_vec, label_score, label_bin, ref_name) {
    truth_vec <- as.integer(truth_vec)
    
    # ROC/AUC (works with binary or continuous predictor)
    roc_obj <- try(pROC::roc(truth_vec, label_score, quiet = TRUE), silent = TRUE)
    if (inherits(roc_obj, "try-error")) roc_obj <- NULL
    
    conf <- compute_conf_metrics(truth_vec, label_bin, positive = positive_label)
    
    metrics <- conf %>%
      mutate(
        AUC = ifelse(is.null(roc_obj), NA_real_, as.numeric(pROC::auc(roc_obj))),
        n_pos = sum(truth_vec == 1L, na.rm = TRUE),
        n_neg = sum(truth_vec == 0L, na.rm = TRUE),
        Reference = ref_name
      ) %>%
      select(Reference, n_pos, n_neg, AUC, sensitivity, specificity, ppv, npv, accuracy, TP, TN, FP, FN)
    
    plots <- list(
      roc_plot = plot_roc_module6(roc_obj, paste0("ROC vs ", ref_name)),
      metrics_plot = plot_metrics_module6(metrics, paste0("Metrics vs ", ref_name)),
      confusion_plot = plot_confusion_module6(conf, paste0("Confusion vs ", ref_name))
    )
    
    list(metrics = metrics, roc = roc_obj, confusion = conf, plots = plots)
  }
  
  res_code <- eval_one(eval_df$ref_code, eval_df$label, eval_df$label_bin, "Target Code")
  res_cui  <- eval_one(eval_df$ref_cui,  eval_df$label, eval_df$label_bin, "Target CUI")
  res_both <- eval_one(eval_df$ref_both, eval_df$label, eval_df$label_bin, "Code + CUI")
  
  list(
    metrics_table = bind_rows(res_code$metrics, res_cui$metrics, res_both$metrics),
    plots = list(
      code = res_code$plots,
      cui  = res_cui$plots,
      both = res_both$plots
    ),
    eval_dataframe = eval_df
  )
}

# Pretty metrics table for Module 6
print_module6_table <- function(results_module6,
                                caption = "Module 6: Label Performance vs Target Code and CUI") {
  
  if (is.null(results_module6$metrics_table)) {
    stop("metrics_table not found in results_module6.")
  }
  
  metrics_display <- results_module6$metrics_table %>%
    dplyr::mutate(
      AUC = round(AUC, 3),
      Sensitivity = scales::percent(sensitivity, accuracy = 0.1),
      Specificity = scales::percent(specificity, accuracy = 0.1),
      PPV = scales::percent(ppv, accuracy = 0.1),
      NPV = scales::percent(npv, accuracy = 0.1),
      Accuracy = scales::percent(accuracy, accuracy = 0.1)
    ) %>%
    dplyr::select(
      Reference,
      `N Positive` = n_pos,
      `N Negative` = n_neg,
      AUC,
      Sensitivity,
      Specificity,
      PPV,
      NPV,
      Accuracy, TP, TN, FP, FN
    )
  
  knitr::kable(
    metrics_display,
    caption = caption
  )
}
