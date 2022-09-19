### Pharmacy Mapper ###
# Takes in a classification term and outputs pharmacy rows that match
#' Retrieve rows from GEMINI pharmacy data matching specified drug term(s).
#'
#' @param db_con PostgreSQL connection class
#' The connection object for the desired database and user. Obtained with odbc::dbConnect()
#' Supports clean_v1_4_0 and up
#' @param class_input string or list of strings
#' Classification keyword or specific ids like ATC are accepted.
#' @param drug_input string or list of strings
#' Generic or brand name is accepted. Spelling is normalized.
#' @param min_match_score (optional, default: 0) integer
#' The minimum match score needed for a match to be considered, from 0 to 100.
#' Recommended to use a lower score to avoid false negatives, then manually review for false positives.
#' @param detailed_search (optional, default: T) logical
#' If true, search for every related concept to each selected drug, instead of just the selected drugs.
#' This will greatly expand the search but will sometimes match related concepts that are not desired.
#' @param sites (optional, default: 'all') string or list of strings
#' Reduces search to only certain sites, speeding up the query if not all are needed.
#' Site names should match the hospital_id code format like 'SMH' or 'UHNTG'. Search all using 'all'.
#'  @param genc_ids (optional, default: 'all') integer or list of integers which defines an encounter identifier. 
#' This parameter reduces search to only certain encounters based on their integer id, speeding up the query if a specified list of encounters is known. 
#' @param validation_mode (optional, default: F) logical
#' Outputs a condensed file of matches made with unique combinations of key fields
#' Enables easier validation of the matches but does not return full pharmacy data
#' @param return_drug_list (optional, default: F) logical
#' Outputs the search drug list instead of searching
#' @return
#' The GEMINI pharmacy dataframe for matched rows, or NA if cancelled.
#' An additional 'rxnorm_match' column is added specifying what rxnorm matched the row to.
#'
#' @examples
#' \dontrun{
#' diabetes_orders_condensed <- gemini_rxnorm_query(db_con = con,
#'                                        drug_input = c("metformin", "insulin"),
#'                                        class_input = "diabetes",
#'                                        validation_mode = T)
#' }
#'
#' @import RCurl odbc httr jsonlite DT getPass data.table dplyr
#'
#' @export

gemini_rxnorm_query <- function(db_con, class_input = NA, drug_input = NA, ...){
  
  # Assert there is atleast one input else stop 
  if(is.na(class_input[1]) & is.na(drug_input[1])){
    stop("At least one of class_input or drug_input must not be NA.")
  }
  # Extract the optional parameters
  opt <- list(...)
  if(!is.integer(opt$min_match_score)) opt$min_match_score <- 0
  if(!is.logical(opt$detailed_search)) opt$detailed_search <- T
  if(is.null(opt$sites)){
    opt$sites <- "all"
  }else if(opt$sites[1] != "all"){
    opt$sites <- paste0("(", paste(paste0("'", toupper(opt$sites), "'"), collapse=","), ")")
  }
  if(is.null(opt$genc_ids)){
    opt$genc_ids <- "all"
  }else if(opt$genc_ids[1] != "all"){
    opt$genc_ids <- paste(paste0("(", unique(opt$genc_ids), ")"), collapse=",")
  }
  if(!is.logical(opt$validation_mode)) opt$validation_mode <- F
  if(!is.logical(opt$return_drug_list)) opt$return_drug_list <- F

  if(!is.na(class_input[1])){
    ###### CLASSIFICATION SEARCH ######
    # Find ATC classes with search_input in name or id
    call <- paste0("https://rxnav.nlm.nih.gov/REST/",
                   "rxclass/allClasses.json?classTypes=","ATC1-4")
    api_call_text <- httr::content(GET(call), "text", encoding = "UTF-8")
    api_call_df_raw <- as.data.frame(jsonlite::fromJSON(api_call_text, flatten = TRUE))
    all_classes <- api_call_df_raw %>%
      rename(class_id = rxclassMinConceptList.rxclassMinConcept.classId,
             class_name = rxclassMinConceptList.rxclassMinConcept.className,
             class_type = rxclassMinConceptList.rxclassMinConcept.classType)


    ###
    class_list <- data.table()
    for(class in class_input){
      skip_to_next <- FALSE
      class_list_i <- all_classes[grep(class, paste(all_classes$class_name, all_classes$class_id),
                                       ignore.case = T), ]

      if(nrow(class_list_i)==0){
        #Get spelling suggestions for class names
        tryCatch({
          call <- paste0("https://rxnav.nlm.nih.gov/REST/",
                         "rxclass/spellingsuggestions.json?type=", "DRUG",
                         "&term=", RCurl::curlEscape(class))

          api_call_text <- httr::content(GET(call), "text", encoding = "UTF-8")
          api_call_df_raw <- as.data.frame(jsonlite::fromJSON(api_call_text, flatten = TRUE))

          spelling_suggestions <- api_call_df_raw %>%
            select(suggestion)

          cat("\nWARNING:", class, "was skipped because no classes have this keyword.",
              "\n\tPotential spelling suggestions:", paste(spelling_suggestions$suggestion, collapse=", "), "\n")
        }, error = function(e) {
          cat("\nWARNING:", class, "was skipped because no classes have this keyword.\n")
        }, finally={
          skip_to_next <<- TRUE
        })
      }
      if(skip_to_next) next
      class_list <- rbind(class_list, class_list_i, fill=T)
    }
    if(nrow(class_list)==0){
      cat("\nCould not find any class matches for your search input.\n")
      return(NA)
    }

    class_list <- distinct(class_list) %>%
      arrange(class_id)
    rownames(class_list) <- NULL

    # Prompt user to select an ATC class
    print(datatable(class_list))
    cat("\nThe table in the viewer displays the drug classes containing your search term.
        Press enter to confirm all, or enter the indexes to remove separated by commas.
        (or enter c to cancel)")
    class_input <- readline()
    if(class_input == "c"){
      cat("Cancelled.")
      return(NA)
    }
    if(class_input != ""){
      remove_indexes <- as.numeric(unlist(strsplit(class_input, ",")))
      cat("Removing", length(remove_indexes), "classes from search list...")
      class_list <- class_list %>% filter(!row_number() %in% remove_indexes)
    }
    cat("\nFinding drugs belonging to selected classes...\n")

    # Find drugs belonging to ATC classes
    class_id_list <- class_list$class_id

    drug_list <- data.table()
    term_types <- "BN+BPCK+DF+DFG+GPCK+IN+MIN+PIN+SBD+SBDC+SBDF+SBDG+SCD+SCDC+SCDF+SCDG"
    for(class_id in class_id_list){
      skip_to_next <- FALSE
      drug_list_i <- data.table()
      tryCatch({
        call <- paste0("https://rxnav.nlm.nih.gov/REST/",
                       "rxclass/classMembers.json?classId=",RCurl::curlEscape(class_id),
                       "&relaSource=", "ATC",
                       "&trans=", 0,
                       "&ttys=",term_types)
        api_call_text <- httr::content(GET(call), "text", encoding = "UTF-8")
        api_call_df_raw <- as.data.frame(jsonlite::fromJSON(api_call_text, flatten = TRUE))

        drug_list_i <- api_call_df_raw %>%
          rename(rxcui = drugMemberGroup.drugMember.minConcept.rxcui,
                 drug_name = drugMemberGroup.drugMember.minConcept.name) %>%
          select(rxcui, drug_name) %>%
          mutate(atc_class = class_id)
      }, error = function(e) {
        skip_to_next <<- TRUE
      })
      drug_list <- rbind(drug_list, drug_list_i, fill=T)
    }
    if(nrow(drug_list) == 0){
      stop("Could not find any active drugs belonging to the selected class(es) in RxNorm")
    }
    drug_list_class <- distinct(drug_list) %>%
      group_by_at(setdiff(names(drug_list), "atc_class")) %>%
      summarise(atc_class = paste(atc_class, collapse = ", "), .groups = 'drop')

  }

  if(!is.na(drug_input[1])){
    ###### DRUG NAME SEARCH ######
    # Find drugs with the names in list
    drug_list <- data.table()
    for(drug in drug_input){
      skip_to_next <- FALSE
      drug_list_i <- data.table()
      tryCatch({
        call <- paste0("https://rxnav.nlm.nih.gov/REST/",
                       "rxcui.json?name=",RCurl::curlEscape(drug),
                       "&allsrc=",0,
                       "&search=",2)
        api_call_text <- httr::content(GET(call), "text", encoding = "UTF-8")
        drug_list_i <- as.data.frame(fromJSON(api_call_text, flatten = TRUE)) %>%
          dplyr::rename("rxcui" = "rxnormId") %>%
          mutate(atc_class = "MANUAL",
                 drug_name = drug) %>%
          select(rxcui, drug_name, atc_class) %>%
          head(1)
      }, error = function(e) {
        #Get spelling suggestions for drug names
        tryCatch({
          call <- paste0("https://rxnav.nlm.nih.gov/REST/",
                         "rxclass/spellingsuggestions.json?type=", "DRUG",
                         "&term=", RCurl::curlEscape(drug))

          api_call_text <- httr::content(GET(call), "text", encoding = "UTF-8")
          api_call_df_raw <- as.data.frame(jsonlite::fromJSON(api_call_text, flatten = TRUE))

          spelling_suggestions <- api_call_df_raw %>%
            select(suggestion)

          cat("\nWARNING:", drug, "was skipped because it could not be found.",
              "\n\tPotential spelling suggestions:", paste(spelling_suggestions$suggestion, collapse=", "), "\n")
        }, error = function(e) {
          cat("\nWARNING:", drug, "was skipped because it could not be found.\n")
        }, finally={
          skip_to_next <<- TRUE
        })
      })
      if(skip_to_next) next
      drug_list <- rbind(drug_list, drug_list_i, fill=T)

    }
    if(nrow(drug_list)==0){
      cat("\nCould not find any matches for your search input.\n")
      return(NA)
    }
    drug_list_drug <- distinct(drug_list)
  }
  if(!is.na(drug_input[1]) & !is.na(class_input[1])){
    rxcui_overlap <- intersect(drug_list_drug$rxcui, drug_list_class$rxcui)
    drug_list_drug <- drug_list_drug %>% filter(!rxcui %in% rxcui_overlap)
    drug_list <- drug_list_class %>%
      mutate(atc_class = ifelse(rxcui %in% rxcui_overlap, paste0(atc_class, ", MANUAL"), atc_class)) %>%
      rbind(drug_list_drug) %>%
      arrange(drug_name)
  } else if(!is.na(drug_input[1])){
    drug_list <- drug_list_drug
  } else if(!is.na(class_input[1])){
    drug_list <- drug_list_class
  } else stop("No drugs found.")


  ###### RXCUI SEARCH ######
  # Prompt user to confirm the drugs
  print(datatable(drug_list))
  cat("\nThe table in the viewer displays the drugs that will be searched.
      Press enter to confirm all, or enter the indexes to remove separated by commas.
      (or enter c to cancel)")
  drug_input <- readline()
  if(drug_input == "c"){
    cat("Cancelled.")
    return(NA)
  }
  if (drug_input != ""){
    remove_indexes <- as.numeric(unlist(strsplit(drug_input, ",")))
    cat("Removing", length(remove_indexes), "rxcui from search list...")
    drug_list <- drug_list %>% filter(!row_number() %in% remove_indexes)
  }
  if(opt$return_drug_list) return(drug_list)
  drug_list <- drug_list %>% select(rxcui, drug_name) %>% as.data.frame()

  if(opt$detailed_search){
    # Optional: search for the rxcui of all concepts related to selected drugs
    cat("\nSearching for selected drugs and their related concepts...\n")
    related_drug_list <- data.table()
    term_types <- "BN+BPCK+DF+DFG+GPCK+IN+MIN+PIN+SBD+SBDC+SBDF+SBDG+SCD+SCDC+SCDF+SCDG"
    for(drug_row in 1:nrow(drug_list)){
      skip_to_next <- FALSE
      related_drug_list_i <- data.table()
      tryCatch({
        call <- paste0("https://rxnav.nlm.nih.gov/REST/",
                       "rxcui/", RCurl::curlEscape(drug_list[drug_row, "rxcui"]),
                       "/related.json?tty=", term_types)
        api_call_text <- httr::content(GET(call), "text", encoding = "UTF-8")

        api_call_df_raw <- as.data.frame(rbindlist(fromJSON(api_call_text, flatten = TRUE)$relatedGroup$conceptGroup$conceptProperties))
        related_drug_list_i <- api_call_df_raw %>%
          select(rxcui) %>%
          mutate(drug_name = drug_list[drug_row, "drug_name"])
      }, error = function(e) {
        skip_to_next <<- TRUE
      })
      if(skip_to_next) next
      related_drug_list <- rbind(related_drug_list, related_drug_list_i, fill=T)
      Sys.sleep(0.05) # So that not too many searches get sent
    }

    drug_list <- related_drug_list %>%
      group_by_at(setdiff(names(related_drug_list), c("drug_name"))) %>%
      dplyr::summarise(drug_name = paste(drug_name, collapse = " | "), .groups = 'drop') %>%
      distinct()

  }else cat("\nSearching for selected drugs....\n")

  search_rxcui_list <- paste(paste0("(", drug_list$rxcui, ")"), collapse = ",")

  # Pull the pharmacy rows that match the rxcui list
  # NOTE: This will have to be changed to your database!

  query_str <- paste0(
    paste0("WITH rxcui_selected (rxcui) AS (VALUES", search_rxcui_list, ")"),
    ifelse(opt$genc_ids[1] != "all", paste0(", genc_ids_selected (genc_id) AS (VALUES", opt$genc_ids,")"), ""),

   " SELECT p.*, ndc_rxcui, din_rxcui, gen_rxcui, bran_rxcui, hos_rxcui, iv_rxcui,",
   " GREATEST(ndc_score, din_score, gen_score, bran_score, hos_score, iv_score) as rxnorm_top_score",
   " FROM pharmacy p",
   ifelse(opt$genc_ids[1] != "all", paste0(" INNER JOIN genc_ids_selected g ON g.genc_id = p.genc_id"), ""),

   " LEFT JOIN (SELECT rc.rxcui as din_rxcui, raw_input, score as din_score FROM rxnorm_cache rc",
   " INNER JOIN rxcui_selected rs ON rs.rxcui = rc.rxcui",
   " WHERE search_type = 'med_id_din'", " AND active = TRUE",
   ifelse(opt$min_match_score != 0, paste0(" AND score >= ", opt$min_match_score, ") rx_din"), ") rx_din"),
   " ON p.med_id_din = rx_din.raw_input",

   " LEFT JOIN (SELECT rc.rxcui as ndc_rxcui, raw_input, score as ndc_score FROM rxnorm_cache rc",
   " INNER JOIN rxcui_selected rs ON rs.rxcui = rc.rxcui",
   " WHERE search_type = 'med_id_ndc'", " AND active = TRUE",
   ifelse(opt$min_match_score != 0, paste0(" AND score >= ", opt$min_match_score, ") rx_ndc"), ") rx_ndc"),
   " ON p.med_id_ndc = rx_ndc.raw_input",

   " LEFT JOIN (SELECT rc.rxcui as gen_rxcui, raw_input, score as gen_score FROM rxnorm_cache rc",
   " INNER JOIN rxcui_selected rs ON rs.rxcui = rc.rxcui",
   " WHERE search_type = 'med_id_generic_name_raw'", " AND active = TRUE",
   ifelse(opt$min_match_score != 0, paste0(" AND score >= ", opt$min_match_score, ") rx_gen"), ") rx_gen"),
   " ON p.med_id_generic_name_raw = rx_gen.raw_input",

   " LEFT JOIN (SELECT rc.rxcui as bran_rxcui, raw_input, score as bran_score FROM rxnorm_cache rc",
   " INNER JOIN rxcui_selected rs ON rs.rxcui = rc.rxcui",
   " WHERE search_type = 'med_id_brand_name_raw'", " AND active = TRUE",
   ifelse(opt$min_match_score != 0, paste0(" AND score >= ", opt$min_match_score, ") rx_bran"), ") rx_bran"),
   " ON p.med_id_brand_name_raw = rx_bran.raw_input",

   " LEFT JOIN (SELECT rc.rxcui as hos_rxcui, raw_input, score as hos_score FROM rxnorm_cache rc",
   " INNER JOIN rxcui_selected rs ON rs.rxcui = rc.rxcui",
   " WHERE search_type = 'med_id_hospital_code_raw'", " AND active = TRUE",
   ifelse(opt$min_match_score != 0, paste0(" AND score >= ", opt$min_match_score, ") rx_hos"), ") rx_hos"),
   " ON p.med_id_hospital_code_raw = rx_hos.raw_input AND p.hospital_id IN ('MSH', 'THPC')", #Only search hospital code for MSH and THPC

   " LEFT JOIN (SELECT rc.rxcui as iv_rxcui, raw_input, score as iv_score FROM rxnorm_cache rc",
   " INNER JOIN rxcui_selected rs ON rs.rxcui = rc.rxcui",
   " WHERE search_type = 'iv_component_type'", " AND active = TRUE",
   ifelse(opt$min_match_score != 0, paste0(" AND score >= ", opt$min_match_score, ") rx_iv"), ") rx_iv"),
   " ON p.med_id_hospital_code_raw = rx_iv.raw_input AND p.hospital_id IN ('THPM', 'THPC')", #Only search iv_component for THPC and THPM

   " WHERE COALESCE(din_rxcui, ndc_rxcui, gen_rxcui, bran_rxcui, hos_rxcui, iv_rxcui) IS NOT NULL",
   ifelse(opt$sites[1] != "all", paste0(" AND hospital_id IN ", opt$sites, ";"), ";")
    )
    cat("Searching database. This may take a few minutes...\n")
    # Pull the data from the databases
    tryCatch({
      pharm_matches <- odbc::dbGetQuery(db_con, query_str)
    }, error = function(e) {
      print(e)
      stop("Error occured when querying database. Please re-create your database connection and try again.")
    })

    cat("Finalizing. This may take a few minutes...\n")
    # Attach the name of matched drugs, combined from all matches
    final_matches <- pharm_matches %>%
      merge(drug_list, by.x="din_rxcui", by.y="rxcui", all.x=T) %>%
      dplyr::rename(din_name = drug_name) %>%
      merge(drug_list, by.x="ndc_rxcui", by.y="rxcui", all.x=T) %>%
      dplyr::rename(ndc_name = drug_name) %>%
      merge(drug_list, by.x="gen_rxcui", by.y="rxcui", all.x=T) %>%
      dplyr::rename(gen_name = drug_name) %>%
      merge(drug_list, by.x="bran_rxcui", by.y="rxcui", all.x=T) %>%
      dplyr::rename(bran_name = drug_name) %>%
      merge(drug_list, by.x="hos_rxcui", by.y="rxcui", all.x=T) %>%
      dplyr::rename(hos_name = drug_name) %>%
      merge(drug_list, by.x="iv_rxcui", by.y="rxcui", all.x=T) %>%
      dplyr::rename(iv_name = drug_name) %>%
      group_by_at(setdiff(names(pharm_matches), c("din_name", "ndc_name", "gen_name", "bran_name", "hos_name", "iv_name",
                                                  "din_rxcui", "ndc_rxcui", "gen_rxcui", "bran_rxcui",
                                                  "hos_rxcui", "iv_rxcui", "rxnorm_top_score"))) %>%
      summarise(din_name = paste(din_name, collapse = " | "),
                ndc_name = paste(ndc_name, collapse = " | "),
                gen_name = paste(gen_name, collapse = " | "),
                bran_name = paste(bran_name, collapse = " | "),
                hos_name = paste(hos_name, collapse = " | "),
                iv_name = paste(iv_name, collapse = " | "),
                rxnorm_top_score = max(rxnorm_top_score)
                , .groups = 'drop') %>%
      mutate(rxnorm_match = paste(din_name, ndc_name, gen_name, bran_name, hos_name, iv_name, sep = " | ")) %>%
      relocate(c(rxnorm_top_score, rxnorm_match), .before=med_id_generic_name_raw) %>%
      select(-din_name, -ndc_name, -gen_name, -bran_name, -hos_name, -iv_name) %>%
      as.data.table()


    final_matches$rxnorm_match <- unname(sapply(final_matches$rxnorm_match, function(x) {
      trimws(sub(paste(sort(trimws(unique(strsplit(x, split=" | ", fixed=T)[[1]]))), collapse=' || '), pattern = "NA \\|\\||\\|\\| NA", replacement = ""))} ))

    final_matches <- distinct(final_matches)
    final_matches_n <- nrow(final_matches)

    if(opt$validation_mode){
      final_matches <- final_matches %>%
        group_by(hospital_id, rxnorm_top_score, rxnorm_match, med_id_generic_name_raw, med_id_brand_name_raw,
                 med_id_hospital_code_raw, iv_component_type, med_id_din, med_id_ndc) %>%
        summarise(occurences = n(), .groups="drop") %>%
        ungroup() %>%
        relocate(occurences, .before=rxnorm_top_score) %>%
        arrange(desc(occurences))
    }

    cat("Found",final_matches_n,"matches.")
    return(final_matches)
}


# con <- dbConnect(odbc::odbc(), Driver = "PostgreSQL ODBC Driver(UNICODE)",
#                  Server = "geminidb", Port = "PORT_NUMBER",
#                  Database = "DB_NAME", UID = "USER",
#                  PWD = getPass::getPass("Database password: "))
#
#
#
# ####### TESTING AREA #######
# diabetes_orders_condensed <- gemini_rxnorm_query(db_con = con,
#                                        drug_input = c("furosemide"),
#                                        validation_mode = T)
