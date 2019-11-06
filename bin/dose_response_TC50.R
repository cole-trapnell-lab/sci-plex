library(tidymodels)
library(purrr)
library(dplyr)
library(drc)
library(glmnet)


measure_pseudodose = function(cds_pdata, normalization_terms=c("cell_type", "culture_plate"), normalization_formula_str=NULL, pseudodose=0){
  sci_chem_counts = cds_pdata %>% 
    dplyr::filter(is.finite(Pseudotime)) %>%
    dplyr::group_by(cell_type, sample, replicate, culture_plate, well_oligo, treatment, dose, vehicle) %>% 
    dplyr::summarize(pseudodose = mean(Pseudotime)) %>% 
    dplyr::ungroup()
  
  sci_chem_counts = tibble::rownames_to_column(sci_chem_counts)
  
  #sci_chem_counts = sci_chem_counts %>% mutate(vehicle = dose %in% c(0)) # count_cells should do this automatically using the vehicle column in pData
  vehicle_counts = sci_chem_counts %>% dplyr::filter(vehicle) 
  
  sci_chem_counts = sci_chem_counts %>% dplyr::select(cell_type, sample, replicate, culture_plate, well_oligo, treatment, dose, pseudodose, vehicle)
  
  vehicle_counts = vehicle_counts %>% dplyr::select(cell_type, sample, replicate, culture_plate, well_oligo, pseudodose, vehicle)
  vehicle_count_rows = unique(sci_chem_counts$treatment) %>% furrr::future_map_dfr(.f = function(tmnt) {
    vehicle_count_rows = vehicle_counts
    vehicle_count_rows$treatment = tmnt
    vehicle_count_rows$dose = pseudodose
    vehicle_count_rows
  })
  # vehicle_count_rows = full_join(vehicle_counts, 
  #                                dplyr::select(sci_chem_counts, cell_type, sample, replicate, culture_plate, well_oligo,  treatment, vehicle) %>% distinct(),
  #                                by=c("cell_type", "sample", "replicate")) %>% distinct()
  # vehicle_count_rows$dose = pseudodose
  vehicle_count_rows  = vehicle_count_rows %>% dplyr::select(cell_type, sample, replicate, culture_plate, well_oligo, treatment, dose, vehicle, pseudodose)
  # FIXME: put this back!
  sci_chem_counts = rbind(sci_chem_counts, vehicle_count_rows)
  
  return(sci_chem_counts)
}


fit_drm_helper = function(df, max_pseudodose=NA, verbose=FALSE){
    tryCatch({
            #drc::drm(norm_cells ~ dose, data = df, fct=drc::LL.4(), type="Poisson" )
            drc::drm(pseudodose ~ dose, data = df, fct=drc::LL.4(fixed=c(NA,NA,max_pseudodose,NA)))
        }, 
        error = function(e) { 
            if (verbose)
                print (e); 
            NA 
        })
}

fit_dose_response_models = function(cell_count_table, max_pseudodose=NA, verbose=FALSE){
    sci_chem_models = cell_count_table %>% 
    tidyr::nest(-cell_type, -treatment) %>% 
    mutate(
        drc_model = purrr::map(.f = fit_drm_helper, .x = data, max_pseudodose=max_pseudodose, verbose=verbose)
    ) 
    sci_chem_models
}


extract_terms_helper = function(model){
     if (class(model) == "drc"){
      coefs = drc:::coef.drc(model)
      ed = tryCatch( { drc:::ED.drc(model, 50, interval="delta", display=FALSE) } , error = function(e) { ed = rep(NA, 4) })
      # Isolate the individual parameters and store in a small tibble
      tmp = tibble(
            steepness = coefs[1],
            lower_response_limit = coefs[2],
            upper_response_limit = coefs[3],
            #ec_param = coefs[3],
            ec50 = ed[1],
            ec50_lower = ed[3],
            ec50_upper = ed[4]
      )
      tmp
     }else{
         return (tibble(
            steepness = NA,
            lower_response_limit = NA,
            upper_response_limit = NA,
            ec50 = NA,
            ec50_lower = NA,
            ec50_upper = NA
      ))
     }
            
}

extract_terms = function(sci_chem_models){
    sci_chem_models = sci_chem_models %>%
        dplyr::mutate(drc.terms = furrr::future_map(.f = extract_terms_helper, .x = drc_model)) %>% 
        tidyr::unnest(drc.terms)
    return(sci_chem_models)
}



test_reduction_helper = function(model){
     if (class(model) == "drc"){
        tryCatch( {  
            cfs = drc:::coef.drc(model)
            #lower_response_limit = max(0, cfs[2])
            #upper_response_limit = max(0, cfs[3])
            lower_response_limit = predict(model, newdata=data.frame(dose=c(min(model$origData$dose))))
            upper_response_limit = predict(model, newdata=data.frame(dose=c(max(model$origData$dose))))
            fold_change = upper_response_limit / lower_response_limit
            #print (paste(upper_response_limit, lower_response_limit, fold_change))
            #print (noEffect(model))
            tbl = drc::noEffect(model)
            tbl = c(fold_change, tbl)
            
            #print (class(tbl))
            names(tbl) = c("fold_change", 
                              "Chi-square test",
                              "Df", 
                              "p-value") 
            tbl = as_tibble(t(tbl))
                              
            #print (tbl)
            return(tbl)
        }, error = function(e) { 
            print (e)
            tibble("fold_change" = NA,
                    "Chi-square test" = NA,
                    "Df" = NA, 
                    "p-value" = NA) 
        })
     }else{
        tibble("fold_change" = NA,
               "Chi-square test" = NA,
               "Df" = NA, 
               "p-value" = NA) 
     }     
}

test_dose_curve = function(sci_chem_models){
    sci_chem_sig_response = sci_chem_models %>%
        mutate(fit_test_res = furrr::future_map(.f = test_reduction_helper, .x = drc_model)) %>% 
        tidyr::unnest(fit_test_res) %>%
        dplyr::select(cell_type, treatment, "fold_change", "Chi-square test", "Df", "p-value")
    return(sci_chem_sig_response)
}


extract_dose_model_stats = function(sci_chem_counts, sci_chem_models){
    dose_range_df = sci_chem_counts %>% group_by(cell_type, treatment) %>% summarize(min_dose = min(dose), max_dose = max(dose))
    sci_chem_models = extract_terms(sci_chem_models)
    coefficient_df = sci_chem_models %>% dplyr::select(cell_type, treatment, ec50, ec50_lower, ec50_upper, upper_response_limit, lower_response_limit)
    coefficient_df$lower_response_limit = pmax(0, coefficient_df$lower_response_limit)
    coefficient_df$upper_response_limit = pmax(0, coefficient_df$upper_response_limit)
    coefficient_df = coefficient_df %>% mutate(fold_change = (upper_response_limit - lower_response_limit)/ upper_response_limit) 
    coefficient_df = left_join(coefficient_df, dose_range_df)
    coefficient_df$ec50[coefficient_df$ec50 > coefficient_df$max_dose] = NA
    coefficient_df$ec50_lower[coefficient_df$ec50_lower > coefficient_df$max_dose] = max(coefficient_df$max_dose) # FIXME: this isn't generic. Could be different max doses for different drugs
    coefficient_df$ec50_lower[coefficient_df$ec50_lower < coefficient_df$min_dose] = 0
    coefficient_df$ec50[coefficient_df$ec50 < coefficient_df$min_dose] = NA
    coefficient_df$ec50_upper[coefficient_df$ec50_upper > coefficient_df$max_dose] = max(coefficient_df$max_dose) # FIXME: this isn't generic. Could be different max doses for different drugs
    coefficient_df$ec50_upper[coefficient_df$ec50_upper < coefficient_df$min_dose] = 0
    return(coefficient_df)
}


gen_dose_range_helper = function(model){
     tryCatch({
                if (class(model) == "drc"){
                    eps_dose = min(model$origData$dose[model$origData$dose>0], na.rm=T)
                    exp(seq(log(min(model$origData$dose, na.rm=T) + eps_dose/2), log(max(model$origData$dose, na.rm=T)), length=100))
                }else{
                    return (rep(NA, 100))
                }
     }, error = function(e) { print(e); return (rep(NA, 100)) }
     )
}


extract_curve_helper = function(model){
     tryCatch({
                if (class(model) == "drc"){
                    eps_dose = min(model$origData$dose[model$origData$dose>0], na.rm=T)
                    newdata =  expand.grid(dose=exp(seq(log(min(model$origData$dose, na.rm=T) + eps_dose/2), log(max(model$origData$dose, na.rm=T)), length=100)))
                    #newdata = data.frame(dose=seq(min(model$origData$dose), max(model$origData$dose),length.out=100))
                    res = predict(model, newdata=newdata)
                }else{
                    return (rep(NA, 100)) 
                }
     }, 
     warning = function(e) { return (rep(NA, 100)) },
     error = function(e) { return (rep(NA, 100)) }
     )
}

extract_curves = function(sci_chem_models){
    sci_chem_curves = sci_chem_models %>%
        mutate(dose = furrr::future_map(.f = gen_dose_range_helper, .x = drc_model)) %>%
        mutate(fit = furrr::future_map(.f = extract_curve_helper, .x = drc_model)) %>%
        dplyr::select(cell_type, treatment,  dose, fit) %>% 
        tidyr::unnest(dose, fit)
    return(sci_chem_curves)
}

plot_dose_response_curves = function(sci_chem_data, sci_chem_models, color_by="replicate", ncol=NULL){
    # sci_chem_data = sci_chem_models %>%
    #     #mutate(fit = furrr::future_map(.f = extract_curve, .x = drc_model)) %>%
    #     dplyr::select(cell_type, treatment, data) %>% 
    #     unnest(data)
    coefficient_df = extract_dose_model_stats(sci_chem_data, sci_chem_models)

    sci_chem_curves = extract_curves(sci_chem_models)
    min_plot_dose = 0.05
    ggplot(aes(dose + min_plot_dose, pseudodose), data=sci_chem_data) + 
      geom_point(aes_string(color=color_by), position="jitter", stoke = 0, size = .75) +
        geom_line(aes(dose, fit), data=sci_chem_curves) + 
        geom_rect(aes(x=ec50, y=0, xmin=ec50_lower + min_plot_dose, xmax=ec50_upper, ymin=0, ymax=Inf), alpha=I(0.1), data=coefficient_df) + 
        geom_vline(aes(xintercept=ec50), linetype="dotted", data=coefficient_df) + 
        scale_x_log10() + 
        scale_color_brewer(palette="Set1") + monocle:::monocle_theme_opts() + 
        facet_wrap(treatment~cell_type, scales="free_y", ncol=ncol)
}
