## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)


## ----setup--------------------------------------------------------------------
# devtools::load_all()
# library(dplyr)
# library(parallel)

## -----------------------------------------------------------------------------
# #select data for a single chemical and species
# my_data <- subset(cvt,
#                  analyte_dtxsid %in% c("DTXSID8031865",
#                                             "DTXSID0020442")
#                  )

## -----------------------------------------------------------------------------
# my_pk <- pk(data = my_data)

## -----------------------------------------------------------------------------
#   my_pk <- my_pk +
#   settings_preprocess(suppress.messages = FALSE,
#                       keep_data_original = TRUE)

## -----------------------------------------------------------------------------
# #preprocess data
# my_pk <- do_preprocess(my_pk)

## -----------------------------------------------------------------------------
# #get summary data info
# my_pk <- do_data_info(my_pk)

## -----------------------------------------------------------------------------
# #pre-fitting (get model parameter bounds & starting values)
# my_pk <- do_prefit(my_pk)

## -----------------------------------------------------------------------------
# #model fitting
# my_pk <- do_fit(my_pk)

## -----------------------------------------------------------------------------
# system.time(my_pk <- do_fit(my_pk)) #without parallel processing

## -----------------------------------------------------------------------------
# system.time(my_pk <- do_fit(my_pk,
#                 n_cores = parallel::detectCores()-1))

## ----eval = FALSE-------------------------------------------------------------
# my_pk <- pk(data = my_data)
# 
# my_pk <- do_fit(my_pk)

## -----------------------------------------------------------------------------
# coef(my_pk)

## -----------------------------------------------------------------------------
# coef(my_pk) %>% slice(1) %>%  mutate(as.data.frame(as.list(coefs_vector[[1]]))) %>%
#   glimpse()

## -----------------------------------------------------------------------------
# coef(my_pk) %>% slice(5) %>%  mutate(as.data.frame(as.list(coefs_vector[[1]]))) %>%
#   glimpse()

## -----------------------------------------------------------------------------
# my_resids <- residuals(my_pk)
# my_resids %>% glimpse()

## -----------------------------------------------------------------------------
# my_preds <- predict(my_pk)
# my_preds %>% glimpse()

## -----------------------------------------------------------------------------
# predict(my_pk, newdata = data.frame(
#   Chemical = "DTXSID8031865",
#   Species = "rat",
#   Time = seq(0, 5, by = 0.5),
#   Time.Units = "hours",
#   Conc.Units = "mg/kg",
#   Dose = 1,
#   Route = "iv",
#   Media = "plasma",
#   exclude = FALSE)) %>% glimpse()

## -----------------------------------------------------------------------------
# p <- plot(x = my_pk)
# print(p$final_plot)

## -----------------------------------------------------------------------------
# p2 <- plot(my_pk,
#      use_scale_conc = list(dose_norm = TRUE,
#                            log10_trans = FALSE)
#      )
# print(p2$final_plot)

## -----------------------------------------------------------------------------
# logLik(my_pk)

## -----------------------------------------------------------------------------
# AIC(my_pk)

## -----------------------------------------------------------------------------
# BIC(my_pk)

## -----------------------------------------------------------------------------
# rmse(my_pk)

## -----------------------------------------------------------------------------
# rsq(my_pk)

## -----------------------------------------------------------------------------
# AFE(my_pk)

## -----------------------------------------------------------------------------
# AAFE(my_pk)

## -----------------------------------------------------------------------------
# data_proc <- get_data(my_pk)
# data_proc %>% glimpse()

## -----------------------------------------------------------------------------
# my_data_info <- get_data_info(my_pk)
# names(my_data_info)
# 
# my_data_info$data_flags %>% glimpse()

## -----------------------------------------------------------------------------
# my_tkstats <- get_tkstats(my_pk, suppress.messages = TRUE)
# my_tkstats %>% glimpse()

## -----------------------------------------------------------------------------
# wm <- get_winning_model(my_pk)
# wm

## -----------------------------------------------------------------------------
# rsq_df <- rsq(my_pk)
# 
# rsq_win <- wm %>% dplyr::left_join(rsq_df)
# 
# rsq_win

## -----------------------------------------------------------------------------
# plot(my_pk, best_fit = TRUE, use_scale_conc = list(dose_norm = TRUE,
#                                                    log10_trans = FALSE)) %>%
#   pull(final_plot)

## -----------------------------------------------------------------------------
# my_pk <- pk(data = my_data)

## -----------------------------------------------------------------------------
# names(my_pk)

## -----------------------------------------------------------------------------
# head(my_pk$data_original)
# print(my_pk)

## ----eval = FALSE-------------------------------------------------------------
# ggplot(data = my_data,
#        aes(x = time_hr,
#            y = conc,
#            color = dose_level_normalized_corrected,
#            shape = as.factor(document_id)
#        )
# )

## ----eval = FALSE-------------------------------------------------------------
# pk(data = my_data,
#    mapping = ggplot2::aes(
#                   mapping = ggplot2::aes(
#                    Chemical = analyte_dtxsid,
#                    Chemical_Name = analyte_name_original,
#                    DTXSID = analyte_dtxsid,
#                    CASRN = analyte_casrn,
#                    Species = species,
#                    Reference = fk_extraction_document_id,
#                    Media = conc_medium_normalized,
#                    Route = administration_route_normalized,
#                    Dose = invivPK_dose_level,
#                    Dose.Units = "mg/kg",
#                    Subject_ID = fk_subject_id,
#                    Series_ID = fk_series_id,
#                    Study_ID = fk_study_id,
#                    ConcTime_ID = conc_time_id,
#                    N_Subjects = n_subjects_normalized,
#                    Weight = weight_kg,
#                    Weight.Units = "kg",
#                    Time = time_hr,
#                    Time.Units = "hours",
#                    Value = invivPK_conc,
#                    Value.Units = "mg/L",
#                    Value_SD = invivPK_conc_sd,
#                    LOQ = invivPK_loq
#                  ))

## -----------------------------------------------------------------------------
# names(cvt)

## ----eval = FALSE-------------------------------------------------------------
# Reference = as.character(
#   ifelse(
#     is.na(
#       documents_reference.id
#     ),
#     documents_extraction.id,
#     documents_reference.id
#   )
# )
# 

## ----eval=FALSE---------------------------------------------------------------
# Value.Units = "mg/L"

## -----------------------------------------------------------------------------
# get_status(my_pk)

## -----------------------------------------------------------------------------
# my_pk <- pk(my_data)

## ----eval = FALSE-------------------------------------------------------------
# my_pk <-   pk(my_data) +
#   facet_data(facets = dplyr::vars(Chemical, Species))

## ----eval = FALSE-------------------------------------------------------------
# my_pk <- pk(my_data)

## ----eval = FALSE-------------------------------------------------------------
# my_pk <- pk(my_data) + settings_preprocess()

## -----------------------------------------------------------------------------
# formals(settings_preprocess)

## ----eval = FALSE-------------------------------------------------------------
# my_pk <- pk(my_data)

## ----eval = FALSE-------------------------------------------------------------
# my_pk <- pk(my_data) +
#   settings_preprocess() +
#   settings_data_info()

## -----------------------------------------------------------------------------
# formals(settings_data_info)

## ----eval = FALSE-------------------------------------------------------------
# my_pk <- pk(my_data) +
#   settings_preprocess() +
#   settings_data_info() +
#   settings_optimx()

## -----------------------------------------------------------------------------
# formals(settings_optimx)

## ----eval = FALSE-------------------------------------------------------------
# my_pk <- pk(my_data) +
#   settings_preprocess() +
#   settings_data_info() +
#   settings_optimx() +
#   scale_conc()

## -----------------------------------------------------------------------------
# formals(scale_conc)

## ----eval = FALSE-------------------------------------------------------------
# my_pk <- pk(my_data) +
#   settings_preprocess() +
#   settings_data_info() +
#   settings_optimx() +
#   scale_conc() +
#   scale_time()

## -----------------------------------------------------------------------------
# formals(scale_time)

## ----eval = FALSE-------------------------------------------------------------
# my_pk <- pk(my_data) +
#   settings_preprocess() +
#   settings_data_info() +
#   settings_optimx() +
#   scale_conc() +
#   scale_time() +
#   stat_model()

## -----------------------------------------------------------------------------
# formals(stat_model)

## ----eval = FALSE-------------------------------------------------------------
# my_pk <-pk(my_data) +
#   settings_preprocess() +
#   settings_data_info() +
#   settings_optimx() +
#   scale_conc() +
#   scale_time() +
#   stat_model() +
#   stat_error_model()

## -----------------------------------------------------------------------------
# formals(stat_error_model)

## ----eval=FALSE---------------------------------------------------------------
# hier_pk <- my_pk +
#   facet_data(ggplot2::vars(Chemical, Species)) +
#   stat_error_model(error_group = vars(Chemical, Species, Reference))

## ----eval=FALSE---------------------------------------------------------------
# pooled_pk <- my_pk +
#   facet_data(dplyr::vars(Chemical, Species)) +
#   stat_error_model(error_group = vars(Chemical, Species))

## ----eval=FALSE---------------------------------------------------------------
# separate_pk <- my_pk +
#   facet_data(dplyr::vars(Chemical, Species, Reference)) +
#   stat_error_model(error_group = vars(Chemical, Species, Reference))

## ----eval = FALSE-------------------------------------------------------------
# my_pk <- pk(data = my_data)
# my_pk <- do_fit(my_pk)

## -----------------------------------------------------------------------------
# my_pk <- pk(my_data) +
#     #instructions to use an error model that puts all observations in the same group
#   stat_error_model(error_group = ggplot2::vars(Chemical, Species)) +
#   #instructions for concentration scaling/transformation
#   scale_conc(dose_norm = TRUE,
#              log10_trans = TRUE) +
#   #instructions for time rescaling
#   scale_time(new_units = "auto") +
#   #instructions to use only one method for optimx::optimx()
#   settings_optimx(method = "L-BFGS-B") +
#   #instructions to impute missing LOQs slightly differently
#   settings_preprocess(calc_loq_factor = 0.5)

## -----------------------------------------------------------------------------
# get_data_original(my_pk)

## -----------------------------------------------------------------------------
# get_mapping(my_pk)

## -----------------------------------------------------------------------------
# get_status(my_pk)

## -----------------------------------------------------------------------------
# get_data_group(my_pk)

## -----------------------------------------------------------------------------
# get_settings_preprocess(my_pk)

## -----------------------------------------------------------------------------
# get_settings_data_info(my_pk)

## -----------------------------------------------------------------------------
# get_settings_optimx(my_pk)

## -----------------------------------------------------------------------------
# get_scale_conc(my_pk)

## -----------------------------------------------------------------------------
# get_scale_time(my_pk)

## -----------------------------------------------------------------------------
# get_stat_model(my_pk)

## -----------------------------------------------------------------------------
# get_stat_error_model(my_pk)

## -----------------------------------------------------------------------------
# my_pk <- my_pk +
#   scale_conc(dose_norm = TRUE, log10_trans = FALSE)

## -----------------------------------------------------------------------------
# get_scale_conc(my_pk)

## -----------------------------------------------------------------------------
# my_pk <- my_pk + stat_model(model = "model_1comp")

## -----------------------------------------------------------------------------
# get_stat_model(my_pk)

## -----------------------------------------------------------------------------
# my_pk <- pk(data = my_data)
# names(my_pk)

## -----------------------------------------------------------------------------
# my_pk <- do_preprocess(my_pk)

## -----------------------------------------------------------------------------
# names(my_pk)

## -----------------------------------------------------------------------------
# get_data(my_pk) %>% glimpse()

## -----------------------------------------------------------------------------
# my_pk <- do_data_info(my_pk)

## -----------------------------------------------------------------------------
# names(my_pk)

## -----------------------------------------------------------------------------
# my_data_info <- get_data_info(my_pk) #extracts the `data_info` element as a named list
# names(my_data_info)

## -----------------------------------------------------------------------------
# my_nca <- get_nca(my_pk) #extracts the `data_info$nca` element specifically

## -----------------------------------------------------------------------------
# my_pk <- do_prefit(my_pk)

## -----------------------------------------------------------------------------
# names(my_pk)

## -----------------------------------------------------------------------------
# get_prefit(my_pk) %>% names()

## -----------------------------------------------------------------------------
# my_pk$prefit$stat_error_model$sigma_DF %>% glimpse()

## -----------------------------------------------------------------------------
# my_pk$prefit$par_DF %>% glimpse()

## -----------------------------------------------------------------------------
# my_pk$prefit$fit_check %>% glimpse()

## -----------------------------------------------------------------------------
# my_pk <- do_fit(my_pk)

## -----------------------------------------------------------------------------
# names(my_pk)

## -----------------------------------------------------------------------------
# get_fit(my_pk) %>% glimpse()

## -----------------------------------------------------------------------------
# my_pk <- pk(data = my_data) + #initialize a `pk` object
#   stat_model(model = c("model_flat",
#                        "model_1comp",
#                        "model_2comp")) + #add PK models to fit
#   settings_optimx(method = "L-BFGS-B") #use only this optimx::optimx() algorithm
# 
# get_status(my_pk) #status is 1

## -----------------------------------------------------------------------------
# my_pk <- do_fit(my_pk)
# get_status(my_pk)

## -----------------------------------------------------------------------------
# AIC(my_pk, suppress.messages = TRUE)

## -----------------------------------------------------------------------------
# my_pk <- my_pk + scale_conc(dose_norm = TRUE)

## -----------------------------------------------------------------------------
# get_status(my_pk)

## ----error = TRUE-------------------------------------------------------------
try({
# coef(my_pk) #throws an error
})

## -----------------------------------------------------------------------------
# my_pk <- do_fit(my_pk)
# get_status(my_pk)

## -----------------------------------------------------------------------------
# AIC(my_pk, suppress.messages = TRUE)

## -----------------------------------------------------------------------------
# #fit a pk object
# my_pk <- pk(data = my_data) +
#   settings_preprocess(suppress.messages = TRUE)
# 
# my_pk <- do_fit(my_pk)
# 
# get_status(my_pk) #status is now 5
# 
# #copy it to a new variable
# my_pk_new <- my_pk
# 
# #and modify scale_conc() on the new copy
# my_pk_new <- my_pk + scale_conc(dose_norm = TRUE)
# 
# get_status(my_pk_new) #status has been reset to 1 for the new copy
# get_status(my_pk) #but status of the original is still 5

## -----------------------------------------------------------------------------
# #suppress messages
# my_pk <- my_pk + settings_preprocess(suppress.messages = TRUE)
# 
# #dose normalization, no log transformation
# pk1 <- my_pk + scale_conc(dose_norm = TRUE, log10_trans = FALSE)
# 
# #log transformation, no dose normalization
# pk2 <- my_pk + scale_conc(dose_norm = FALSE, log10_trans = TRUE)
# 
# #both normalization and log transformation
# pk3 <- my_pk + scale_conc(dose_norm = TRUE, log10_trans = TRUE)
# 
# #do fits
# pk1 <- do_fit(pk1, n_cores = parallel::detectCores()-1)
# pk2 <- do_fit(pk2, n_cores = parallel::detectCores()-1)
# pk3 <- do_fit(pk3, n_cores = parallel::detectCores()-1)

## -----------------------------------------------------------------------------
# plot(pk1) %>% pull(final_plot)

## -----------------------------------------------------------------------------
# plot(pk2) %>% pull(final_plot)

## -----------------------------------------------------------------------------
# plot(pk3) %>% pull(final_plot)

## -----------------------------------------------------------------------------
# plot(pk1, use_scale_conc = list(dose_norm = TRUE,
#                                 log10_trans = FALSE)) %>% pull(final_plot)

## -----------------------------------------------------------------------------
# plot(pk2, use_scale_conc = list(dose_norm = TRUE,
#                                 log10_trans = FALSE)) %>% pull(final_plot)

## -----------------------------------------------------------------------------
# plot(pk3, use_scale_conc = list(dose_norm = TRUE,
#                                 log10_trans = FALSE)) %>% pull(final_plot)

## -----------------------------------------------------------------------------
# pk_hier <- my_pk +
#   stat_error_model(error_group = vars(Chemical, Species, Reference))
# 
# pk_pool <- my_pk +
#   stat_error_model(error_group = vars(Chemical, Species))
# 
# 
# pk_hier <- do_fit(pk_hier, n_cores = parallel::detectCores()-1)
# pk_pool <- do_fit(pk_pool, n_cores = parallel::detectCores()-1)

## -----------------------------------------------------------------------------
# plot(pk_hier, use_scale_conc = list(dose_norm = TRUE,
#                                 log10_trans = FALSE)) %>% pull(final_plot)

## -----------------------------------------------------------------------------
# plot(pk_pool, use_scale_conc = list(dose_norm = TRUE,
#                                 log10_trans = FALSE)) %>% pull(final_plot)

## -----------------------------------------------------------------------------
# `model_1comp`

## -----------------------------------------------------------------------------
# `model_2comp`

