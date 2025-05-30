test_that("creating a pk object works for CvT data for DTXSID3061635",
          {
            expect_no_error(my_pk <- pk(
              data = subset(cvt, analyte_dtxsid %in% "DTXSID3061635")
            )
            )
          }
)

test_that("pk object has the expected default mapping",
          {
            my_pk <- pk(
              data = subset(cvt,
                            analyte_dtxsid %in% "DTXSID3061635"
              )
            )
            expect_equal(my_pk$mapping,
                         ggplot2::aes(
                           Chemical = analyte_dtxsid,
                           Chemical_Name = analyte_name_original,
                           DTXSID = analyte_dtxsid,
                           CASRN = analyte_casrn,
                           Species = species,
                           Reference = fk_extraction_document_id,
                           Media = conc_medium_normalized,
                           Route = administration_route_normalized,
                           Dose = invivPK_dose_level,
                           Dose.Units = "mg/kg",
                           Subject_ID = fk_subject_id,
                           Series_ID = fk_series_id,
                           Study_ID = fk_study_id,
                           ConcTime_ID = conc_time_id,
                           N_Subjects = n_subjects_normalized,
                           Weight = weight_kg,
                           Weight.Units = "kg",
                           Time = time_hr,
                           Time.Units = "hours",
                           Value = invivPK_conc,
                           Value.Units = "mg/L",
                           Value_SD = invivPK_conc_sd,
                           LOQ = invivPK_loq
                         ),
                         ignore_attr = TRUE)
          })

test_that("creating a pk object with a mapping missing all required variables throws the expected warning",
          {
            expect_warning(
              pk(
                data = subset(cvt,
                              analyte_dtxsid %in% "DTXSID3061635"
                ),
                mapping = NULL
              )
            )
          }
)

test_that("A newly-created pk object has the expected status of 1",
          {
            my_pk <- pk(
              data = subset(cvt,
                            analyte_dtxsid %in% "DTXSID3061635"
              )
            )
            expect_equal(my_pk$status, 1L)

          }
          )


