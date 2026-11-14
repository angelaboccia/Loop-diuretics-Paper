##Adverse Events to Loop Diuretics combined with guideline-directed medical 
# therapies for Heart Failure: a pharmacovigilance analysis
# Authors: Angela Boccia, Jacopo Palombarini

# Setup ------------------------------------------------------------------------
# Download ---------------------------------------------------------------------
# Download DiAna package
# devtools::install_github("fusarolimichele/DiAna_package")

# Packages ---------------------------------------------------------------------
library(DiAna)
library(readxl)
library(here)
library(dplyr)
library(writexl)
library(tidyverse)
library(patchwork)
library(data.table)

## Data ------------------------------------------------------------------------
# setup_DiAna("25Q3") # only if you don't have 25Q3 FAERS data yet.

# import data
import("DEMO","25Q3")
import("DRUG", "25Q3")
import("REAC","25Q3")
import("INDI","25Q3")
import("OUTC","25Q3")
import("THER","25Q3")
import_MedDRA()
import_ATC()

# De-duplication by event date, sex, age, weight, adverse events, and suspected drugs
Demo <- Demo[RB_duplicates_only_susp == FALSE]

#path <- paste0(here::here(), "/projects/Diuretici paper/")

# Add continents to Demo and keep only US reports

# countries <- setDT(read.csv2(paste0(here::here(), "/external_sources/Countries.csv")))
# countries <- countries[, .(Continent_Name, Country_Name)]
# colnames(countries) <- c("continent", "reporter_country")
# countries <- countries[!is.na(reporter_country)]
# countries <- distinct(countries)
# Demo <- left_join(Demo, countries)
# cont <- c("North America", "Africa", "Oceania", "Europe", "Asia", "South America", "Antarctica")

Reac <- Reac[primaryid %in% Demo$primaryid]
Drug <- Drug[primaryid %in% Demo$primaryid]
Indi <- Indi[primaryid %in% Demo$primaryid]
Outc <- Outc[primaryid %in% Demo$primaryid]
Ther <- Ther[primaryid %in% Demo$primaryid]

## Cases selection -------------------------------------------------------------
# Heart failure indication -----------------------------------------------------
smqs <- extractSMQ(Narrow = FALSE)
SMQs_selected <- list("cardiac failure" = as.list(smqs[["Cardiac failure (SMQ)"]]))
indi_pt_selected <- unique(unlist(flatten(SMQs_selected)))
pids_indi <- unique(Indi[indi_pt %in% indi_pt_selected]$primaryid)

# Drugs selection --------------------------------------------------------------
## Diuretics -------------------------------------------------------------------
diuretics <- c("furosemide",
               "piretanide",
               "torasemide",
               "etacrynic acid",
               "bumetanide",
               "tienilic acid")

# Raw
pids_diuretics      <- unique(Drug[substance %in% diuretics]$primaryid)
# Indi
pids_diuretics_indi <- unique(Drug[substance %in% diuretics & primaryid %in% pids_indi]$primaryid)

## Guideline-directed medical therapy (GDMT)------------------------------------
# Defination of ARNI class
pids_sacubitril <- unique(Drug[substance == "sacubitril"]$primaryid)
pids_valsartan  <- unique(Drug[substance == "valsartan"]$primaryid)
pids_arni       <- intersect(pids_sacubitril, pids_valsartan)

Drug_GDMT_source <- copy(Drug)
Drug_GDMT_source[, substance := as.character(substance)]
Drug_GDMT_source[primaryid %in% pids_arni & substance %in% c("sacubitril", "valsartan"),
                   substance := "sacubitril/valsartan"]

ace_inhibitors <- c("captopril", "enalapril", "lisinopril", "ramipril", "trandolapril", "perindopril")
arb            <- c("candesartan", "valsartan", "losartan", "irbesartan", "telmisartan")
beta_blockers  <- c("bisoprolol", "carvedilol", "metoprolol", "nebivolol")
mra            <- c("eplerenone", "spironolactone", "canrenone")
arni           <- c("sacubitril/valsartan")
sglt           <- c("dapagliflozin", "empagliflozin", "canagliflozin")

all_GDMT <- c(ace_inhibitors, arb, beta_blockers, mra, arni, sglt)

pids_GDMT               <- unique(Drug[substance %in% all_GDMT, ]$primaryid)
pids_GDMT_indi          <- unique(intersect(pids_GDMT, pids_indi))

### Pids classes GDMT and pids diuretics and only one GDMT ---------------------
Drug_GDMT      <- copy(Drug_GDMT_source)[primaryid %in% pids_GDMT & substance %in% all_GDMT, ]
Drug_GDMT_indi <- copy(Drug_GDMT_source)[primaryid %in% pids_GDMT_indi & substance %in% all_GDMT, ]

GDMT_classes <- list("ACE inhibitors"    = ace_inhibitors, 
                       "ARBs"              = arb, 
                       "Beta-blockers"     = beta_blockers,
                       "MRAs"              = mra,
                       "ARNI"              = arni,
                       "SGLT-2 inhibitors" = sglt)

for (class_name in names(GDMT_classes)) {
  Drug_GDMT[substance %in% GDMT_classes[[class_name]], class := class_name]
}
Drug_GDMT <- Drug_GDMT%>%
  na.omit()
for (class_name in names(GDMT_classes)) {
  Drug_GDMT_indi[substance %in% GDMT_classes[[class_name]], class := class_name]
}
Drug_GDMT_indi <- Drug_GDMT_indi %>%
  na.omit()

# Raw
GDMT_class_count <- Drug_GDMT[, .(n_classes = uniqueN(class)), by = primaryid]
# One class
pids_one_class <- unique(GDMT_class_count[n_classes == 1, primaryid])
# Indi
GDMT_class_count_indi <- Drug_GDMT_indi[, .(n_classes = uniqueN(class)), by = primaryid]
# One class
pids_one_class_indi <- unique(GDMT_class_count_indi[n_classes == 1, primaryid])

# Raw
# Diuretics + ACEi
pids_ace                    <- unique(Drug_GDMT[substance %in% ace_inhibitors,]$primaryid)
pids_diuretic_ace           <- intersect(pids_diuretics, pids_ace)
pids_diuretic_ace_one_class <- intersect(pids_diuretic_ace, pids_one_class)
# Diuretics + ARB
pids_arb                    <- unique(Drug_GDMT[substance %in% arb]$primaryid)
pids_diuretic_arb           <- intersect(pids_diuretics, pids_arb)
pids_diuretic_arb_one_class <- intersect(pids_diuretic_arb, pids_one_class)
# Diuretics + beta-bloccanti
pids_beta                    <- unique(Drug_GDMT[substance %in% beta_blockers]$primaryid)
pids_diuretic_beta           <- intersect(pids_diuretics, pids_beta)
pids_diuretic_beta_one_class <- intersect(pids_diuretic_beta, pids_one_class)
# Diuretics + MRA
pids_mra                    <- unique(Drug_GDMT[substance %in% mra]$primaryid)
pids_diuretic_mra           <- intersect(pids_diuretics, pids_mra)
pids_diuretic_mra_one_class <- intersect(pids_diuretic_mra, pids_one_class)
# Diuretics + ARNI
pids_arni                    <- unique(Drug_GDMT[substance %in% arni]$primaryid)
pids_diuretic_arni           <- intersect(pids_diuretics, pids_arni)
pids_diuretic_arni_one_class <- intersect(pids_diuretic_arni, pids_one_class)
# Diuretics + SGLT-2
pids_sglt                    <- unique(Drug_GDMT[substance %in% sglt]$primaryid)
pids_diuretic_sglt           <- intersect(pids_diuretics, pids_sglt)
pids_diuretic_sglt_one_class <- intersect(pids_diuretic_sglt, pids_one_class)

# Indi
# Diuretics + ACEi
pids_diuretic_ace_indi           <- intersect(pids_diuretics_indi, pids_ace)
pids_diuretic_ace_one_class_indi <- intersect(pids_diuretic_ace_indi, pids_one_class_indi)
# Diuretics + ARB
pids_diuretic_arb_indi           <- intersect(pids_diuretics_indi, pids_arb)
pids_diuretic_arb_one_class_indi <- intersect(pids_diuretic_arb_indi, pids_one_class_indi)
# Diuretics + beta-bloccanti
pids_diuretic_beta_indi           <- intersect(pids_diuretics_indi, pids_beta)
pids_diuretic_beta_one_class_indi <- intersect(pids_diuretic_beta_indi, pids_one_class_indi)
# Diuretics + MRA
pids_diuretic_mra_indi           <- intersect(pids_diuretics_indi, pids_mra)
pids_diuretic_mra_one_class_indi <- intersect(pids_diuretic_mra_indi, pids_one_class_indi)
# Diuretics + ARNI
pids_diuretic_arni_indi           <- intersect(pids_diuretics_indi, pids_arni)
pids_diuretic_arni_one_class_indi <- intersect(pids_diuretic_arni_indi, pids_one_class_indi)
# Diuretics + SGLT-2
pids_diuretic_sglt_indi           <- intersect(pids_diuretics_indi, pids_sglt)
pids_diuretic_sglt_one_class_indi <- intersect(pids_diuretic_sglt_indi, pids_one_class_indi)

# Diuretic + GDMTs 
all_drugs <- c(diuretics, all_GDMT)

pids_diuretic_GDMT      <- unique(intersect(pids_diuretics, pids_GDMT))
pids_diuretic_GDMT_indi <- unique(intersect(pids_diuretic_GDMT, pids_indi))

## Reac selection --------------------------------------------------------------
smqs_broad  <- extractSMQ(Narrow = FALSE)
smqs_narrow <- extractSMQ(Narrow = TRUE)

## Cardiac arrhythmias ---------------------------------------------------------
# Broad extraction with exclusions ---------------------------------------------
filter_SMQ_arr <- as.list(smqs_broad$`Cardiac arrhythmias (SMQ)`)

# exclude
palpitations <- as.list(grep("palpitation", filter_SMQ_arr, value = TRUE))
heart_rate   <- as.list(grep("heart rate", filter_SMQ_arr, value = TRUE)[1:3])
ecg          <- as.list(grep("electro", filter_SMQ_arr, value = TRUE))
loc          <- as.list(grep("loss of", filter_SMQ_arr, value = TRUE))
sync         <- as.list(grep("syncope", filter_SMQ_arr, value = TRUE))
foetal       <- as.list(grep("foetal", filter_SMQ_arr, value = TRUE))
neonatal     <- as.list(grep("neonatal", filter_SMQ_arr, value = TRUE))
reflex       <- as.list(grep("reflex", filter_SMQ_arr, value = TRUE))
cerebro      <- as.list(grep("cerebro", filter_SMQ_arr, value = TRUE))
rebound      <- as.list(grep("rebound", filter_SMQ_arr, value = TRUE))
sinus        <- as.list(grep("respiratory sinus arrh", filter_SMQ_arr, value = TRUE))

exclude_arr <- list(palpitations, heart_rate, ecg, loc, sync, foetal, neonatal, reflex, cerebro, rebound, sinus)

# include
bradycardia              <- as.list(grep("bradycardia", filter_SMQ_arr, value = TRUE)[1])
cardiac_arrest           <- as.list(grep("cardiac arrest", filter_SMQ_arr, value = TRUE)[1])
cardiorespiratory_arrest <- as.list(grep("cardio-respiratory arrest", filter_SMQ_arr, value = TRUE)[1])
sudden_death             <- as.list(grep("sudden death", filter_SMQ_arr, value = TRUE))
cardiac_death            <- as.list(grep("cardiac death", filter_SMQ_arr, value = TRUE))
tachycardia              <- as.list(grep("tachycardia", filter_SMQ_arr, value = TRUE)[c(2,3)])

include_arr <- list(bradycardia, cardiac_arrest, cardiorespiratory_arrest, sudden_death, cardiac_death, tachycardia)

SMQs_selected_A_broad <- list("arrhythmias broad" = as.list(include_arr))

# Narrow SMQ extraction 
SMQs_selected_A_narrow <- list("arrhythmias narrow"                      = as.list(smqs_narrow[["Cardiac arrhythmias (SMQ)"]]),
                               "arrhythmias terms"                       = as.list(smqs_narrow[["Cardiac arrhythmia terms (incl bradyarrhythmias and tachyarrhythmias) (SMQ)"]]),
                               "arrhythmias symptoms"                    = as.list(smqs_narrow[["Arrhythmia related investigations, signs and symptoms (SMQ)"]]),
                               "bradyarrhythmias"                        = as.list(smqs_narrow[["Bradyarrhythmias (incl conduction defects and disorders of sinus node function) (SMQ)"]]),
                               "tachyarrhythmias"                        = as.list(smqs_narrow[["Tachyarrhythmias (incl supraventricular and ventricular tachyarrhythmias) (SMQ)"]]),
                               "supraventricular tachyarrhythmias"       = as.list(smqs_narrow[["Supraventricular tachyarrhythmias (SMQ)"]]),
                               "ventricular tachyarrhythmias"            = as.list(smqs_narrow[["Ventricular tachyarrhythmias (SMQ)"]]),
                               "QT prolongation"                         = as.list(smqs_narrow[["Torsade de pointes/QT prolongation (SMQ)"]]))

SMQs_selected_A <- c(SMQs_selected_A_broad, SMQs_selected_A_narrow)
event_arr <- setdiff(unique(unlist(flatten(SMQs_selected_A))), exclude_arr)

# list_pt <- data.frame(event_arr)
# write_xlsx(list_pt, paste0(path, "list_pt_ARR.xlsx"))

### Acute renal failure -----
filter_SMQ_ren <- as.list(smqs_broad$`Acute renal failure (SMQ)`)

#### Broad extraction with exclusions -----
# exclude
prote <- as.list(grep("protein", filter_SMQ_ren, value = TRUE))
albume <- as.list(grep("albumin", filter_SMQ_ren, value = TRUE))
crystal <- as.list(grep("crystal", filter_SMQ_ren, value = TRUE))
transplant <- as.list(grep("transplant", filter_SMQ_ren, value = TRUE))
nutrition <- as.list(grep("nutrition", filter_SMQ_ren, value = TRUE))
kidney <- as.list(grep("kidney injury molec", filter_SMQ_ren, value = TRUE))
NGAL <- as.list(grep("neutrophil", filter_SMQ_ren, value = TRUE))

exclude_ren <- list(prote, albume, NGAL, crystal, transplant, nutrition, kidney)

# include
blood_creatinine_increased            <- as.list(grep("blood creatinine increased", filter_SMQ_ren, value = TRUE))
blood_creatinine_abnormal             <- as.list(grep("blood creatinine abnormal", filter_SMQ_ren, value = TRUE))
glomerular_filtration_rate_decreased  <- as.list(grep("glomerular filtration rate decreased", filter_SMQ_ren, value = TRUE))
glomerular_filtration_rate_abnormal   <- as.list(grep("glomerular filtration rate abnormal", filter_SMQ_ren, value = TRUE))
blood_urea_increased                  <- as.list(grep("blood urea increased", filter_SMQ_ren, value = TRUE))
blood_urea_abnormal                   <- as.list(grep("blood urea abnormal", filter_SMQ_ren, value = TRUE))
creatinine_ratio_increased            <- as.list(grep("creatinine ratio increased", filter_SMQ_ren, value = TRUE))
urine_output_decreased                <- as.list(grep("urine output decreased", filter_SMQ_ren, value = TRUE))
renal_tubular_necrosis                <- as.list(grep("renal tubular necrosis", filter_SMQ_ren, value = TRUE))
renal_tubular_dysfunction             <- as.list(grep("renal tubular dysfunction", filter_SMQ_ren, value = TRUE))
renal_tubular_injury                  <- as.list(grep("renal tubular injury", filter_SMQ_ren, value = TRUE))
nephritis                             <- as.list(grep("nephritis", filter_SMQ_ren, value = TRUE))

include_ren <- list(blood_creatinine_increased, blood_creatinine_abnormal, glomerular_filtration_rate_decreased, 
                    glomerular_filtration_rate_abnormal, blood_urea_increased, blood_urea_abnormal, 
                    creatinine_ratio_increased, urine_output_decreased, renal_tubular_necrosis, 
                    renal_tubular_dysfunction, renal_tubular_injury, nephritis)

#### Narrow SMQ extraction -----
SMQs_selected_B <- list("renal failure broad" = as.list(include_ren),
                        "renal failure narrow" = as.list(smqs_narrow[["Acute renal failure (SMQ)"]]))

event_ren       <- setdiff(unique(unlist(flatten(SMQs_selected_B))), exclude_ren)

#list_pt <- data.frame(event_ren)
#write_xlsx(list_pt, paste0(path, "list_pt_REN.xlsx"))

### Syncope -----
SMQs_selected_C <- list("syncope" = as.list(c("syncope", 
                                            "presyncope", 
                                            "loss of consciousness", 
                                            "orthostatic hypotension", 
                                            "blood pressure orthostatic decreased", 
                                            "postural orthostatic tachycardia syndrome",
                                            "circulatory collapse")))

event_syn       <- unique(unlist(flatten(SMQs_selected_C)))
#list_pt <- data.frame(event_syn)
#write_xlsx(list_pt, paste0(path, "list_pt_syn.xlsx"))

all_pt <- unique(c(event_arr, event_ren, event_syn))

# Raw
temp         <- left_join(Reac[primaryid %in% pids_diuretic_GDMT], 
                          Drug[primaryid %in% pids_diuretic_GDMT])
pt_reac      <- unique(temp[pt %in% all_pt][, .N,by= pt][order(-N)])
# Indi
temp_indi    <- left_join(Reac[primaryid %in% pids_diuretic_GDMT_indi],
                          Drug[primaryid %in% pids_diuretic_GDMT_indi])
pt_reac_indi <- unique(temp_indi[pt %in% all_pt][,.N,by = pt][order(-N)])

## Descriptive analysis --------------------------------------------------------
# Descriptive over cases per HF ------------------------------------------------
#path <- paste0(here::here(), "/projects/Diuretici paper/results/")
descriptive(pids_cases = pids_diuretic_GDMT_indi, drug = diuretics, 
            file_name = paste0(path, "descriptive_diur_GDMT_indi.xlsx"))

# GDMT classes combinations frequencies ----------------------------------------
# Raw
GDMT_class_list      <- Drug_GDMT[!is.na(class), .(n_classes = uniqueN(class),
                                                   classes = paste(sort(unique(class)), collapse = " + ")),
                                  by = primaryid]

combo_counts           <- GDMT_class_list[, .(n_primaryid = .N), by = .(n_classes, classes)][
  order(n_classes, -n_primaryid)][, percentage := round((n_primaryid / sum(n_primaryid)) * 100, 2)]

# Indi
GDMT_class_list_indi <- Drug_GDMT_indi[!is.na(class), .(n_classes = uniqueN(class),
                                                        classes = paste(sort(unique(class)), collapse = " + ")),
                                       by = primaryid]

combo_counts_indi      <- GDMT_class_list_indi[, .(n_primaryid = .N), by = .(n_classes, classes)][
  order(n_classes, -n_primaryid)][, percentage := round((n_primaryid / sum(n_primaryid)) * 100, 2)]


# GDMT classes combinations with diuretics ----
# Raw
GDMT_combinations <- GDMT_class_list[
  primaryid %in% pids_diuretic_GDMT][
    , .(n_primaryid = .N), by = .(n_classes, classes)][
      , classes := paste0("Diuretics + ", classes)][
        order(n_classes, -n_primaryid)][
          , percentage := round((n_primaryid / sum(n_primaryid)) * 100, 2)]

# Indi
GDMT_combinations_indi <- GDMT_class_list_indi[
  primaryid %in% pids_diuretic_GDMT_indi][
    , .(n_primaryid = .N), by = .(n_classes, classes)][
      , classes := paste0("Diuretics + ", classes)][
        order(n_classes, -n_primaryid)][
          , percentage := round((n_primaryid / sum(n_primaryid)) * 100, 2)]

# write_xlsx(GDMT_combinations, paste0(path, "diuretic_concomitant_GDMT_classes_combinations.xlsx"))
# write_xlsx(GDMT_combinations_indi, paste0(path, "diuretic_concomitant_GDMT_classes_combinations_indi.xlsx"))

# Full GDMT classes concomitant with diuretics ----
# Raw
GDMT_classes_full <- Drug_GDMT[
  primaryid %in% pids_diuretic_GDMT &
    substance %in% all_GDMT &
    !(substance %in% diuretics) &
    !is.na(class),
  .(N = uniqueN(primaryid)),
  by = class][
    order(-N)][
      , percentage := round((N / uniqueN(pids_diuretic_GDMT)) * 100, 2)]

# Indi
GDMT_classes_full_indi <- Drug_GDMT_indi[
  primaryid %in% pids_diuretic_GDMT_indi &
    substance %in% all_GDMT &
    !(substance %in% diuretics) &
    !is.na(class),
  .(N = uniqueN(primaryid)),
  by = class][
    order(-N)][
      , percentage := round((N / uniqueN(pids_diuretic_GDMT_indi)) * 100, 2)]

# write_xlsx(GDMT_classes_full, paste0(path, "diuretic_concomitant_GDMT_classes_full.xlsx"))
# write_xlsx(GDMT_classes_full_indi, paste0(path, "diuretic_concomitant_GDMT_classes_full_indi.xlsx"))

# GDMT substances concomitant with diuretics ----
# Raw
concomitant <- Drug_GDMT[
  primaryid %in% pids_diuretic_GDMT &
    substance %in% all_GDMT &
    !(substance %in% diuretics) &
    !is.na(substance)]

total_concomitant_primaryid <- uniqueN(concomitant$primaryid)

concomitant_perc <- concomitant[
  , .(N = uniqueN(primaryid)), by = substance][
    , percentage := round((N / total_concomitant_primaryid) * 100, 2)][
      order(-N)]

# Indi
concomitant_indi <- Drug_GDMT_indi[
  primaryid %in% pids_diuretic_GDMT_indi &
    substance %in% all_GDMT &
    !(substance %in% diuretics) &
    !is.na(substance)]

total_concomitant_primaryid_indi <- uniqueN(concomitant_indi$primaryid)

concomitant_perc_indi <- concomitant_indi[
  , .(N = uniqueN(primaryid)), by = substance][
    , percentage := round((N / total_concomitant_primaryid_indi) * 100, 2)][
      order(-N)]

# write_xlsx(concomitant_perc, paste0(path, "diuretic_concomitant_GDMT_drugs.xlsx"))
# write_xlsx(concomitant_perc_indi, paste0(path, "diuretic_concomitant_GDMT_drugs_indi.xlsx"))

# Number of concomitant GDMT classes with diuretics ----
# Raw
class_counts <- concomitant[!is.na(class), .(n_classes = uniqueN(class)), by = primaryid]
total_patients <- uniqueN(class_counts$primaryid)

summary_counts <- data.table(
  category = c("1 class", "2 classes", "3-4 classes", ">4 classes"),
  N = c(sum(class_counts$n_classes == 1),
        sum(class_counts$n_classes == 2),
        sum(class_counts$n_classes >= 3 & class_counts$n_classes <= 4),
        sum(class_counts$n_classes > 4)))[
          , percentage := round((N / total_patients) * 100, 2)]

# Indi
class_counts_indi <- concomitant_indi[!is.na(class), .(n_classes = uniqueN(class)), by = primaryid]
total_patients_indi <- uniqueN(class_counts_indi$primaryid)

summary_counts_indi <- data.table(
  category = c("1 class", "2 classes", "3-4 classes", ">4 classes"),
  N = c(sum(class_counts_indi$n_classes == 1),
        sum(class_counts_indi$n_classes == 2),
        sum(class_counts_indi$n_classes >= 3 & class_counts_indi$n_classes <= 4),
        sum(class_counts_indi$n_classes > 4)))[
          , percentage := round((N / total_patients_indi) * 100, 2)]


write_xlsx(combo_counts, paste0(path,"combo_classes_frequencies.xlsx"))
write_xlsx(combo_counts_indi, paste0(path, "combo_classes_frequencies_indi.xlsx"))
write_xlsx(GDMT_combinations, paste0(path, "GDMT_combinations_with_diuretics.xlsx"))
write_xlsx(GDMT_combinations_indi, paste0(path, "GDMT_combinations_with_diuretics_indi.xlsx"))
write_xlsx(concomitant_perc, paste0(path, "diuretic_concomitant_GDMT_substances.xlsx"))
write_xlsx(concomitant_perc_indi, paste0(path, "diuretic_concomitant_GDMT_substances_indi.xlsx"))
write_xlsx(summary_counts, paste0(path, "diuretic_number_GDMT_classes.xlsx"))
write_xlsx(summary_counts_indi, paste0(path, "diuretic_number_GDMT_classes_indi.xlsx"))

# Reports with at least one selected event ----
report_selected_events_indi <- unique(Reac[primaryid %in% pids_diuretic_GDMT_indi & pt %in% all_pt]$primaryid)

## Percentage of reports with the three classes of events ----
Reac_temp_indi <- copy(Reac)
Reac_temp_indi <- Reac_temp_indi[primaryid %in% pids_diuretic_GDMT_indi]

Reac_temp_indi[pt %in% event_arr, event := "Cardiac Arrhythmias"]
Reac_temp_indi[pt %in% event_ren, event := "Acute renal failure"]
Reac_temp_indi[pt %in% event_syn, event := "Syncope"]

all_reac <- c("Cardiac Arrhythmias", "Acute renal failure", "Syncope")
Reac_temp_indi[is.na(event) | event %notin% all_reac, event := "Other event"]
report_arr <- unique(Reac[primaryid %in% pids_diuretic_GDMT_indi & pt %in% event_arr]$primaryid)
report_ren <- unique(Reac[primaryid %in% pids_diuretic_GDMT_indi & pt %in% event_ren]$primaryid)
report_syn <- unique(Reac[primaryid %in% pids_diuretic_GDMT_indi & pt %in% event_syn]$primaryid)

## Dataset with demographic, reaction and drug information ----
temp_indi <- left_join(Demo[primaryid %in% pids_diuretic_GDMT_indi], Reac_temp_indi, by = "primaryid")
temp_indi <- left_join(temp_indi, Drug_GDMT_indi, by = "primaryid")
temp_indi <- distinct(temp_indi)

## Percentage of events of interest over total reaction terms ----
all_event        <- unique(Reac[primaryid %in% pids_diuretic_GDMT_indi]$pt)
n_selected       <- sum(all_event %in% all_pt)
percent_selected <- (n_selected / length(all_event)) * 100

## Summary of selected events ----
df_events <- data.table(
  event = c("Cardiac Arrhythmias", "Acute renal failure", "Syncope"),
  N = c(length(report_arr),
        length(report_ren),
        length(report_syn)),
  Perc = c(round((length(report_arr) / length(report_selected_events_indi)) * 100, 2),
           round((length(report_ren) / length(report_selected_events_indi)) * 100, 2),
           round((length(report_syn) / length(report_selected_events_indi)) * 100, 2)))


## Frequencies for reactions ----
# Raw
pids_AR  <- unique(Reac[pt %in% event_arr]$primaryid)
pids_REN <- unique(Reac[pt %in% event_ren]$primaryid)
pids_SYN <- unique(Reac[pt %in% event_syn]$primaryid)

# Indi
pids_AR_indi  <- unique(Reac[primaryid %in% pids_indi & pt %in% event_arr]$primaryid)
pids_REN_indi <- unique(Reac[primaryid %in% pids_indi & pt %in% event_ren]$primaryid)
pids_SYN_indi <- unique(Reac[primaryid %in% pids_indi & pt %in% event_syn]$primaryid)

n_reports      <- uniqueN(Reac$primaryid)
n_reports_indi <- uniqueN(Reac[primaryid %in% pids_indi]$primaryid)

df_freq_ADR <- data.table(
  reaction = c("Cardiac Arrhythmias", "Acute renal failure", "Syncope"),
  N = c(length(pids_AR),
        length(pids_REN),
        length(pids_SYN)),
  Perc = c(round((length(pids_AR) / n_reports) * 100, 2),
           round((length(pids_REN) / n_reports) * 100, 2),
           round((length(pids_SYN) / n_reports) * 100, 2)),
  N_indi = c(length(pids_AR_indi),
             length(pids_REN_indi),
             length(pids_SYN_indi)),
  Perc_indi = c(round((length(pids_AR_indi) / n_reports_indi) * 100, 2),
                round((length(pids_REN_indi) / n_reports_indi) * 100, 2),
                round((length(pids_SYN_indi) / n_reports_indi) * 100, 2)))


## Type of reporter for each reaction reported ----
temp <- left_join(Demo[primaryid %in% pids_diuretic_GDMT_indi],
                  Reac_temp_indi[primaryid %in% pids_diuretic_GDMT_indi],
                  by = "primaryid")
temp <- left_join(temp, MedDRA, by = "pt")
reac_reporter <- temp[!is.na(hlgt) & !is.na(occp_cod), .(N = .N),by = .(hlgt, occp_cod)]

reac_reporter[, total := sum(N), by = hlgt]
reac_reporter[, perc := round((N / total) * 100, 2)]
setorder(reac_reporter, hlgt, -N)

## Type of reporter for each selected event ----
reac_reporter_arr <- temp[
  primaryid %in% report_arr & event == "Cardiac Arrhythmias",
  .(N = uniqueN(primaryid)),
  by = occp_cod][order(-N)][, percentage := round((N / length(report_arr)) * 100, 2)]

reac_reporter_ren <- temp[
  primaryid %in% report_ren & event == "Acute renal failure",
  .(N = uniqueN(primaryid)),
  by = occp_cod][order(-N)][, percentage := round((N / length(report_ren)) * 100, 2)]

reac_reporter_syn <- temp[
  primaryid %in% report_syn & event == "Syncope",
  .(N = uniqueN(primaryid)),
  by = occp_cod][order(-N)][, percentage := round((N / length(report_syn)) * 100, 2)]


## General event summary ----
events_summary <- data.table(
  measure = c("Reports with at least one selected event",
              "Selected PTs over all reported PTs"),
  N = c(length(report_selected_events_indi),
        n_selected),
  denominator = c(length(pids_diuretic_GDMT_indi),
                  length(all_event)),
  percentage = c(round((length(report_selected_events_indi) / length(pids_diuretic_GDMT_indi)) * 100, 2),
                 round(percent_selected, 2)))

write_xlsx(df_events, paste0(path, "percentage_of_events_of_interest.xlsx"))
write_xlsx(df_freq_ADR, paste0(path, "frequency_of_ADRs.xlsx"))
write_xlsx(reac_reporter, paste0(path, "type_of_reporter_per_reaction.xlsx"))
write_xlsx(reac_reporter_arr, paste0(path, "type_of_reporter_arrhythmias.xlsx"))
write_xlsx(reac_reporter_ren, paste0(path, "type_of_reporter_renal_failure.xlsx"))
write_xlsx(reac_reporter_syn, paste0(path, "type_of_reporter_syncope.xlsx"))
write_xlsx(events_summary, paste0(path, "events_summary.xlsx"))

# Full GDMT classes concomitant with diuretics, by reaction ----
# Raw
get_GDMT_classes_full <- function(data, pids, denominator) {
  data[
    primaryid %in% pids &
      substance %in% all_GDMT &
      !(substance %in% diuretics) &
      !is.na(class),
    .(N = uniqueN(primaryid)),
    by = class][
      order(-N)][
        , percentage := round((N / uniqueN(denominator)) * 100, 2)]
}

GDMT_classes_full_AR <- get_GDMT_classes_full(
  Drug_GDMT,
  intersect(pids_diuretic_GDMT, pids_AR),
  intersect(pids_diuretic_GDMT, pids_AR))

GDMT_classes_full_REN <- get_GDMT_classes_full(
  Drug_GDMT,
  intersect(pids_diuretic_GDMT, pids_REN),
  intersect(pids_diuretic_GDMT, pids_REN))

GDMT_classes_full_SYN <- get_GDMT_classes_full(
  Drug_GDMT,
  intersect(pids_diuretic_GDMT, pids_SYN),
  intersect(pids_diuretic_GDMT, pids_SYN))

# Indi
GDMT_classes_full_AR_indi <- get_GDMT_classes_full(
  Drug_GDMT_indi,
  intersect(pids_diuretic_GDMT_indi, pids_AR_indi),
  intersect(pids_diuretic_GDMT_indi, pids_AR_indi))

GDMT_classes_full_REN_indi <- get_GDMT_classes_full(
  Drug_GDMT_indi,
  intersect(pids_diuretic_GDMT_indi, pids_REN_indi),
  intersect(pids_diuretic_GDMT_indi, pids_REN_indi))

GDMT_classes_full_SYN_indi <- get_GDMT_classes_full(
  Drug_GDMT_indi,
  intersect(pids_diuretic_GDMT_indi, pids_SYN_indi),
  intersect(pids_diuretic_GDMT_indi, pids_SYN_indi))


write_xlsx(list("Cardiac Arrhythmias" = GDMT_classes_full_AR,
                "Acute renal failure" = GDMT_classes_full_REN,
                "Syncope" = GDMT_classes_full_SYN), 
           paste0(path, "diuretic_concomitant_GDMT_classes_full_by_reaction.xlsx"))

write_xlsx(list("Cardiac Arrhythmias" = GDMT_classes_full_AR_indi,
                "Acute renal failure" = GDMT_classes_full_REN_indi,
                "Syncope" = GDMT_classes_full_SYN_indi), 
           paste0 (path, "diuretic_concomitant_GDMT_classes_full_by_reaction_indi.xlsx"))

# Comorbidities ----------------------------------------------------------------
temp <- left_join(Drug[primaryid %in% pids_diuretic_GDMT_indi],
                  Indi[primaryid %in% pids_diuretic_GDMT_indi],
                  by = "primaryid")

other_indi <- temp[
  !substance %in% c(diuretics, all_GDMT) &
    !(indi_pt %in% indi_pt_selected) &
    !is.na(indi_pt),
  .(N = uniqueN(primaryid)),
  by = indi_pt][
    order(-N)][
      , percentage := round((N / uniqueN(temp$primaryid)) * 100, 2)]


# Distribution of reports over time -----
pids_diuretics_indi_ <- setdiff(pids_diuretics_indi, pids_diuretic_GDMT_indi)
pids_GDMT_indi_    <- setdiff(pids_GDMT_indi, pids_diuretic_GDMT_indi)

Reac_temp_indi2 <- copy(Reac)
Reac_temp_indi2 <- Reac_temp_indi2[primaryid %in% pids_indi]

Reac_temp_indi2[pt %in% event_arr, event := "Cardiac Arrhythmias"]
Reac_temp_indi2[pt %in% event_ren, event := "Acute renal failure"]
Reac_temp_indi2[pt %in% event_syn, event := "Syncope"]

# Trend datasets -----
get_year_counts <- function(pids, label) {
  copy(Demo)[primaryid %in% pids,
             year := as.numeric(substr(init_fda_dt, 1, 4))
  ][!is.na(year),
    .N, by = year
  ][, group := label]
}

get_trend_data <- function(event_name = NULL) {
  
  if (is.null(event_name)) {
    
    pids_diuretics <- pids_diuretics_indi_
    pids_GDMTs   <- pids_GDMT_indi_
    pids_combined  <- pids_diuretic_GDMT_indi
    
  } else {
    
    pids_diuretics <- unique(Reac_temp_indi2[
      event == event_name & primaryid %in% pids_diuretics_indi_]$primaryid)
    
    pids_GDMTs <- unique(Reac_temp_indi2[
      event == event_name & primaryid %in% pids_GDMT_indi_]$primaryid)
    
    pids_combined <- unique(Reac_temp_indi2[
      event == event_name & primaryid %in% pids_diuretic_GDMT_indi]$primaryid)
    
  }
  
  rbindlist(list(
    get_year_counts(pids_diuretics, "Diuretics"),
    get_year_counts(pids_GDMTs, "GDMTs"),
    get_year_counts(pids_combined, "Diuretics + GDMTs")
  ))
}

df_drugs    <- get_trend_data()
df_reac_AR  <- get_trend_data("Cardiac Arrhythmias")
df_reac_REN <- get_trend_data("Acute renal failure")
df_reac_SYN <- get_trend_data("Syncope")

# Plots -----
plot_trend <- function(data, colors, hide_x = FALSE) {
  
  approval_years <- data.frame(
    class = c("MRA", "Loop diuretics", "Beta-blockers", "ACE inhibitors", "ARBs", "SGLT-2 inhibitors", "ARNI"),
    year = c(1960, 1966, 1978, 1981, 1996, 2013, 2015)
  )
  approval_years <- approval_years[approval_years$year >= min(data$year) & approval_years$year <= max(data$year), ]
  
  p <- ggplot(data, aes(x = year, y = N, color = group)) +
    geom_vline(data = approval_years, aes(xintercept = year), inherit.aes = FALSE,
               linetype = "dashed", linewidth = .7, colour = "gray40") +
    geom_text(data = approval_years, aes(x = year, y = Inf, label = paste0(class, " (", year, ")")),
              inherit.aes = FALSE, angle = 90, vjust = 1.2, hjust = 1, size = 3.5, fontface = "bold") +
    geom_line(linewidth = 2.6) +
    geom_point(size = 4) +
    scale_color_manual(values = colors) +
    scale_y_continuous(limits = c(0, max(data$N) * 1.1)) +
    labs(
      x = "Year of submission",
      y = "N° of Reports",
      color = "Treatment group"
    ) +
    theme_bw(base_size = 15) +
    theme(
      text = element_text(size = 17, face = "bold"),
      axis.title = element_text(size = 17, face = "bold"),
      axis.text = element_text(size = 15, face = "bold", angle = 45),
      legend.text = element_text(size = 15, face = "bold"),
      legend.title = element_text(size = 17, face = "bold"),
      legend.position = c(0.2, 0.8),
      legend.background = element_rect(fill = "white", color = "black", size = 0.5),
      legend.key.size = unit(2, "lines")
    )
  
  if (hide_x) {
    
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank()
      )
    
  }
  
  p
}
plot_drugs <- plot_trend(df_drugs, c("Diuretics" = "green3",
                                     "GDMTs" = "seagreen3", "Diuretics + GDMTs" = "darkgreen"),
                         hide_x = TRUE)

plot_reac_AR <- plot_trend(df_reac_AR, c("Diuretics" = "blue3",
                                         "GDMTs" = "steelblue3", "Diuretics + GDMTs" = "grey50"),
                           hide_x = TRUE)

plot_reac_REN <- plot_trend(df_reac_REN, c("Diuretics" = "coral",
                                           "GDMTs" = "indianred3", "Diuretics + GDMTs" = "purple"))

plot_reac_SYN <- plot_trend(df_reac_SYN, c("Diuretics" = "goldenrod1",
                                           "GDMTs" = "goldenrod4", "Diuretics + GDMTs" = "orangered"))

# Combined trend plot -----
combined_trend_plot <- (plot_drugs + plot_reac_AR) / (plot_reac_REN + plot_reac_SYN) +
  plot_annotation(
    title = "Distribution of reports over time",
    subtitle = "Full available period, conditioned by indication (Cardiac failure)",
    theme = theme(
      plot.title = element_text(size = 22, face = "bold"),
      plot.subtitle = element_text(size = 16, face = "bold")))

ggsave(paste0(path, "combined_trends.tiff"),
       plot = combined_trend_plot,
       width = 20, height = 16, units = "in",
       dpi = 1000, compression = "lzw")

# Summary of reports over time -----
trend_summary <- rbindlist(list(
  copy(df_drugs)[, event := "All ADRs"],
  copy(df_reac_AR)[, event := "Cardiac Arrhythmias"],
  copy(df_reac_REN)[, event := "Acute renal failure"],
  copy(df_reac_SYN)[, event := "Syncope"]
), use.names = TRUE, fill = TRUE)

setcolorder(trend_summary, c("event", "group", "year", "N"))
setorder(trend_summary, event, group, year)

# Total reports by treatment group and event -----
trend_totals <- trend_summary[
  , .(N = sum(N)),
  by = .(event, group)][
    , percentage := round((N / sum(N)) * 100, 2),
    by = event][
      order(event, -N)]

# Most frequently reported reactions -----
event_reported <- merge(
  Reac[primaryid %in% pids_diuretic_GDMT_indi,
       .(N = uniqueN(primaryid)),
       by = pt],
  MedDRA[, .(pt, hlt)],
  by = "pt",
  all.x = TRUE)[
    !is.na(hlt)][
      , N_hlt := sum(N),
      by = hlt][
        order(-N_hlt, -N)][
          , .(pt = paste0(pt, " (N = ", N, ")", collapse = "; ")),
          by = .(hlt, N_hlt)][
            order(-N_hlt)]

# Most frequently reported individual reactions -----
event_reported_pt <- Reac[
  primaryid %in% pids_diuretic_GDMT_indi,
  .(N = uniqueN(primaryid)),
  by = pt][
    order(-N)][
      , percentage := round((N / uniqueN(pids_diuretic_GDMT_indi)) * 100, 2)]

write_xlsx(other_indi, paste0(path, "Most_reported_comorbidities.xlsx"))
write_xlsx(trend_summary, paste0(path, "reports_over_time.xlsx"))
write_xlsx(trend_totals, paste0(path, "reports_over_time_summary.xlsx"))
write_xlsx(event_reported, paste0(path, "Most_reported_reac.xlsx"))
write_xlsx(event_reported_pt, paste0(path, "Most_reported_reac_PT.xlsx"))

# ----------------------------------------------------------------

# Clear & save environment -----
rm(Drug_GDMT, Drug_GDMT_indi, NGAL, Outc, Reac_temp_indi, Reac_temp_indi2, 
   SMQs_selected_A_broad, SMQs_selected_A_narrow, Ther, albume, all_drugs, 
   all_event, all_reac, cerebro, class_counts, class_counts_indi, combined_trend_plot, 
   combo_counts, combo_counts_indi, concomitant, concomitant_indi, concomitant_perc, 
   concomitant_perc_indi, crystal, df_drugs, df_events, df_freq_ADR, df_reac_AR, 
   df_reac_REN, df_reac_SYN, ecg, event_reported, event_reported_pt, events_summary, 
   exclude_arr, filter_SMQ_arr, filter_SMQ_ren, flowchart, foetal, heart_rate, kidney, 
   loc, n_reports, n_reports_indi, n_selected, neonatal, nutrition, other_indi, 
   palpitations, perc_arr, perc_raw, perc_ren, perc_syn, percent_selected, pids_AR, 
   pids_AR_indi, pids_REN, pids_REN_indi, pids_SYN, pids_SYN_indi, 
   pids_diuretic_ace_one_class, pids_diuretic_ace_one_class_indi, 
   pids_diuretic_arb_one_class, pids_diuretic_arb_one_class_indi, 
   pids_diuretic_arni_one_class, pids_diuretic_arni_one_class_indi, 
   pids_diuretic_beta_one_class, pids_diuretic_beta_one_class_indi, 
   pids_diuretic_mra_one_class, pids_diuretic_mra_one_class_indi, pids_diuretic_GDMT, 
   pids_diuretic_GDMT_indi, pids_diuretic_sglt_one_class, 
   pids_diuretic_sglt_one_class_indi, pids_diuretics_indi, pids_diuretics_indi_, 
   pids_one_class, pids_one_class_indi, pids_GDMT_indi_, GDMT_class_list, 
   GDMT_class_list_indi, GDMT_combinations, GDMT_combinations_indi, plot_drugs, 
   plot_reac_AR, plot_reac_REN, plot_reac_SYN, prote, reac_reporter, 
   reac_reporter_arr, reac_reporter_ren, reac_reporter_syn, rebound, reflex, 
   report_arr, report_ren, report_selected_events, report_selected_events_indi, 
   report_syn, sinus, smqs_broad, smqs_narrow, summary_counts, summary_counts_indi, 
   svg_txt, sync, temp, temp_indi, total_concomitant_primaryid, 
   total_concomitant_primaryid_indi, total_patients, total_patients_indi, transplant, 
   trend_summary, trend_totals)

# Save environment
save.image(paste0(path, "Descriptive_analysis.RData"))
