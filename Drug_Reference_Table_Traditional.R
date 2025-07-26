

# ATC Drug Reference Data

ATCCode <- c("N06AA09", "N06AB04", "N06AB05", "N06AX16", "N06AX11",
             "N06AB10", "N06AB06","N06AB03", "N06AX21", "N06AX23")
DrugName <- c("Amitriptyline", "Citalopram", "Paroxetine", "Venlafaxine", "Mirtazapine", 
              "Escitalopram", "Sertraline", "Fluoxetine", "Duloxetine", "Desvenlafaxine")
ATCClass <- c("TCA", "SSRI", "SSRI", "SNRI", "TeCA",
              "SSRI", "SSRI", "SSRI", "SNRI", "SNRI"
)

# Create the data frame
drug_ref <- data.frame(
  ATCCode = ATCCode,
  Drug = DrugName,
  DrugName = paste0(ATCClass, ":", DrugName),
  DrugClass = ATCClass,
  stringsAsFactors = FALSE
)

