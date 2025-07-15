

# ATC Drug Reference Data

ATCCode <- c("N06AA09", "N06AB04", "N06AB05", "N06AX16", "N06AX11", "N06AB10", "N06AB06",
             "N06AB03", "N06AX21", "N06AA10", "N06AX12", "N06AA04", "N06AB08",
             "N06AX05", "N06AF04", "N06AX26", "N06AX22", "N06AA12", "N06AA02",
             "N06AA21", "N06AA16", "N06AF03", "N06AG02", "N06AX03", "N06AX25", "N05AN01",
             "N06AX23")
DrugName <- c("Amitriptyline", "Citalopram", "Paroxetine", "Venlafaxine", "Mirtazapine", "Escitalopram",
              "Sertraline","Fluoxetine", "Duloxetine", "Nortriptyline", "Bupropion", "Clomipramine",
              "Fluvoaxmine", "Trazodone", "Tranylcypromine", "Vortioxetine", "Agomelatine", "Doxepine",
              "Imipramine", "Maprotiline", "Dosulepin", "Phenelzine", "Moclobemide", "Mianserin",
              "HypericiHerba", "Lithium", "Desvenlafaxine")
ATCClass <- c("AIVc", "AIIa", "AIIa", "AIIb","NaSSA", "AIIa",
              "AIIa", "AIIa", "AIIb", "AIVc", "AIIc", "AIVc",
              "AIIa", "AIVa", "AIa", "AIVa", "B", "AIVc",
              "AIVc", "AIVb", "AIVc", "AIa", "AIc", "AIVb",
              "HypericiHerba", "Lithium", "AIIb"
)
            
# Create the data frame
drug_ref <- data.frame(
  ATCCode = ATCCode,
  Drug = DrugName,
  DrugName = paste0(ATCClass, ":", DrugName),
  DrugClass = ATCClass,
  stringsAsFactors = FALSE
)

