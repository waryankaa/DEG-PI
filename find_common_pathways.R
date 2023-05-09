
library(WebGestaltR)

interestGeneFile = "/Users/abbywaryanka/seniorComp/Streamlit/webgestaltData.rnk"

enrichrestults <- WebGestaltR(
  enrichMethod = "GSEA",
  organism = "hsapiens",
  enrichDatabase = "pathway_KEGG",
  interestGeneFile = interestGeneFile,
  interestGeneType = "genesymbol",
  collapseMethod = "mean",
  minNum = 3,
  maxNum = 1000,
  sigMethod = "top",
  fdrMethod = "BH",
  fdrThr = 0.05,
  topThr = 10,
  reportNum = 20,
  perNum = 1000,
  gseaP = 1,
  isOutput = TRUE,
  outputDirectory = "/Users/abbywaryanka/seniorComp/Streamlit",
  projectName = "test",
  dagColor = "continuous",
  saveRawGseaResult = TRUE,
  gseaPlotFormat = c("png", "svg"),
  setCoverNum = 10,
  neighborNum = 10,
  highlightType = "Seeds",
  highlightSeedNum = 10,
  nThreads = 1,
  hostName = "https://www.webgestalt.org/"
)

 

