#load the libaries required
library(shiny)
library(readr)
library(tidyverse)
library(seqRFLP)
library(ggplot2)
library(Rcpp)
library(dplyr)
library(DT)
library(shinythemes)

#read in the datasets
MyTable <- read_csv("data/pfam_mitominer_tissue.csv")
MyOrtho <- read_csv("data/orthology_mitominer_domcomb_tissue.csv")
MyPre <- read_csv("data/presequence_mitominer_domcomb_tissue.csv")
MyGO <- read_csv("data/go_mitominer.csv")
MyGene <- read_csv("data/gene.csv")
MyOntology <- read_csv("data/ontology.csv")
MyDomOnt <- read_csv("data/pfam2go.csv")

#extract colum informations
col1 <-  MyTable$species
col2 <-  MyTable$domname
col11 <- MyTable$tp
col12 <- MyTable$mf_pred
col30 <- MyTable$HumanProteinAtlas
col31 <- MyTable$GFP
col32 <- MyTable$MassSpecStudies
col33 <- MyTable$tissue

col3 <- MyOrtho$og
col4 <- MyOrtho$species
col5 <- MyOrtho$mito
col20 <- MyOrtho$tp
col21 <- MyOrtho$mf_pred
col22 <- MyOrtho$GO
col23 <- MyOrtho$HumanProteinAtlas
col24 <- MyOrtho$GFP
col25 <- MyOrtho$MassSpecStudies
col34 <- MyOrtho$tissue

col10 <- MyPre$species
col7 <- MyPre$tp
col8 <- MyPre$rc
col9 <- MyPre$mf
col26 <- MyPre$OG_present
col27 <- MyPre$HumanProteinAtlas
col28 <- MyPre$GFP
col29 <- MyPre$MassSpecStudies
col35 <- MyPre$tissue

col15 <-  MyGO$species
col16 <-  MyGO$goid
col18 <-  MyGO$tp
col19 <-  MyGO$mf_pred
col36 <-  MyGO$HumanProteinAtlas
col37 <-  MyGO$GFP
col38 <-  MyGO$MassSpecStudies

###########################################################################################################################
########################## part 2 Define ui logic #########################################################################
###########################################################################################################################

ui <- fluidPage(
  theme = shinytheme("sandstone"),
  
  titlePanel("Metazoan Mitochondrial Proteome Database (MMPdb)"),
  
  #Define the tab panels
  mainPanel(
    tabsetPanel(
      type = "tabs",
      tabPanel(
        "About",
        br(),
        p(
          "MMPdb is a database for facilitating comparative analysis of experimentally-characterized mitochondrial proteomes of bilaterian animals. Publicly 
           available mitochondrial proteomes were downloaded for four animals - Homo sapiens [1], Mus musculus [1], Caenorhabditis elegans [2], Drosophila
           melanogaster [3] and two outgroups - Acanthamoeba castellanii [4] and Saccharomyces cerevisiae [5]. Each species is denoted by a four letter 
           abbreviation listed below."
        ),
        p("acas: Acanthamoeba castellanii"),
        p("scer: Saccharomyces cerevisiae"),
        p("hsap: Homo sapiens"),
        p("mmus: Mus musculus"),
        p("cele: Caenorhabditis elegans"),
        p("dmel: Drosophila melanogaster"),
        br(),
        p(
          "In this tabset, there are instructions regarding the navigation and usage of the database."
        ),
        strong("Usage of the Orthology Tabset"),
        p(
          "Proteinortho v5.16b [6] was used to identify Orthologous Groups (OGs) in the four bilaterian species and the two outgroups.
          Some OGs include both mitochondrial and non-mitochondrial proteins. Each OG is assigned a unique identifier, which is also the
          primary query for this tabset. The OG number for a protein of interest can be identified as follows. Use the search tool below 
          to find the OG ID for querying the Orthology Tabset. "
        ),
        p(
          "Different search identifiers need to be used to search for proteins from different species."
        ),
        p(
          "For human, mouse, A.castellanii and D.melanogaster, enter the Primary accession number from https://www.uniprot.org/. For example, 
          if you want to extract orthologs of human fumarate hydratase, then enter its Uniprot identifier 'P07954' in the search tool below.
          The search tool will display the OG number 'OG1339' the protein belongs to. Enter that OG number 'OG1339' in the search bar in the 
          Orthology tabset to explore the fumarate hydratase OG."
        ),
        p(
          "For the nematode C. elegans, enter the Protein Wormbase ID from http://intermine.wormbase.org/. For example: CE11580 for fumarate hydratase."
        ),
        p(
          "For yeast, enter the ID from the Saccharomyces Cerevisiae Genome database at https://www.yeastgenome.org/. For example: YPL262W for fumarate hydratase."
        ),
        p(
          "Once you have the OG ID, use it as the search query in the Orthology tabset. (e.g. enter OG1339) in the search tool in the Orthology
          tabset to fetch the OG.)"
        ),
        br(),
        br(),
        DT::dataTableOutput(outputId = "view6"),
        br(),
        strong("Usage of the MTS Tabset"),
        p(
          "N-terminal mitochondrial targeting presequences (MTS) were identified using TargetP [7] and MitoFates [8]. This tabset allows for 
          exploration of the MTS-prediction results for all the animal and outgroup proteins."
        ),
        br(),
        strong("Usage of the Domain Tabset"),
        p(
          "Protein domain analysis was performed using PFAM libraries and pfam_scan.pl [9]. In tabset, a user can explore and download proteins
           containing a PFAM domain of interest. For example, to fetch data for all proteins containing the mitochondrial carrier protein domain,
           enter the domain name 'Mito_carr' in the search bar in the Domain tabset."
        ),
        p("Below we provide the the pfam2go mappings which might aid in selection of protein domains based
           on the GO term they are mapped to. It is important to note that if a domain is not present in the pfam2go mapping file, it does not mean
           that it is not identified in any of the six eukaryotes in the database."),
        br(),
        br(),
        DT::dataTableOutput(outputId = "view8"),
        br(),
        strong("Gene Ontology Tabset"),
        p(
          "Gene Ontology analysis was performed using Pannzer2.0 [10]. GO terms can be searched using the search tool below."
        ),
        br(),
        br(),
        DT::dataTableOutput(outputId = "view7"),
        br(),
        br(),
        p(
          "Complete database of results, mitochondrial proteomes can be downloaded from : V., Lavrov, D. & Kandoi, G. Data for Metazoan Mitochondrial Proteome database (MMPdb). (2019). Available at: osf.io/gfyq9.  ."
        ),
        br(),
        strong("Additional information for mammalian mitochondrial proteomes"),
        p(
          "The most extensively studied mitochondrial proteomes from animals belong to mammals. There is additional information for human and mouse mitochondrial
           proteomes in this database:"
        ),
        p(
          "Mitochondrial localization evidence from MitoMiner: 1] Mass-spectrometry and 2] GFP-analysis. This data is available for both mammals."
        ),
        p(
          "1] Mitochondrial localization evidence and 2] Tissue-enrichment evidence from Human Protein Atlas (https://www.proteinatlas.org/humanproteome/tissue/tissue+specific).
           This data is available for just human proteins."
        ),
        br(),
        strong("References"),
        p(
          "[1] Calvo, Sarah E., Karl R. Clauser, and Vamsi K. Mootha. MitoCarta2. 0: an updated inventory of mammalian mitochondrial proteins. Nucleic acids research 44, no. D1 (2015): D1251-D1257."
        ),
        p(
          "[2] Li, J., Cai, T., Wu, P., Cui, Z., Chen, X., Hou, J., Xie, Z., Xue, P., Shi, L., Liu, P. and Yates III, J.R., 2009. Proteomic analysis of mitochondria from Caenorhabditis elegans. Proteomics, 9(19), pp.4539-4553."
        ),
        p(
          "[3] Hu, Y., Comjean, A., Perkins, L.A., Perrimon, N. and Mohr, S.E., 2015. GLAD: an online database of gene list annotation for Drosophila. Journal of genomics, 3, p.75."
        ),
        p(
          "[4] Gawryluk, R.M., Chisholm, K.A., Pinto, D.M. and Gray, M.W., 2014. Compositional complexity of the mitochondrial proteome of a unicellular eukaryote (Acanthamoeba castellanii, supergroup Amoebozoa) rivals that of animals, fungi, and plants. Journal of proteomics, 109, pp.400-416."
        ),
        p(
          "[5] Cherry, J.M., Adler, C., Ball, C., Chervitz, S.A., Dwight, S.S., Hester, E.T., Jia, Y., Juvik, G., Roe, T., Schroeder, M. and Weng, S., 1998. SGD: Saccharomyces genome database. Nucleic acids research, 26(1), pp.73-79."
        ),
        p(
          "[6] Lechner, M., Findei, S., Steiner, L., Marz, M., Stadler, P.F. and Prohaska, S.J., 2011. Proteinortho: detection of (co-) orthologs in large-scale analysis. BMC bioinformatics, 12(1), p.124."
        ),
        p(
          "[7] Emanuelsson, O., Brunak, S., Von Heijne, G. and Nielsen, H., 2007. Locating proteins in the cell using TargetP, SignalP and related tools. Nature protocols, 2(4), p.953."
        ),
        p(
          "[8] Fukasawa, Y., Tsuji, J., Fu, S.C., Tomii, K., Horton, P. and Imai, K., 2015. MitoFates: improved prediction of mitochondrial targeting sequences and their cleavage sites. Molecular & Cellular Proteomics, 14(4), pp.1113-1126."
        ),
        p(
          "[9] Bateman, A., Coin, L., Durbin, R., Finn, R.D., Hollich, V., Griffiths-Jones, S., Khanna, A., Marshall, M., Moxon, S., Sonnhammer, E.L. and Studholme, D.J., 2004. The Pfam protein families database. Nucleic acids research, 32(suppl_1), pp.D138-D141."
        ),
        p(
          "[10] Medlar, A.J., Toronen, P., Zosa, E. and Holm, L., PANNZER 2: Annotate a complete proteome in minutes!. Nucl. Acids Res, 43, pp.W24-W29."
        )
        ),
      
      ######################################################################################################################################333
      
      tabPanel(
        "Orthology",
        sidebarPanel(
          helpText("All field are required fields."),
          checkboxGroupInput(
            "species",
            label = "Select species:",
            choices = sort(unique(col4)),
            selected = sort(unique(col4))
          ),
          helpText(
            "Each species is denoted by a four letter abbreviation. Check the 'About tab' for details about the abbreviations."
          ),
          checkboxGroupInput(
            "mitolabel",
            label = "Mitochondrial / Non-mitochondrial proteins:",
            c("Mitochondrial" = "Y", "Non-mitochondrial" = "N"),
            selected = c("Y")
          ),
          helpText(
            "Some Orthologous Groups (OGs) include both mitochondrial and non-mitochondrial proteins. To fetch only mitochondrial proteins, select 'Mitochondrial',
            and to select both mitochondrial and non-mitochondrial, select both options"
          ),
          checkboxGroupInput(
            "tp_3",
            label = "TargetP prediction",
            choices = sort(unique(col20)),
            selected = sort(unique(col20))
          ),
          helpText("MTS:MTS detected by TargetP, No MTS:No MTS detected by TargetP"),
          checkboxGroupInput(
            "mf_pred_3",
            label = "MitoFates prediction",
            choices = sort(unique(col21)),
            selected = sort(unique(col21))
          ),
          helpText(
            "MTS:MTS detected by MitoFates, No MTS:No MTS detected by MitoFates"
          ),
          checkboxGroupInput(
            "GFP_3",
            label = "GFP evidence (Mammals only)",
            choices = sort(unique(col24)),
            selected = sort(unique(col24))
          ),
          helpText(
            "1: GFP evidence of mitochondrial localization, 0: No GFP evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "HumanProteinAtlas_3",
            label = "Human Protein Atlas (HPA) evidence (Human only)",
            choices = sort(unique(col23)),
            selected = sort(unique(col23))
          ),
          helpText(
            "TRUE: HPA evidence of mitochondrial localization, FALSE: No HPA evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "MassSpecStudies_3",
            label = "Mass Spectrometry evidence (Mammals only)",
            choices = sort(unique(col25)),
            selected = sort(unique(col25))
          ),
          helpText(
            "TRUE: MS evidence of mitochondrial localization, FALSE: No MS evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "tissue_3",
            label = "Tissue-enriched (HPA)(Human only)",
            choices = sort(unique(col34)),
            selected = sort(unique(col34))
          ),
          helpText(
            "Tissue-enrichment evidence from Human-Protein atlas."
          ),
          downloadButton(outputId = "download_data2", label = "Download sequences"),
          helpText("Download selected sequences is a FASTA format"),
          br(),
          strong("Key for table"),
          br(),
          strong("og"),
          p("Unique Orthologous Group (OG) identifier."),
          strong("mito"),
          p("Mitochondrial/Non-mitochondrial protein"),
          strong("tp"),
          p("TargetP prediction"),
          strong("rc"),
          p("Reliability Class for TargetP prediction"),
          strong("mf_pred"),
          p("MitoFates prediction"),
          strong("mf_prob"),
          p("MitoFates prediction probability"),
          strong("mppsite"),
          p("Mitochondrial Processing Peptidase cleavage site"),
          strong("GO"),
          p("Mitochondrial localization evidence from MitoMiner (Gene Ontology)"),
          strong("HumanProteinAtlas"),
          p("Mitochondrial localization evidence from MitoMiner (Human Protein Atlas)"),
          strong("GFP"),
          p("Mitochondrial localization evidence from MitoMiner (GFP analysis)"),
          strong("MassSpecStudies"),
          p("Mitochondrial localization evidence from MitoMiner (Mass spectrometry studies)"),
          strong("domainComb"),
          p("Protein domain combination"),
          strong("tissue"),
          p("Tissue-enrichent evidence from MitoMiner"),
          br()
        ),
        
        mainPanel(
          selectInput(
            inputId = "og",
            label = "Select OG:",
            choices = sort(unique(col3)),
            multiple = TRUE
          ),
          plotOutput(outputId = "plot2", width = "100%"),
          br(),
          br(),
          DT::dataTableOutput(outputId = "view2")
        )
        
      ),
      
 
      ###############################################################################################################################      
      tabPanel(
        "MTS",
        sidebarPanel(
          helpText("All field are required fields."),
          checkboxGroupInput(
            "species_3",
            label = "Select species:",
            choices = sort(unique(col10)),
            selected = sort(unique(col10))
          ),
          helpText(
            "Each species is denoted by a four letter abbreviation. Check the About tab for details about the abbreviations."
          ),
          checkboxGroupInput(
            "mito",
            label = "Mitochondrial / Non-mitochondrial proteins:",
            c("Mitochondrial" = "Y", "Non-mitochondrial" = "N"),
            selected = c("Y")
          ),
          helpText(
            "Mitochondrial: Mitochondrial proteins; Non-mitochondrial: Non-mitochondrial proteins"
          ),
          checkboxGroupInput(
            "OG_present",
            label = "Present/absent in OG",
            choices = sort(unique(col26)),
            selected = sort(unique(col26))
          ),
          helpText(
            "Is the protein present in any OG, i.e. does it have an ortholog in any other species."
          ),
          checkboxGroupInput(
            "tp",
            label = "TargetP prediction",
            choices = sort(unique(col7)),
            selected = sort(unique(col7))
          ),
          helpText(
            "TargetP prediction result, M:Mitochondrial,S:Secreted,_:No prediction"
          ),
          checkboxGroupInput(
            "rc",
            label = "TargetP RC",
            choices = sort(unique(col8)),
            selected = sort(unique(col8))
          ),
          helpText(
            "Reliability Class of TargetP prediction. This ranges from 1 (strongest prediction) to 5 (weakest prediction)"
          ),
          checkboxGroupInput(
            "mf",
            label = "MitoFates prediction",
            choices = sort(unique(col9)),
            selected = sort(unique(col9))
          ),
          helpText(
            "Prediction of MitoFates, No_mitochondrial_presequence: No MTS detected, Possessing_mitochondrial_presequence: MTS detected"
          ),
          
          checkboxGroupInput(
            "GFP",
            label = "GFP evidence (Mammals only)",
            choices = sort(unique(col28)),
            selected = sort(unique(col28))
          ),
          helpText(
            "1: GFP evidence of mitochondrial localization, 0: No GFP evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "HumanProteinAtlas",
            label = "Human Protein Atlas (HPA) evidence (Human only)",
            choices = sort(unique(col27)),
            selected = sort(unique(col27))
          ),
          helpText(
            "TRUE: HPA evidence of mitochondrial localization, FALSE: No HPA evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "MassSpecStudies",
            label = "Mass Spectrometry evidence (Mammals only)",
            choices = sort(unique(col29)),
            selected = sort(unique(col29))
          ),
          helpText(
            "0-15: Number of MS studies with evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "tissue",
            label = "Tissue-enriched (HPA) (Human only)",
            choices = sort(unique(col35)),
            selected = sort(unique(col35))
          ),
          helpText(
            "Tissue-enrichment evidence from Human-Protein atlas"
          ),
          downloadButton(outputId = "download_data_3", label = "Download sequences"),
          helpText("Download selected sequences is a FASTA format"),
          br(),
          br(),
          downloadButton(outputId = "download_data_4", label = "Download table"),
          helpText("Download resulting table for selected sequences in .csv format"),
          br(),
          strong("Key for table"),
          br(),
          strong("og"),
          p("Outgroup ID (gene name)"),
          strong("mito"),
          p("Mitochondrial/Non-mitochondrial protein"),
          strong("tp"),
          p("TargetP prediction"),
          strong("rc"),
          p("Reliability Class for TargetP prediction"),
          strong("mf"),
          p("MitoFates prediction"),
          strong("mf_prob"),
          p("MitoFates prediction probability"),
          strong("mppsite"),
          p("Mitochondrial Processing Peptidase cleavage site"),
          p("Mitochondrial localization evidence from MitoMiner (Gene Ontology)"),
          strong("HumanProteinAtlas"),
          p("Mitochondrial localization evidence from MitoMiner (Human Protein Atlas)"),
          strong("GFP"),
          p("Mitochondrial localization evidence from MitoMiner (GFP analysis)"),
          strong("MassSpecStudies"),
          p("Mitochondrial localization evidence from MitoMiner (Mass spectrometry studies)"),
          strong("domainComb"),
          p("Protein domain combination"),
          strong("tissue"),
          p("Tissue-enrichent evidence from MitoMiner"),
          br()
        ),
        mainPanel(
          plotOutput(outputId = "plot3"),
          br(),
          br(),
          DT::dataTableOutput(outputId = "view3")
        )
        
      ),
 
      ###############################################################################################################################           
      tabPanel(
        "Domain",
        sidebarPanel(
          helpText("All field are required fields."),
          checkboxGroupInput(
            "species_2",
            label = "Select species:",
            choices = sort(unique(col1)),
            selected = sort(unique(col1))
          ),
          helpText(
            "Each species is denoted by a four letter abbreviation. Check the About tab for details about the abbreviations."
          ),
          checkboxGroupInput(
            "mito_2",
            label = "Mitochondrial / Non-mitochondrial proteins:",
            c("Mitochondrial" = "Y", "Non-mitochondrial" = "N"),
            selected = c("Y", "N")
          ),
          helpText(
            "Mitochondrial: Mitochondrial proteins; Non-mitochondrial: Non-mitochondrial proteins"
          ),
          checkboxGroupInput(
            "tp_2",
            label = "TargetP:",
            choices = sort(unique(col11)),
            selected = sort(unique(col11))
          ),
          helpText(
            "TargetP prediction result, M:Mitochondrial,S:Secreted,_:No prediction"
          ),
          checkboxGroupInput(
            "mf_pred_2",
            label = "MitoFates:",
            choices = sort(unique(col12)),
            selected = sort(unique(col12))
          ),
          helpText(
            "Prediction of MitoFates, No_mitochondrial_presequence: No MTS detected, Possessing_mitochondrial_presequence: MTS detected"
          ),
          
          checkboxGroupInput(
            "GFP_2",
            label = "GFP evidence (Mammals only)",
            choices = sort(unique(col31)),
            selected = sort(unique(col31))
          ),
          helpText(
            "1: GFP evidence of mitochondrial localization, 0: No GFP evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "HumanProteinAtlas_2",
            label = "Human Protein Atlas (HPA) evidence (Human only)",
            choices = sort(unique(col30)),
            selected = sort(unique(col30))
          ),
          helpText(
            "TRUE: HPA evidence of mitochondrial localization, FALSE: No HPA evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "MassSpecStudies_2",
            label = "Mass Spectrometry evidence (Mammals only)",
            choices = sort(unique(col32)),
            selected = sort(unique(col32))
          ),
          helpText(
            "0-15: Number of MS studies with evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "tissue_2",
            label = "Tissue-enriched (HPA) (Human only)",
            choices = sort(unique(col33)),
            selected = sort(unique(col33))
          ),
          helpText(
            "Tissue-enrichment evidence from Human-Protein atlas"
          ),
          downloadButton(outputId = "download_data", label = "Download sequences"),
          helpText("Download selected sequences is a FASTA format"),
          br(),
          br(),
          downloadButton(outputId = "download_data_6", label = "Download table"),
          helpText("Download resulting table for selected sequences in .csv format"),
          br(),
          strong("Key for table"),
          br(),
          strong("og"),
          p("Outgroup ID (gene name)"),
          strong("mito"),
          p("Mitochondrial/Non-mitochondrial protein"),
          strong("domid"),
          p("PFAM domain Id"),
          strong("domname"),
          p("PFAM domain name"),
          strong("comb"),
          p("Protein domain combination"),
          strong("tp"),
          p("TargetP prediction"),
          strong("rc"),
          p("Reliability Class for TargetP prediction"),
          strong("mf_pred"),
          p("MitoFates prediction"),
          strong("mf_prob"),
          p("MitoFates prediction probability"),
          strong("mppsite"),
          p("Mitochondrial Processing Peptidase cleavage site"),
          p("Mitochondrial localization evidence from MitoMiner (Gene Ontology)"),
          strong("HumanProteinAtlas"),
          p("Mitochondrial localization evidence from MitoMiner (Human Protein Atlas)"),
          strong("GFP"),
          p("Mitochondrial localization evidence from MitoMiner (GFP analysis)"),
          strong("MassSpecStudies"),
          p("Mitochondrial localization evidence from MitoMiner (Mass spectrometry studies)"),
          strong("tissue"),
          p("Tissue-enrichent evidence from MitoMiner"),
          br()
          
        ),
        mainPanel(
          selectInput(
            inputId = "domname",
            label = "Select domain:",
            choices = sort(unique(col2)),
            multiple = TRUE
          ),
          plotOutput(outputId = "plot"),
          br(),
          br(),
          plotOutput(outputId = "plot5"),
          br(),
          br(),
          DT::dataTableOutput(outputId = "view")
        )
        
      ),

      ###############################################################################################################################      
      tabPanel(
        "GO analysis",
        sidebarPanel(
          helpText("All field are required fields."),
          checkboxGroupInput(
            "species_5",
            label = "Select species:",
            choices = sort(unique(col15)),
            selected = sort(unique(col15))
          ),
          helpText(
            "Each species is denoted by a four letter abbreviation. Check the About tab for details about the abbreviations."
          ),
          checkboxGroupInput(
            "mito_5",
            label = "Mitochondrial / Non-mitochondrial proteins:",
            c("Mitochondrial" = "Y", "Non-mitochondrial" = "N"),
            selected = c("Y", "N")
          ),
          helpText(
            "Mitochondrial: Mitochondrial proteins; Non-mitochondrial: Non-mitochondrial proteins"
          ),
          checkboxGroupInput(
            "tp_5",
            label = "TargetP:",
            choices = sort(unique(col18)),
            selected = sort(unique(col18))
          ),
          helpText(
            "TargetP prediction result, M:Mitochondrial,S:Secreted,_:No prediction"
          ),
          checkboxGroupInput(
            "mf_pred_5",
            label = "MitoFates:",
            choices = sort(unique(col19)),
            selected = sort(unique(col19))
          ),
          helpText(
            "Prediction of MitoFates, No_mitochondrial_presequence: No MTS detected, Possessing_mitochondrial_presequence: MTS detected"
          ),
          checkboxGroupInput(
            "GFP_5",
            label = "GFP evidence (Mammals only)",
            choices = sort(unique(col37)),
            selected = sort(unique(col37))
          ),
          helpText(
            "1: GFP evidence of mitochondrial localization, 0: No GFP evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "HumanProteinAtlas_5",
            label = "Human Protein Atlas (HPA) evidence (Human only)",
            choices = sort(unique(col36)),
            selected = sort(unique(col36))
          ),
          helpText(
            "TRUE: HPA evidence of mitochondrial localization, FALSE: No HPA evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "MassSpecStudies_5",
            label = "Mass Spectrometry evidence (Mammals only)",
            choices = sort(unique(col38)),
            selected = sort(unique(col38))
          ),
          helpText(
            "0-15: Number of MS studies with evidence, nd: No Data"
          ),
          downloadButton(outputId = "download_data_8", label = "Download sequences"),
          helpText("Download selected sequences is a FASTA format"),
          br(),
          br(),
          downloadButton(outputId = "download_data_7", label = "Download table"),
          helpText("Download resulting table for selected sequences in .csv format"),
          br(),
          strong("Key for table"),
          br(),
          strong("og"),
          p("Outgroup ID (gene name)"),
          strong("mito"),
          p("Mitochondrial/Non-mitochondrial protein"),
          strong("goid"),
          p("Gene Ontology term ID"),
          strong("tp"),
          p("TargetP prediction"),
          strong("rc"),
          p("Reliability Class for TargetP prediction"),
          strong("mf"),
          p("MitoFates prediction"),
          strong("mf_prob"),
          p("MitoFates prediction probability"),
          strong("mppsite"),
          p("Mitochondrial Processing Peptidase cleavage site"),
          p("Mitochondrial localization evidence from MitoMiner (Gene Ontology)"),
          strong("HumanProteinAtlas"),
          p("Mitochondrial localization evidence from MitoMiner (Human Protein Atlas)"),
          strong("GFP"),
          p("Mitochondrial localization evidence from MitoMiner (GFP analysis)"),
          strong("MassSpecStudies"),
          p("Mitochondrial localization evidence from MitoMiner (Mass spectrometry studies)"),
          strong("tissue"),
          p("Tissue-enrichent evidence from MitoMiner"),
          br()
          
        ),
        mainPanel(
          selectInput(
            inputId = "goid",
            label = "Select GO term:",
            choices = sort(unique(col16)),
            multiple = TRUE
          ),
          plotOutput(outputId = "plot6"),
          br(),
          br(),
          DT::dataTableOutput(outputId = "view5")
        )
        
      )
        )
      )
)

###############################################################################################################################
########################## part 2 Define server logic #########################################################################
###############################################################################################################################

server <- function(input, output) {
  output$view6 <- DT::renderDataTable({
    DT::datatable(
      data = MyGene,
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })
  
  output$view7 <- DT::renderDataTable({
    DT::datatable(
      data = MyOntology,
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })
  
  output$view8 <- DT::renderDataTable({
    DT::datatable(
      data = MyDomOnt,
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })
  
  
  ################################# domain tabset #######################################
  # Create data table
  output$view <- DT::renderDataTable({
    req(input$domname)
    req(input$species_2)
    
    dom_select <-
      MyTable %>% filter(domname %in% input$domname) %>% filter(mito %in% input$mito_2) %>% filter(species %in% input$species_2) %>% filter(tp %in% input$tp_2) %>% filter(mf_pred %in% input$mf_pred_2) %>% filter(GFP %in% input$GFP_2) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_2) %>% filter(MassSpecStudies %in% input$MassSpecStudies_2) %>% filter(tissue %in% input$tissue_2)  
    
    DT::datatable(
      data = dom_select,
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })
  
  df_subset_2 <- reactive({
    a <-
      MyTable[col1 %in% input$species_2, ] %>% filter(domname %in% input$domname) %>% filter(species %in% input$species_2) %>% filter(mito %in% input$mito_2) %>% filter(tp %in% input$tp_2)  %>% filter(mf_pred %in% input$mf_pred_2) %>% filter(GFP %in% input$GFP_2) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_2)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_2) %>% filter(tissue %in% input$tissue_2) %>% select(protein)
      MyPre %>% filter(protein %in% a$protein) %>% select(protein, seq)
  })
  
  df_subset_3 <- reactive({
    MyTable[col1 %in% input$species_2, ] %>% filter(domname %in% input$domname) %>% filter(species %in% input$species_2) %>% filter(mito %in% input$mito_2)  %>% filter(tp %in% input$tp_2) %>% filter(mf_pred %in% input$mf_pred_2) %>% filter(GFP %in% input$GFP_2) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_2)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_2) %>% filter(tissue %in% input$tissue_2)  
  })
  
  output$plot = renderPlot({
    p <- ggplot(df_subset_3(), aes(factor(species)))
    p + geom_bar(fill = "#2bd47b") + stat_count(
      aes(label = ..count..),
      vjust = -1,
      geom = "text",
      position = "identity"
    ) + theme_minimal() + theme(
      axis.text.x = element_text(size = rel(1.5)),
      axis.text.y = element_text(size = rel(1.5)),
      axis.title.x = element_text(size = rel(1.5)),
      axis.title.y = element_text(
        size = rel(1.5),
        margin = margin(
          t = 0,
          r = 20,
          b = 0,
          l = 0
        )
      )
    ) + xlab("Species") + ylab("Number of proteins with domain")
  })
  
  output$plot5 = renderPlot({
    p <- ggplot(df_subset_3(), aes(factor(comb)))
    p + geom_bar(fill = "#2bd47b") + stat_count(
      aes(label = ..count..),
      vjust = -1,
      geom = "text",
      position = "identity"
    ) + theme_minimal() + theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = rel(1.5)),
      axis.title.x = element_text(size = rel(1.5)),
      axis.title.y = element_text(
        size = rel(1.5),
        margin = margin(
          t = 0,
          r = 20,
          b = 0,
          l = 0
        )
      )
    ) + xlab("Domain combinations")  + ylab("Number of proteins with combination")
  })
  
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste("domain.fasta")
    },
    
    content = function(file) {
      dataframe2fas(data.frame(df_subset_2()), file)
    }
  )
  
  output$download_data_6 <- downloadHandler(
    filename = function() {
      paste("domain.csv")
    },
    content = function(file) {
      write.csv(df_subset_3(), file, quote = FALSE)
    }
  )

  ################################# gene ontology tabset #######################################
  # Create data table
  output$view5 <- DT::renderDataTable({
    go_select <-
      MyGO %>% filter(goid %in% input$goid)  %>% filter(mito %in% input$mito_5) %>% filter(species %in% input$species_5) %>% filter(tp %in% input$tp_5) %>% filter(mf_pred %in% input$mf_pred_5) %>% filter(GFP %in% input$GFP_5) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_5)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_5) 
    
    DT::datatable(
      data = go_select,
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })
  
  df_subset_12 <- reactive({
    a <-
      MyGO[col15 %in% input$species_5, ] %>% filter(goid %in% input$goid) %>% filter(species %in% input$species_5)  %>% filter(mito %in% input$mito_5)  %>% filter(tp %in% input$tp_5)  %>% filter(mf_pred %in% input$mf_pred_5) %>% filter(GFP %in% input$GFP_5) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_5)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_5) %>% select(protein)
      MyPre %>% filter(protein %in% a$protein) %>% select(protein, seq)
  })
  
  df_subset_13 <- reactive({
    MyGO[col15 %in% input$species_5, ] %>% filter(goid %in% input$goid) %>% filter(species %in% input$species_5)  %>% filter(mito %in% input$mito_5) %>% filter(tp %in% input$tp_5) %>% filter(mf_pred %in% input$mf_pred_5) %>% filter(GFP %in% input$GFP_5) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_5)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_5)
  })
  
  output$plot6 = renderPlot({
    p <- ggplot(df_subset_13(), aes(factor(species)))
    p + geom_bar(fill = "#FFA07A") + stat_count(
      aes(label = ..count..),
      vjust = -1,
      geom = "text",
      position = "identity"
    ) + theme_minimal() + theme(
      axis.text.x = element_text(size = rel(1.5)),
      axis.text.y = element_text(size = rel(1.5)),
      axis.title.x = element_text(size = rel(1.5)),
      axis.title.y = element_text(
        size = rel(1.5),
        margin = margin(
          t = 0,
          r = 20,
          b = 0,
          l = 0
        )
      )
    ) +
      xlab("Species") + ylab("Number of proteins")
  })
  
  output$download_data_8 <- downloadHandler(
    filename = function() {
      paste("go.fasta")
    },
    
    content = function(file) {
      dataframe2fas(data.frame(df_subset_12()), file)
    }
  )
  
  output$download_data_7 <- downloadHandler(
    filename = function() {
      paste("go.csv")
    },
    content = function(file) {
      write.csv(df_subset_13(), file, quote = FALSE)
    }
  )
  
  ################################# presequence tabset #####################################
  output$view3 <- DT::renderDataTable({
    req(input$species_3)
    
    pre_select <-
      MyPre %>% filter(species %in% input$species_3) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf) %>% filter(OG_present %in% input$OG_present)  %>% filter(GFP %in% input$GFP) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas)  %>% filter(MassSpecStudies %in% input$MassSpecStudies) %>% filter(tissue %in% input$tissue)  
    
    DT::datatable(
      data = pre_select[, -11],
      options = list(pageLength = 10),
      rownames = FALSE
    )
    
  })
  
  df_subset_5 <- reactive({
    MyPre %>% filter(species %in% input$species_3) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf) %>% filter(OG_present %in% input$OG_present)  %>% filter(GFP %in% input$GFP) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas)  %>% filter(MassSpecStudies %in% input$MassSpecStudies) %>% filter(tissue %in% input$tissue) %>% select(protein, seq)
  })
  
  df_subset_6 <- reactive({
    MyPre %>% filter(species %in% input$species_3) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf) %>% filter(OG_present %in% input$OG_present)  %>% filter(GFP %in% input$GFP) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas)  %>% filter(MassSpecStudies %in% input$MassSpecStudies) %>% filter(tissue %in% input$tissue)
  })
  
  df_subset_7 <- reactive({
    MyPre %>% filter(species %in% input$species_3) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf) %>% filter(OG_present %in% input$OG_present) %>% filter(GFP %in% input$GFP) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas)  %>% filter(MassSpecStudies %in% input$MassSpecStudies) %>% filter(tissue %in% input$tissue)
  })
  
  df_subset_9 <- reactive({
    MyPre %>% filter(species %in% input$species_3) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf) %>% filter(OG_present %in% input$OG_present) %>% select(mppsite) %>% filter(GFP %in% input$GFP) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas)  %>% filter(MassSpecStudies %in% input$MassSpecStudies) %>% filter(tissue %in% input$tissue)
  })
  
  output$plot3 = renderPlot({
    p <- ggplot(df_subset_6(), aes(factor(species)))
    p + geom_bar(fill = "#adbce6") + stat_count(
      aes(label = ..count..),
      vjust = -1,
      geom = "text",
      position = "identity"
    ) + theme_minimal() + theme(
      axis.text.x = element_text(size = rel(1.5)),
      axis.text.y = element_text(size = rel(1.5)),
      axis.title.x = element_text(size = rel(1.5)),
      axis.title.y = element_text(
        size = rel(1.5),
        margin = margin(
          t = 0,
          r = 20,
          b = 0,
          l = 0
        )
      )
    ) +
      xlab("Species") + ylab("Number of proteins")
  })
  
  output$download_data_3 <- downloadHandler(
    filename = function() {
      paste("mts.fasta")
    },
    content = function(file) {
      dataframe2fas(data.frame(df_subset_5()), file)
    }
  )
  
  output$download_data_4 <- downloadHandler(
    filename = function() {
      paste("mts.csv")
    },
    content = function(file) {
      write.csv(df_subset_7(), file, quote = FALSE)
    }
  )
  
  
  ################################# orthology tabset #######################################
  output$view2 <- DT::renderDataTable({
    req(input$og)
    req(input$species)
    
    ortho_select <-
      MyOrtho %>% filter(og %in% input$og) %>% filter(species %in% input$species) %>% filter(mito %in% input$mitolabel) %>% filter(tp %in% input$tp_3) %>% filter(mf_pred %in% input$mf_pred_3) %>% filter(GFP %in% input$GFP_3) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_3)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_3) %>% filter(tissue %in% input$tissue_3)
    
    DT::datatable(
      data = ortho_select,
      options = list(pageLength = 10),
      rownames = FALSE
    )
  })
  
  df_subset <- reactive({
    s <-
      MyOrtho[col4 %in% input$species, ] %>% filter(og %in% input$og) %>% filter(species %in% input$species) %>% filter(mito %in% input$mitolabel) %>% filter(tp %in% input$tp_3) %>% filter(mf_pred %in% input$mf_pred_3) %>% filter(GFP %in% input$GFP_3) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_3)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_3) %>% filter(tissue %in% input$tissue_3) %>% select(protein)
      MyPre %>% filter(protein %in% s$protein) %>% select(protein, seq)
  })
  df_subset_4 <- reactive({
    MyOrtho[col4 %in% input$species, ] %>% filter(og %in% input$og) %>% filter(species %in% input$species) %>% filter(mito %in% input$mitolabel) %>% filter(tp %in% input$tp_3) %>% filter(mf_pred %in% input$mf_pred_3) %>% filter(GFP %in% input$GFP_3) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_3)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_3) %>% filter(tissue %in% input$tissue_3)
  })
  
  output$plot2 = renderPlot({
    p <- ggplot(df_subset_4(), aes(factor(species)))
    p + geom_bar(fill = "#e6add8") + stat_count(
      aes(label = ..count..),
      vjust = -1,
      geom = "text",
      position = "identity"
    ) + theme_minimal() + theme(
      axis.text.x = element_text(size = rel(1.5)),
      axis.text.y = element_text(size = rel(1.5)),
      axis.title.x = element_text(size = rel(1.5)),
      axis.title.y = element_text(
        size = rel(1.5),
        margin = margin(
          t = 0,
          r = 20,
          b = 0,
          l = 0
        )
      )
    ) +
      xlab("Species") + ylab("Number of proteins")
  })
  
  output$download_data2 <- downloadHandler(
    filename = function() {
      paste("ortho.fasta")
    },
    content = function(file) {
      dataframe2fas(data.frame(df_subset()), file)
    }
  )
}

# part 3 Bind ui and server together
shinyApp(ui = ui, server = server)
