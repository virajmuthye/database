library(shiny)
library(readr)
library(tidyverse)
library(seqRFLP)
library(ggplot2)
library(Rcpp)
library(dplyr)
library(DT)
library(shinythemes)
library(data.table)

#read in the datasets
MyTable <- fread("pfam_mitominer_tissue.csv")
MyOrtho <- fread("orthology_mitominer_domcomb_tissue.csv")
MyPre <- fread("presequence_mitominer_domcomb_tissue.csv")
MyGO <- fread("go_mitominer.csv")
MyGene <- fread("gene.csv", fill=TRUE)
MyOntology <- fread("ontology.csv")

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
          "The Metazoan Mitochondrial Proteome Database (MMPdb) is a database for facilitating comparative analysis of experimentally-characterized 
           mitochondrial proteomes of animals - Homo sapiens (hsap), Mus musculus (mmus), Caenorhabditis elegans (cele), Drosophila melanogaster (dmel) 
           and two outgroups - Acanthamoeba castellanii (acas) and Saccharomyces cerevisiae (scer). Each species is denoted by a four letter abbreviation 
          listed in the parentheses."
        ),
        br(),
        p(
          "MMPdb is organized into four tabsets: Orthology, MTS, Domain and Gene Ontology"
        ),
        br(),
        strong("Orthology Tabset"),
        p(
          "The Orthology tabset can be used to analyze and download proteins belonging to a specific Orthology Group (OG) of interest.
          An OG is a group of orthologous proteins in the species mentioned above, as identified by Proteinortho v5.16b."
        ),
        br(),
        strong("Mitochondrial Targeting Signal (MTS) Tabset"),
        p(
          "The majority of mitochondrial proteins are imported into the organelle via the N-terminus Mitochondrial Targeting Signal (MTS)
          pathway. MTS were identified using TargetP and MitoFates. This tabset allows for the exploration, analysis and download of the 
          MTS-prediction results for all the proteins from six species listed above."
        ),
        br(),
        strong("Protein Domain Tabset"),
        p(
          "The Domain tabset can be used to analyze and download proteins containing a protein domain of interest. Protein domains were
          identified using PFAM"
        ),
        br(),
        strong("Gene Ontology Tabset"),
        p(
          "Gene Ontology analysis was performed using Pannzer2.0 [10].This tabset can be used to analyze and download proteins mapped to
          a GO term of interest."
        ),
        br(),
        p(
          "Details regarding construction of the database are given in our manuscript titled 'MMPdb and MitoPredictor- tools for facilitating
          comparative analysis of animal mitochondrial proteomes (under review)'"
        ),
        br(),
        strong(
          "It may take a few minutes for the database to load the database."
        ),
        br(),
        p(
          "For queries or suggestions, contact Viraj Muthye at vrmuthye@iastate.edu"
         )
        ),
      
      ######################################################################################################################################333
      
      tabPanel(
        "Orthology",
        sidebarPanel(
          helpText("All field are required fields"),
          checkboxGroupInput(
            "species",
            label = "Select species",
            choices = sort(unique(col4)),
            selected = sort(unique(col4))
          ),
          helpText(
            "Each species is denoted by a four letter abbreviation. Check the About tab for details about the abbreviations."
          ),
          checkboxGroupInput(
            "mitolabel",
            label = "Known subcellular localization (mitochondrial/non-mitochondrial):",
            c("Mitochondrial" = "Y", "Non-mitochondrial" = "N"),
            selected = c("Y","N")
          ),
          checkboxGroupInput(
            "tp_3",
            label = "TargetP prediction",
            choices = sort(unique(col20)),
            selected = sort(unique(col20))
          ),
          helpText(
            "M:Mitochondrial, S:Secreted, _:No prediction"
            ),
          checkboxGroupInput(
            "mf_pred_3",
            label = "MitoFates prediction",
            choices = sort(unique(col21)),
            selected = sort(unique(col21))
          ),
          helpText(
            "M:MTS detected, N:No MTS detected"
          ),
          checkboxGroupInput(
            "GFP_3",
            label = "GFP evidence of mitochondrial localization (Mammals only)",
            choices = sort(unique(col24)),
            selected = sort(unique(col24))
          ),
          helpText(
            "Y: GFP evidence of mitochondrial localization, N: No GFP evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "HumanProteinAtlas_3",
            label = "Human Protein Atlas (HPA) evidence for mitochondrial localization (Human only)",
            choices = sort(unique(col23)),
            selected = sort(unique(col23))
          ),
          helpText(
            "Y: HPA evidence of mitochondrial localization, N: No HPA evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "MassSpecStudies_3",
            label = "Mass Spectrometry evidence for mitochondrial localization (Mammals only)",
            choices = sort(unique(col25)),
            selected = sort(unique(col25))
          ),
          
          helpText(
            "Y: MS evidence of mitochondrial localization, N: No MS evidence, nd: No Data"
          ),
          downloadButton(outputId = "download_data2", label = "Download sequences"),
          helpText("Download selected sequences is a FASTA format"),
          br(),
          strong("Key for table"),
          br(),
          strong("OG"),
          p("Unique Orthologous Group (OG) identifier."),
          strong("Mitochondrial"),
          p("Mitochondrial/Non-mitochondrial protein"),
          strong("MPP cleavage site"),
          p("Mitochondrial Processing Peptidase cleavage site"),
          strong("GO"),
          p("Mitochondrial localization evidence from MitoMiner (Gene Ontology)"),
          strong("Human Protein Atlas"),
          p("Mitochondrial localization evidence from MitoMiner (Human Protein Atlas)"),
          strong("GFP"),
          p("Mitochondrial localization evidence from MitoMiner (GFP analysis)"),
          strong("Mass Spectrometry"),
          p("Mitochondrial localization evidence from MitoMiner (Mass spectrometry studies)"),
          strong("Domain combination"),
          p("List of all protein domains"),
          br()
        ),
        
        
        mainPanel(
          br(),
          strong(
            "Search for Orthologous Group (OG) number"
          ),
          p(
            "The Orthologous Group (OG) number is the primary search query for this tabset. To fetch the OG number of a desired protein,
            enter the Uniprot protein accession (ID) / protein identification code (protein name) for human, mouse, and D. melanogaster proteins.
            For example, the Uniprot protein accession (ID) for fumarate hydratase in humans is 'P07954' and the protein name is 'FUMH'."
          ),
          p("
            For C.elegans, enter the WormBase protein ID or the Wormbase protein name. For the C. elegans protein fumarate hydratase isoform a,
            enter either 'CE11580' or 'fum-1 isoform a'."
          ),
          br(),
          DT::dataTableOutput(outputId = "view6"),
          br(),
          p(
            "Enter the OG number in the search bar below to analyze proteins belonging to that OG:"
          ),
          br(),
          selectInput(
            inputId = "og",
            label = "Enter OG number",
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
          helpText("All field are required fields"),
          checkboxGroupInput(
            "species_3",
            label = "Select species",
            choices = sort(unique(col10)),
            selected = sort(unique(col10))
          ),
          helpText(
            "Each species is denoted by a four letter abbreviation. Check the About tab for details about the abbreviations."
          ),
          checkboxGroupInput(
            "mito",
            label = "Known subcellular localization (mitochondrial/non-mitochondrial)",
            c("Mitochondrial" = "Y", "Non-mitochondrial" = "N"),
            selected = c("Y")
          ),
         checkboxGroupInput(
            "OG_present",
            label = "Present/absent in Orthologous Group (OG)",
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
            "M:Mitochondrial, S:Secreted, _:No prediction"
          ),
          checkboxGroupInput(
            "rc",
            label = "Reliability Class (RC)",
            choices = sort(unique(col8)),
            selected = sort(unique(col8))
          ),
          helpText(
            "RC ranges from 1 (strongest prediction) to 5 (weakest prediction)"
          ),
          checkboxGroupInput(
            "mf",
            label = "MitoFates prediction",
            choices = sort(unique(col9)),
            selected = sort(unique(col9))
          ),
          helpText(
            "N: No MTS detected, Y: MTS detected"
          ),
          checkboxGroupInput(
            "GFP",
            label = "GFP evidence of mitochondrial localization (Mammals only)",
            choices = sort(unique(col28)),
            selected = sort(unique(col28))
          ),
          helpText(
            "Y: GFP evidence of mitochondrial localization, N: No GFP evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "HumanProteinAtlas",
            label = "Human Protein Atlas (HPA) evidence of mitochondrial localization (Human only)",
            choices = sort(unique(col27)),
            selected = sort(unique(col27))
          ),
          helpText(
            "Y: HPA evidence of mitochondrial localization, N: No evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "MassSpecStudies",
            label = "Mass Spectrometry evidence of mitochondrial localization (Mammals only)",
            choices = sort(unique(col29)),
            selected = sort(unique(col29))
          ),
          helpText(
            "Y: MS evidence of mitochondrial localization, N: No MS evidence, nd: No Data"
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
          strong("OG number"),
          p("Unique Orthologous Group (OG) identifier."),
          strong("Mitochondrial"),
          p("Mitochondrial/Non-mitochondrial protein"),
          strong("MPP cleavage site"),
          p("Mitochondrial Processing Peptidase cleavage site"),
          strong("GO"),
          p("Mitochondrial localization evidence from MitoMiner (Gene Ontology)"),
          strong("Human Protein Atlas"),
          p("Mitochondrial localization evidence from MitoMiner (Human Protein Atlas)"),
          strong("GFP"),
          p("Mitochondrial localization evidence from MitoMiner (GFP analysis)"),
          strong("Mass Spectromtery"),
          p("Mitochondrial localization evidence from MitoMiner (Mass spectrometry studies)"),
          strong("Domain combination"),
          p("List of all protein domains"),
          br()
        ),
        mainPanel(
          br(),
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
          helpText("All field are required fields"),
          checkboxGroupInput(
            "species_2",
            label = "Select species",
            choices = sort(unique(col1)),
            selected = sort(unique(col1))
          ),
          helpText(
            "Each species is denoted by a four letter abbreviation. Check the About tab for details about the abbreviations."
          ),
          checkboxGroupInput(
            "mito_2",
            label = "Known subcellular localization (mitochondrial/non-mitochondrial)",
            c("Mitochondrial" = "Y", "Non-mitochondrial" = "N"),
            selected = c("Y", "N")
          ),
          checkboxGroupInput(
            "tp_2",
            label = "TargetP prediction",
            choices = sort(unique(col11)),
            selected = sort(unique(col11))
          ),
          helpText(
            "M:Mitochondrial,S:Secreted,_:No prediction"
          ),
          checkboxGroupInput(
            "mf_pred_2",
            label = "MitoFates prediction",
            choices = sort(unique(col12)),
            selected = sort(unique(col12))
          ),
          helpText(
            "N: No MTS detected, Y: MTS detected"
          ),
          
          checkboxGroupInput(
            "GFP_2",
            label = "GFP evidence for mitochondrial localization (Mammals only)",
            choices = sort(unique(col31)),
            selected = sort(unique(col31))
          ),
          helpText(
            "Y: GFP evidence of mitochondrial localization, N: No GFP evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "HumanProteinAtlas_2",
            label = "Human Protein Atlas (HPA) evidence for mitochondrial localization (Human only)",
            choices = sort(unique(col30)),
            selected = sort(unique(col30))
          ),
          helpText(
            "Y: HPA evidence of mitochondrial localization, N: No HPA evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "MassSpecStudies_2",
            label = "Mass Spectrometry evidence for mitochondrial localization (Mammals only)",
            choices = sort(unique(col32)),
            selected = sort(unique(col32))
          ),
          helpText(
            "Y: MS evidence of mitochondrial localization, N: No MS evidence, nd: No Data"
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
          strong("OG number"),
          p("Unique Orthologous Group (OG) identifier."),
          strong("Mitochondrial"),
          p("Mitochondrial/Non-mitochondrial protein"),
          strong("Domain combination"),
          p("List of all protein domains"),
          strong("MPP cleavage site"),
          p("Mitochondrial Processing Peptidase cleavage site"),
          strong("GO"),
          p("Mitochondrial localization evidence from MitoMiner (Gene Ontology)"),
          strong("Human Protein Atlas"),
          p("Mitochondrial localization evidence from MitoMiner (Human Protein Atlas)"),
          strong("GFP"),
          p("Mitochondrial localization evidence from MitoMiner (GFP analysis)"),
          strong("Mass Spectrometry"),
          p("Mitochondrial localization evidence from MitoMiner (Mass spectrometry studies)"),
          br()
        ),
        mainPanel(
          br(),
          selectInput(
            inputId = "domname",
            label = "Enter protein domain name",
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
          helpText("All field are required fields"),
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
            label = "Known subcellular localization (mitochondrial/non-mitochondrial)",
            c("Mitochondrial" = "Y", "Non-mitochondrial" = "N"),
            selected = c("Y", "N")
          ),
          checkboxGroupInput(
            "tp_5",
            label = "TargetP prediction",
            choices = sort(unique(col18)),
            selected = sort(unique(col18))
          ),
          helpText(
            "M:Mitochondrial,S:Secreted,_:No prediction"
          ),
          checkboxGroupInput(
            "mf_pred_5",
            label = "MitoFates prediction",
            choices = sort(unique(col19)),
            selected = sort(unique(col19))
          ),
          helpText(
            "N: No MTS detected, Y: MTS detected"
          ),
          checkboxGroupInput(
            "GFP_5",
            label = "GFP evidence (Mammals only)",
            choices = sort(unique(col37)),
            selected = sort(unique(col37))
          ),
          helpText(
            "Y: GFP evidence of mitochondrial localization, N: No GFP evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "HumanProteinAtlas_5",
            label = "Human Protein Atlas (HPA) evidence of mitochondrial localization(Human only)",
            choices = sort(unique(col36)),
            selected = sort(unique(col36))
          ),
          helpText(
            "Y: HPA evidence of mitochondrial localization, N: No HPA evidence, nd: No Data"
          ),
          checkboxGroupInput(
            "MassSpecStudies_5",
            label = "Mass Spectrometry evidence of mitochondrial localization(Mammals only)",
            choices = sort(unique(col38)),
            selected = sort(unique(col38))
          ),
          helpText(
            "Y: MS evidence of mitochondrial localization, N: No MS evidence, nd: No Data"
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
          strong("OG number"),
          p("Unique Orthologous Group (OG) identifier."),
          strong("Mitochondrial"),
          p("Mitochondrial/Non-mitochondrial protein"),
          strong("GO ID"),
          p("Gene Ontology term ID"),
          strong("MPP cleavage site"),
          p("Mitochondrial Processing Peptidase cleavage site"),
          strong("GO"),
          p("Mitochondrial localization evidence from MitoMiner (Gene Ontology)"),
          strong("Human Protein Atlas"),
          p("Mitochondrial localization evidence from MitoMiner (Human Protein Atlas)"),
          strong("GFP"),
          p("Mitochondrial localization evidence from MitoMiner (GFP analysis)"),
          strong("Mass Spectromtery"),
          p("Mitochondrial localization evidence from MitoMiner (Mass spectrometry studies)"),
          br()
        ),
        mainPanel(
          br(),
          strong(
            "Search for the GO ID"
            ),
          p(
            "The primary search query for this tabset is the GO ID number. To search for the GO ID, start typing the 
            Go term (mitochondrion) or keywords of interest (mitochondrial inner membrane) in the search bar below"
          ),
          br(),
          DT::dataTableOutput(outputId = "view61"),
          br(),
          selectInput(
            inputId = "goid",
            label = "Enter Gene Ontology (GO) ID",
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
      options = list(pageLength = 1),
      rownames = FALSE
    )
  })
  
  output$view61 <- DT::renderDataTable({
    DT::datatable(
      data = MyOntology,
      options = list(pageLength = 1),
      rownames = FALSE
    )
  })

  ################################# domain tabset #######################################
  # Create data table
  output$view <- DT::renderDataTable({
    req(input$domname)
    req(input$species_2)
    
    dom_select <-
      MyTable %>% filter(domname %in% input$domname) %>% filter(mito %in% input$mito_2) %>% filter(species %in% input$species_2) %>% filter(tp %in% input$tp_2) %>% filter(mf_pred %in% input$mf_pred_2) %>% filter(GFP %in% input$GFP_2) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_2) %>% filter(MassSpecStudies %in% input$MassSpecStudies_2)
    
    DT::datatable(
      data = dom_select[, -c(1,10,15)],
      options = list(pageLength = 50),
      colnames = c('Renamed proteinID','UniprotID','ProteinID','PFAM ID','Domain name','Domain combination','OG number','Mitochondrial','TargetP prediction','TargetP RC','MitoFates probability','MitoFates prediction','GO','Human Protein Atlas','GFP','Mass Spectromtery','Tissue')
    )
  })
  
  df_subset_2 <- reactive({
    a <-
      MyTable[col1 %in% input$species_2, ] %>% filter(domname %in% input$domname) %>% filter(species %in% input$species_2) %>% filter(mito %in% input$mito_2) %>% filter(tp %in% input$tp_2)  %>% filter(mf_pred %in% input$mf_pred_2) %>% filter(GFP %in% input$GFP_2) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_2)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_2) %>% select(protein)
      MyPre %>% filter(protein %in% a$protein) %>% select(protein, seq)
  })
  
  df_subset_3 <- reactive({
    MyTable[col1 %in% input$species_2, ] %>% filter(domname %in% input$domname) %>% filter(species %in% input$species_2) %>% filter(mito %in% input$mito_2)  %>% filter(tp %in% input$tp_2) %>% filter(mf_pred %in% input$mf_pred_2) %>% filter(GFP %in% input$GFP_2) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_2)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_2) 
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
      data = go_select[, -c(1,7,12)],
      options = list(pageLength = 50),
      colnames = c('Renamed proteinID','UniprotID','ProteinID','GO ID','Mitochondrial','TargetP prediction','TargetP RC','MitoFates probability','MitoFates prediction','GO','Human Protein Atlas','GFP','Mass Spectromtery')
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
      MyPre %>% filter(species %in% input$species_3) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf) %>% filter(OG_present %in% input$OG_present)  %>% filter(GFP %in% input$GFP) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas)  %>% filter(MassSpecStudies %in% input$MassSpecStudies)
    
    DT::datatable(
      data = pre_select[, -c(1,5,7,12,13,14)],
      options = list(pageLength = 50),
      colnames= c('Renamed proteinID','UniprotID','ProteinID','Mitochondrial','TargetP prediction','TargetP RC','MitoFates probability','MitoFates prediction','GO','Human Protein Atlas','GFP','Mass Spectromtery','Domain combination','Tissue')
    )
    
  })
  
  df_subset_5 <- reactive({
    MyPre %>% filter(species %in% input$species_3) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf) %>% filter(OG_present %in% input$OG_present)  %>% filter(GFP %in% input$GFP) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas)  %>% filter(MassSpecStudies %in% input$MassSpecStudies) %>% select(protein, seq)
  })
  
  df_subset_6 <- reactive({
    MyPre %>% filter(species %in% input$species_3) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf) %>% filter(OG_present %in% input$OG_present)  %>% filter(GFP %in% input$GFP) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas)  %>% filter(MassSpecStudies %in% input$MassSpecStudies)
  })
  
  df_subset_7 <- reactive({
    MyPre %>% filter(species %in% input$species_3) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf) %>% filter(OG_present %in% input$OG_present) %>% filter(GFP %in% input$GFP) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas)  %>% filter(MassSpecStudies %in% input$MassSpecStudies)
  })
  
  df_subset_9 <- reactive({
    MyPre %>% filter(species %in% input$species_3) %>% filter(mito %in% input$mito) %>% filter(tp %in% input$tp) %>% filter(rc %in% input$rc) %>% filter(mf %in% input$mf) %>% filter(OG_present %in% input$OG_present) %>% select(mppsite) %>% filter(GFP %in% input$GFP) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas)  %>% filter(MassSpecStudies %in% input$MassSpecStudies)
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
      MyOrtho %>% filter(og %in% input$og) %>% filter(species %in% input$species) %>% filter(mito %in% input$mitolabel) %>% filter(tp %in% input$tp_3) %>% filter(mf_pred %in% input$mf_pred_3) %>% filter(GFP %in% input$GFP_3) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_3)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_3)
    
    DT::datatable(
      data = ortho_select[, -c(1,7)],
      options = list(pageLength = 50),
      #rownames = FALSE,
      colnames = c('OG','Renamed ProteinID','ProteinID','Protein name','Mitochondrial','TargetP prediction','TargetP RC','MitoFates probablity','MitoFates prediction','MPP cleavage site','GO','Human Protein Atlas','GFP','Mass Spectrometry','Domain combination','Tissue')
      )
  })
  
  df_subset <- reactive({
    s <-
      MyOrtho[col4 %in% input$species, ] %>% filter(og %in% input$og) %>% filter(species %in% input$species) %>% filter(mito %in% input$mitolabel) %>% filter(tp %in% input$tp_3) %>% filter(mf_pred %in% input$mf_pred_3) %>% filter(GFP %in% input$GFP_3) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_3)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_3) %>% select(protein)
      MyPre %>% filter(protein %in% s$protein) %>% select(protein, seq)
  })
  df_subset_4 <- reactive({
    MyOrtho[col4 %in% input$species, ] %>% filter(og %in% input$og) %>% filter(species %in% input$species) %>% filter(mito %in% input$mitolabel) %>% filter(tp %in% input$tp_3) %>% filter(mf_pred %in% input$mf_pred_3) %>% filter(GFP %in% input$GFP_3) %>% filter(HumanProteinAtlas %in% input$HumanProteinAtlas_3)  %>% filter(MassSpecStudies %in% input$MassSpecStudies_3)
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
