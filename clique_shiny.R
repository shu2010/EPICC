
.libPaths(c( "/home/shu/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths() ) )
library(shiny)
library(gridExtra)
library(ggnet)
library(ggplot2)
#library(sna)
#library(intergraph)
#library(network)

ui <- fluidPage(
  
  # App title ----
  titlePanel("EPICC cliques"),
  
  selectInput("select", label = "Select Clique", 
              choices = c("select clique", c("TCGA_new_OV_2951_C2_60", "TCGA_new_SKCM_2002_C1_109", "TCGA_new_SKCM_1720_C1_72", 
                                               "TCGA_new_SKCM_2002_C1_102", "TCGA_new_GBM_2476_C1_15", "TCGA_new_GBM_2476_C1_31", 
                                               "TCGA_new_GBM_2476_C1_9", "TCGA_new_LUAD_1047_C1_36", "TCGA_new_SKCM_2002_C1_98", 
                                               "TCGA_new_GBM_2476_C1_16", "TCGA_new_LUAD_1047_C1_41", "TCGA_new_GBM_2476_C1_21", 
                                               "TCGA_new_GBM_2476_C1_20", "TCGA_new_GBM_2476_C1_29", "TCGA_new_SKCM_2002_C1_92", 
                                               "TCGA_new_GBM_2476_C1_13", "TCGA_new_GBM_2476_C1_12", "TCGA_new_BRCA_650_C2_6", 
                                               "TCGA_new_SKCM_2002_C1_104", "TCGA_new_SKCM_2002_C1_100", "TCGA_new_SKCM_2002_C2_89", 
                                               "TCGA_new_GBM_2476_C1_18", "TCGA_new_BRCA_650_C2_9", "TCGA_new_SKCM_2002_C1_108", 
                                               "TCGA_new_SKCM_1720_C1_68", "TCGA_new_LUAD_1047_C1_42", "TCGA_new_SKCM_1736_C2_70", 
                                               "TCGA_new_GBM_2476_C1_26")), 
              selected = NULL, multiple = FALSE),
  selectInput("Anno", label = "Annotation", choices = c("select Annotation", "MsigDB", "GO_MF", "GO_BP", "GO_CC"), 
              selected = NULL, multiple = FALSE),
  selectInput("Can", label = "Cancer", choices = c("Cancer type", c("BRCA", "GBM", "LUAD", "OV", "PAAD", "PRAD", "SARC", "SKCM", "All")), 
              selected = NULL, multiple = FALSE),
  sliderInput("in1", "Mutation frequency per patient", min = 1, max = 25, 
              value = 1, step= 1),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Plot",
               fluidRow(
                 plotOutput("distPlot"),
                 tableOutput("view")))),
    tabPanel("Summary",  tableOutput("view1"))
  )
)


##Clique data
filt_cliques <- readRDS("filt_cliques.RDS")

##Clique graphs
cliq_vis <- readRDS("cliq_vis.RDS")

##sample x mutation matrix
samp_mut_mat <- readRDS("samp_mut_mat.rds")

##Pathways
p_hit_msig <- readRDS("p_hit_msig.RDS")
##GO annotations
p_hit_GO_MF <- readRDS("p_hit_GO_MF.RDS")
p_hit_GO_BP <- readRDS("p_hit_GO_BP.RDS")
p_hit_GO_CC <- readRDS("p_hit_GO_CC.RDS")

##survival plots data
p_surv_plots <- readRDS("p_surv_plots.RDS")


# Define server logic ----
server <- function(input, output, session) { 
 

##update number of proteins based on clique inpute
 observeEvent(input$select,  {
   updateSliderInput(session = session, inputId = "in1",
                     max = length(unlist(strsplit(as.character(filt_cliques[filt_cliques$cliq_id == input$select,]$mods), ","))))
 })
  
  output$distPlot <- renderPlot({
    validate(
      need(input$select != "select clique", "Please select a clique!!"),
      need(input$Anno != "select Annotation", "Please select annotation database!!"),
      need(input$Can != "Cancer type", "Please select Cancer type!!"),
      need(length(unlist(strsplit(as.character(filt_cliques[filt_cliques$cliq_id == input$select,]$mods), ","))) != 0, 
           "Please select a new clique")
    )

##cancer histotype frequency
   
    filt_cliques_mod <- unlist(strsplit(as.character(filt_cliques[filt_cliques$cliq_id == input$select,]$mods), ","))
    filt_cliques_mod2 <- samp_mut_mat[,colnames(samp_mut_mat) %in% filt_cliques_mod]
    can_freq <- rownames(filt_cliques_mod2[rowSums(filt_cliques_mod2) >= input$in1,])
    
    if(length(can_freq) == 0){
      plot(1,1,col="white", axes=FALSE)
      text(1,1,"Please choose a lower mutation frequency", cex = 2)
    }else{
      can_freq1 <- as.data.frame(table(unlist(lapply(strsplit(can_freq, split = "_"),function(y)y[3]))))
      x <- can_freq1
          ##p1 condition starts here
      colnames(x) <- c("Cancer", "Patients")
      
      #   plot(x, xlab = "Cancer", main = "Cancer prevalence")
      p <- ggplot(data=x, aes(x=Cancer, y=Patients)) + geom_bar(stat="identity", fill="gray") + 
        ggtitle(paste(input$in1, "Mutation(s) in", sum(x$Patients), "patients", sep = " ")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text=element_text(size=8,face="bold"), plot.title = element_text(size = 10, face = "bold")) 
      p1 <- p + coord_flip() 
      
##protein mutation frequency
      filt_cliques_mod2 <- samp_mut_mat[,colnames(samp_mut_mat) %in% filt_cliques_mod]
      mp1_prot <- colSums(filt_cliques_mod2[rowSums(filt_cliques_mod2) >= input$in1,,drop=FALSE])
      d_mp1_prot <- as.data.frame(mp1_prot)
      colnames(d_mp1_prot) <- "Patients"
      d_mp1_prot$prot <- rownames(d_mp1_prot)
      p2 <- ggplot(data=d_mp1_prot, aes(x=prot, y=Patients)) + geom_bar(stat="identity", fill="gray") + 
        ggtitle(paste("Patients:", sum(x$Patients), ";mut/pat:",input$in1, sep = "")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text=element_text(size=8,face="bold"), plot.title = element_text(size = 10, face = "bold")) 
      p2 <- p2 + coord_flip()
      
##clique graphs
      p4 <- cliq_vis[[input$select]]
###   
      ##survival plots
      #Null plots
      if(is.null(p_surv_plots[[input$select]])){
        df <- data.frame()
        p5 <- ggplot(df) + geom_point() + xlim(0, 1) + ylim(0, 1) + ggtitle("No survival plot available for this module!!")
      }
      else if(is.null(p_surv_plots[[input$select]][[input$Can]])){
        df <- data.frame()
        p5 <- ggplot(df) + geom_point() + xlim(0, 1) + ylim(0, 1) + ggtitle("No survival plot available for this cancer type!!")
      }
      else if(!is.null(p_surv_plots[[input$select]][[input$Can]])){
        p5 <- p_surv_plots[[input$select]][[input$Can]]
      }
      grid.arrange(p1, p2, p4, p5, ncol=4, widths = c(8,8,10,10))
    }
  })
  ##table
  d_tab <- reactive({
    if(input$Anno == "MsigDB"){
      df_tab <- p_hit_msig[[input$select]]
      df_tab[,c(2,3)] <- sapply(df_tab[,c(2,3)], as.integer)
      df_tab
    }
    else if(input$Anno == "GO_MF"){
      df_tab <- p_hit_GO_MF[[input$select]]
      df_tab
    }
    else if(input$Anno == "GO_BP"){
      df_tab <- p_hit_GO_BP[[input$select]]
      df_tab
    }
    else if(input$Anno == "GO_CC"){
      df_tab <- p_hit_GO_CC[[input$select]]
      df_tab
    }
  })
  output$view <- renderTable({
    d_tab()
  }, digits = 5)
  
}

# Run the app ----
shinyApp(ui = ui, server = server)