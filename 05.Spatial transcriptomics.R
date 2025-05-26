# 00 Package --------------------------------------------------------------

library(Seurat)
library(patchwork)
library(tidyverse)
library(viridis)

outdir <- 'Sci_Immunol2023/GSM7097913/analysis'
if(!dir.exists(outdir)) dir.create(outdir,recursive = TRUE)
# 01 Data import ----------------------------------------------------------------------

sy_6w <- Load10X_Spatial('ST/Sci_Immunol2023/GSM7097913/',slice = 'thymus_6wkt1') %>% 
  SCTransform(assay = 'Spatial') %>% 
  RunPCA(assay = 'SCT') %>% 
  FindNeighbors(assay = 'SCT') %>% 
  FindClusters(resolution = 1,assay = 'SCT') %>% 
  RunUMAP(dims = 1:30, reduction = "pca",assay = 'SCT')

# 02 Region identification ----------------------------------------------------------


# 02.1 Image based --------------------------------------------------------

library(EBImage)
library(rjson)
image_region <- function(input_dir,k,outdir){
  library(EBImage)
  library(rjson)
  scale_factors <- fromJSON(file = file.path(input_dir,"spatial/scalefactors_json.json"))
  tissue_lowres_scalef <- scale_factors$tissue_lowres_scalef
  image <- EBImage::readImage(file.path(input_dir,"spatial/tissue_lowres_image.png"))
  spot <- read.csv(file.path(input_dir,"spatial/tissue_positions.csv"), header = TRUE)
  spot$x_lowres <- spot$pxl_col_in_fullres * tissue_lowres_scalef
  spot$y_lowres <- spot$pxl_row_in_fullres * tissue_lowres_scalef
  img_height <- dim(image)[1]
  img_width <- dim(image)[2]
  spot$x_lowres <- pmax(pmin(spot$x_lowres, img_width), 1)
  spot$y_lowres <- pmax(pmin(spot$y_lowres, img_height), 1)
  plot(image)
  points(spot$x_lowres, spot$y_lowres, col = "red", pch = 16, cex = 0.5, main = "Spot Alignment Check")
  gray_img <- EBImage::channel(image, "gray")
  plot(gray_img)
  dev.print(pdf, 
    file = file.path(outdir,paste0(basename(input_dir),"_gray_image.pdf")), 
    width = dim(gray_img)[2]/72,
    height = dim(gray_img)[1]/72)
  denoised_img <- EBImage::normalize(gray_img)
  denoised_img <- gblur(denoised_img, sigma = 1)
  display(denoised_img, title = "Denoised Image")
  dev.print(pdf, 
    file = file.path(outdir,paste0(basename(input_dir),"_denoised_img.pdf")), 
    width = dim(denoised_img)[2]/72,
    height = dim(denoised_img)[1]/72
  )
  spot$gray_value <- sapply(1:nrow(spot), function(i) {
    gray_img[spot$x_lowres[i],spot$y_lowres[i]]
  })
  gray_df <- data.frame(Gray = spot$gray_value)
  set.seed(123)
  kmeans_result <- kmeans(gray_df, centers = k, nstart = 20)
  spot$cluster <- factor(kmeans_result$cluster, labels = paste0("Region", 1:k))
  p <- ggplot(spot, aes(x = x_lowres, y = -y_lowres, color = cluster)) +
    geom_point(size = 1) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Spatial Distribution of Clusters",
         x = "X Coordinate", y = "Y Coordinate") +
    theme_classic()
  print(p)
  return(spot)
}
w6_1 <- image_region('Sci_Immunol2023/GSM7097913',k=4,outdir)
p_density_w6_1 <- w6_1 %>% 
  ggplot(aes(x = gray_value))+
  geom_density(fill='gray',alpha=0.4)+labs(x='Gray value')+theme_light()
p_density_w6_1
ggsave(p_density_w6_1, filename = file.path(outdir,"Figure/Figure1.density_plot.pdf"), width = 8, height = 6)

pw6_1 <- ggplot(w6_1, aes(x = x_lowres, y = -y_lowres, color = cluster)) +
  geom_point(size = 1) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Spatial Distribution of Clusters",
       x = "X Coordinate", y = "Y Coordinate",color='') +
  theme_void()
pw6_1

ggsave(pw6_1, filename = file.path(outdir,"Figure/Figure1.region_identification.pdf")), width = 8, height = 6)

w6_1_image <- w6_1 %>% dplyr::mutate(area_thy=ifelse(cluster=='Region1','Medulla_image',"Cortex_image")) %>% 
  column_to_rownames('barcode') %>% dplyr::select(area_thy)


# 02.2 Krt5 gene based ----------------------------------------------------------------------

krt5_expression <- FetchData(sy_6w, vars = c('orig.ident',"Krt5","Aire"))
# quantile(krt5_expression$Krt5,0.7)
threshold <- quantile(krt5_expression$Krt5,0.7)
p_density_krt5 <- krt5_expression %>% as.data.frame() %>% 
  ggplot(aes(x = Krt5))+
  geom_density(fill='blue',alpha=0.4)+theme_light()#+geom_vline(xintercept = threshold,color = 'red',linetype = 2,size=0.7)
p_density_krt5

p_density_krt5_sample <- ggplot(data.frame(krt5_expression), aes(x = Krt5)) +
  ggiraph::geom_density_interactive(
    aes(fill = ifelse(Krt5 > threshold, paste0(">",round(threshold,2)),paste0("<",round(threshold,2)))), 
    color = NA, alpha = 0.7) +
  geom_vline(xintercept = threshold, color = "red", linetype = 2) +
  scale_fill_manual(values = c("#1B9E77", "#D95F02")) +
  theme_light()+labs(fill='Krt5 threshold')
p_density_krt5_sample
krt5_expression$region <- ifelse(krt5_expression$Krt5 > threshold, "Medulla", "Cortex")

sy_6w$region <- ifelse(krt5_expression$Krt5 > threshold, "Medulla", "Cortex")
sy_6w <- AddMetaData(sy_6w,w6_1_image)
sy_6w@meta.data$region_final <- ifelse(sy_6w$region=='Medulla' & sy_6w$area_thy=='Medulla_image', 'Medulla', 'Cortex')





# 03 CD4 & CD8 ----------------------------------------------------------------------


# 03.1 TEC ----------------------------------------------------------------

df_cd48 <- cbind(FetchData(sy_6w,vars = c('region_final','Cd4','Cd8a','Epcam')),GetTissueCoordinates(sy_6w))
tec_threshold <- median(df_cd48$Epcam)
df_cd48 <- df_cd48 %>% dplyr::mutate(celltype=ifelse(Epcam>tec_threshold,'TEC','Thymocyte'))



# 03.2 Classfication ------------------------------------------------------

df_gene48$classification <- as.character(mod2$classification)
mod <- densityMclust(df_gene48[, c("Cd4","Cd8a")], G = 4,modelNames = "VEV")
classification <- c('1'='CD4+CD8+','2'='CD4-CD8-','3'='CD4-CD8+','4'='CD4+CD8-')
df_gene48$sctype <- classification[df_gene48$classification]
df_cd48 <- merge(df_cd48,df_gene48[,c('cell','sctype')],by='cell',sort=F,all.x=T)
df_cd48 <- dplyr::mutate(df_cd48,sctype=ifelse(celltype=='Thymocyte',sctype,'TEC'))
df_cd48_info <- df_cd48 %>% column_to_rownames('cell') %>% dplyr::select(sctype)
sy_6w <- AddMetaData(sy_6w,df_cd48_info)


# 04 Distance with gene expression ------------------------------------------------------

# 04.1 path ---------------------------------------------------------------

# 07 point selection ------------------------------------------------------


library(shiny)
library(ggplot2)
library(dplyr)

ui <- fluidPage(
  titlePanel("Tool for path"),
  sidebarLayout(
    sidebarPanel(
      actionButton("reset", "RESET", icon = icon("trash")),
      actionButton("undo", "UNDO", icon = icon("undo")),
      downloadButton("download", "import"),
      sliderInput("threshold", "select:", 1, 20, 10)
    ),
    mainPanel(
      plotOutput("plot", 
                 click = "plot_click",
                 width = "800px",
                 height = "600px")
    )
  )
)

server <- function(input, output, session) {
  selected <- reactiveVal(data.frame(x = numeric(), y = numeric()))
  path_points <- reactiveVal(data.frame(x = numeric(), y = numeric()))
  
  output$plot <- renderPlot({
    ggplot(df_medulla, aes(x = xlab, y = ylab)) +
      geom_point(color = '#F781BF', size = 2, alpha = 0.6) +
      {
        if(nrow(selected()) > 0) {
          list(
            geom_point(data = selected(), aes(x, y), 
                       color = "red", size = 4, shape = 18),
            geom_path(data = path_points(), aes(x, y), 
                      color = "darkblue", linewidth = 1.5, alpha = 0.8)
          )
        }
      } +
      theme_void() +
      coord_fixed()
  })
  
  observeEvent(input$plot_click, {
    click <- input$plot_click
    
    near_point <- nearPoints(df_medulla, click, 
                             threshold = input$threshold,
                             maxpoints = 1,
                             xvar = "xlab", 
                             yvar = "ylab")
    
    if(nrow(near_point) > 0) {
      new_point <- data.frame(
        x = near_point$xlab,
        y = near_point$ylab
      )
      
      if(!any(selected()$x == new_point$x & selected()$y == new_point$y)) {
        selected(rbind(selected(), new_point))
        path_points(selected())  
      }
    }
  })
  
  observeEvent(input$undo, {
    if(nrow(selected()) > 0) {
      updated <- selected()[-nrow(selected()), ]
      selected(updated)
      path_points(updated)
    }
  })
  
  observeEvent(input$reset, {
    selected(data.frame(x = numeric(), y = numeric()))
    path_points(data.frame(x = numeric(), y = numeric()))
  })

  output$download <- downloadHandler(
    filename = function() {
      paste0("coordinates_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv")
    },
    content = function(file) {
      write.csv(path_points(), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)



# 04.2 distance -----------------------------------------------------------

library(sf)
df_path <- openxlsx::read.xlsx('ST/Sci_Immunol2023/GSM7097913_Ctrl_bio1_tech1_coordinate_medulla.xlsx')
df_gene <- cbind(FetchData(sy_6w,vars = c('region_final','sctype','Cd44','Il2ra','Rag1','Rag2','Gata3','Runx3','S1pr1','Pecam1','Sphk1')),GetTissueCoordinates(sy_6w)) %>% 
  dplyr::mutate(xlab=y,ylab=-x)

df_gene_filter <- cbind(FetchData(sy_6w,vars = c('region_final','sctype','Cd44','Il2ra','Rag1','Rag2','Gata3','Runx3','S1pr1','Pecam1','Sphk1')),GetTissueCoordinates(sy_6w)) %>% 
  dplyr::mutate(xlab=y,ylab=-x) %>% 
  dplyr::filter(sctype %in% c('CD4+CD8-','CD4+CD8+','CD4-CD8-','CD4-CD8+'))


valid_parts <- df_path %>% group_by(part) %>% dplyr::filter(n() >= 3) %>%ungroup()
polygons_sf <- valid_parts %>%
  group_by(part) %>%
  group_split() %>% 
  lapply(function(sub_df) {
    coords <- sub_df %>% dplyr::select(x, y) %>% as.matrix()
    if (!all(coords[1, ] == coords[nrow(coords), ])) {
      coords <- rbind(coords, coords[1, ])
    }
    st_polygon(list(coords))
  }) %>% 
  st_sfc() %>% 
  st_sf(geometry = ., part = unique(valid_parts$part))

points_sf <- st_as_sf(df_gene_filter, coords = c("xlab", "ylab"), crs = st_crs(polygons_sf))
distance_matrix <- st_distance(points_sf, st_boundary(polygons_sf))
points_sf <- points_sf %>% dplyr::mutate(normalized_distance = (min_distance - min(min_distance)) / 
                                           (max(min_distance) - min(min_distance))) 


points_sf$min_distance <- apply(distance_matrix, 1, min)


points_sf_long <- pivot_longer(points_sf,cols = c('Cd44','Il2ra','Rag1','Rag2','Gata3','Runx3'),names_to = 'gene')%>% 
  group_by(gene) %>% 
  dplyr::mutate(normalized_gene = (value - min(value)) / 
                  (max(value) - min(value))) %>% 
  ungroup() %>% 
  dplyr::mutate(gene=factor(gene,levels=c('Cd44','Il2ra','Rag1','Rag2','Gata3','Runx3')))

p <- ggplot(points_sf_long,aes(x=normalized_distance,y=normalized_gene,color=gene))+
  #geom_point(alpha=1,size=1,shape=1)+
  geom_smooth(show.legend = F,method = 'gam',se=F)+
  #geom_smooth(se = F, formula = y ~ splines::ns(x, 10),show.legend = F) +
  ggsci::scale_color_npg()+
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  facet_wrap(~gene,ncol = 6,scale='free_y',dir = 'v')+
  cowplot::theme_cowplot(font_size = 11)+
  theme(strip.background = element_rect(fill='gray90')) +
  labs(x='Spot order along distance from medulla',y='Gene expression',color='')
p
ggsave(p, filename = file.path(outdir,"Figure/Figure1.distance_gene_expression.pdf"), width = 12, height = 6)