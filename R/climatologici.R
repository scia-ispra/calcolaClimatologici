rm(list=objects())
library("tidyverse")
library("seplyr")
library("skimr")
library("sf")
library("geojson")
source("ClimateData.R")
source("wmo.R")
source("ClimateObjects.R")

annoI<-1991
annoF<-2020

REGIONE<-"FVG"

#se TRUE verifica che gli ultimi anni ULTIMI.ANNI della serie non sia tutta di NA, questo perche'
#e' importante che il climatologico sia rappresentativo degli ultimi anni e non solo dell'inizio
#del trentennio
VERIFICA.ULTIMI.ANNI<-TRUE

#ULTIMI.ANNI: quanti anni da verificare
ULTIMI.ANNI<-2

listaParametri<-c("Prec","Tmax","Tmin")


#QUANTI ANNIMANCANTI accetto per costruire il valore climatologico? Il valore normale standard 
#rihiede al massimo 6 ANNIMANCANTI
listaAnniMancanti<-c(6,10,15)

creaCalendario<-function(annoI,annoF){
  
  if(missing(annoI)) stop("Anno inizio mancante")
  if(missing(annoF)) stop("Anno fine mancante")
  
  as.integer(annoI)->annoI
  as.integer(annoF)->annoF
  
  stopifnot(annoI<=annoF)
  
  seq.Date(from=as.Date(glue::glue("{annoI}-01-01")),to=as.Date(glue::glue("{annoF}-12-31")),by="day")->yymmdd
  
  tibble(yymmdd=yymmdd) %>%
    tidyr::separate(yymmdd,into=c("yy","mm","dd"),sep="-") %>%
    dplyr::mutate(yy=as.integer(yy),mm=as.integer(mm),dd=as.integer(dd))
  
}


#lettura anagrafica
list.files(pattern="^reg.+csv$")->fileAna
stopifnot(length(fileAna)==1)

read_delim(fileAna,delim=";",col_names=TRUE)->ana

#la directory dati contiene il risultato dei controlli spaziali
list.files(path = "./dati",pattern = "^.+txt",full.names = TRUE)->ffile

stopifnot(length(ffile)>0)

creaCalendario(annoI=annoI,annoF=annoF)->calendario

purrr::walk2(listaParametri,listaAnniMancanti,.f=function(PARAMETRO,ANNIMANCANTI){
  
  #Prcp e' il nome della precipitazione nel file dei dati di input
  if(grepl("^P.+",PARAMETRO)) PARAMETRO<-"Prcp"
  
  if(grepl("^Tm..",PARAMETRO)){
    
    MAX.NA.GIORNALIERI<-10
    MAX.SIZE.BLOCK.NA<-4
    
  }else{
    
    MAX.NA.GIORNALIERI<-0
    MAX.SIZE.BLOCK.NA<-0
    
  }

    purrr::map(ffile,.f=function(nomeFile){ 
      
      codice<-str_remove(str_remove(nomeFile,pattern="\\.txt$"),pattern="^.+\\/")
      stopifnot(nchar(codice)>0)
      
      TIPI<-cols(year=col_integer(),month=col_integer(),day=col_integer(),.default=col_double())
      
      read_delim(nomeFile,
                 delim=",",
                 col_names=TRUE,
                 col_types = TIPI)->dati
      
      dati %>%
        seplyr::select_se(c("year","month","day",tolower(PARAMETRO))) %>%
        rename(yy=year,mm=month,dd=day) %>%
        seplyr::rename_se(c("value":=tolower(PARAMETRO))) %>%
        filter(yy>=annoI & yy<=annoF)->subDati
      
      if(!nrow(subDati)) return(NULL)
      
      left_join(calendario,subDati,by=c("yy"="yy","mm"="mm","dd"="dd"))->subDati
    
      ClimateData(x=subDati,param = tolower(PARAMETRO))->clSerie
      aggregaCD(clSerie,max.na = MAX.NA.GIORNALIERI,rle.check = TRUE,max.size.block.na = MAX.SIZE.BLOCK.NA)->clSerieMonthly 
      
      #verifichiamo che gli ultimi anni non siano tutti NA
      if(VERIFICA.ULTIMI.ANNI){
        
        as.data.frame(clSerieMonthly)->temp
        
        
        temp %>%
          mutate(yy=as.integer(yy)) %>%
          filter(yy>=(annoF-ULTIMI.ANNI+1))->temp
        
        if(all(is.na(temp$value))){
          message(glue::glue("Serie {codice} con NA negli ultimi anni, non valida"))
          return()
        }
        
      }
      
      
      climatologiciMensili(clSerieMonthly,yearS = annoI,yearE = annoF,max.na=ANNIMANCANTI,rle.check = FALSE,max.size.block.na = ANNIMANCANTI)->clMensili
      aggregaCD(clMensili,ignore.par = FALSE,max.na = 0,rle.check = TRUE,max.size.block.na = 0,seasonal = TRUE)->clAnnuali
    
      as_tibble(as.matrix(clMensili))->clMensili
      
      if(all(is.na(clMensili))) return(NULL)
      
      as_tibble(as.matrix(clAnnuali))->clAnnuali
    
      names(clMensili)<-codice
      names(clAnnuali)<-codice
      
      list(clMensili,clAnnuali)
        
    })->listaOut
    
    purrr::compact(listaOut)->listaOut
    
    if(!length(listaOut)){
      message(glue::glue("Nessun valore climatologico calcolato: {PARAMETRO}, {ANNIMANCANTI}!")) 
      return()
    } 
    
    purrr::map_dfc(listaOut,1)->valoriMensili
    purrr::map_dfc(listaOut,2)->valoriAnnuali
    
    valoriAnnuali %>%
      mutate(yy=glue::glue("{annoI}-{annoF}")) %>%
      dplyr::select(yy,everything()) %>%
      gather(key="SiteID",value="climatologico",-yy) %>%
      mutate(SiteID=as.numeric(SiteID)) %>%
      filter(!is.na(climatologico))->gvaloriAnnuali
    
    left_join(gvaloriAnnuali,ana)->gvaloriAnnuali
    
    gvaloriAnnuali %>% 
      dplyr::select(yy,SiteID,climatologico,SiteCode,SiteName,Latitude,Longitude,matches("Elevation"))->gvaloriAnnuali
    
    gvaloriAnnuali$regione<-REGIONE
    
    write_delim(gvaloriAnnuali,file=glue::glue("{PARAMETRO}_{REGIONE}_{annoI}_{annoF}_annuali_annimancanti{ANNIMANCANTI}.csv"),delim=";",col_names = TRUE)
    
    
    # st_as_sf(dati,coords=c("Longitude","Latitude"),crs=4326)->puntiStazione
    # geojson::geo_write(as.geojson(as_Spatial(puntiStazione)),file=glue::glue("{PARAMETRO}_annuali.geojson"))
    
    valoriMensili %>%
      mutate(yy=glue::glue("{annoI}-{annoF}")) %>%
      mutate(mm=1:12) %>%
      dplyr::select(yy,mm,everything()) %>%
      gather(key="SiteID",value="climatologico",-yy,-mm) %>%
      mutate(SiteID=as.numeric(SiteID))->gvaloriMensili
    
    left_join(gvaloriMensili,ana)->gvaloriMensili
    
    gvaloriMensili %>% 
      dplyr::select(yy,SiteID,climatologico,SiteCode,SiteName,Latitude,Longitude,matches("Elevation"))->gvaloriMensili
    
    gvaloriMensili$regione<-REGIONE
    
    write_delim(gvaloriMensili,file=glue::glue("{PARAMETRO}_{REGIONE}_{annoI}_{annoF}_mensili_annimancanti{ANNIMANCANTI}.csv"),delim=";",col_names = TRUE)
      
    #st_as_sf(subDati,coords=c("Longitude","Latitude"),crs=4326)->puntiStazione
    #geojson::geo_write(as.geojson(as_Spatial(puntiStazione)),file=glue::glue("{PARAMETRO}_{month.abb[mese]}.geojson"))

}) #fine walk su listaParametri    
