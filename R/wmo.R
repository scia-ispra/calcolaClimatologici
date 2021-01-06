#currying di alcune funzioni in ClimateData.R utilizzando i parametri richiesti dal
#WMO
source("ClimateData.R")

#aggregazione da giornalieri a mensili: wmo_aggregaCD.daily restituisce
#una funzione con i parametri settati in base al parametro precipitazione/temperatura
wmo_aggregaCD<-function(parametro)
{
  
  stopifnot(is.character(parametro))
  stopifnot(parametro %in% c("pr","tasmax","tasmin"))
  
  if(parametro=="pr"){
    
    purrr::partial(.f=aggregaCD.daily,max.na=0,rle.check=FALSE,max.size.block.na=0)    
    
  }else{
    
    purrr::partial(.f=aggregaCD.daily,max.na=5,rle.check=TRUE,max.size.block.na=3)    
    
  }

}#fine wmo_aggregaCD.daily


#calcolo dei climatologici mensili: 
#i parametri si riferiscono al climatologico calcolato su 30 anni
#Il WMO ammette sei mesi mancanti e al massimo 3 NA consecutivi
wmo_climatologiciMensili<-purrr::partial(.f=climatologiciMensili,max.na=6,rle.check=TRUE,max.size.block.na=3)


#calcolo del climatologico annuale secondo i dettami del WMO, ovvero:
#nessun mese mancante.
wmo_climatologiciAnnuali<-function(climatol){
  
  checkClimatol(climatol)    
  aggregaCD.yearly(x=climatol,ignore.par = FALSE,max.na = 0,rle.check = FALSE,max.size.block.na = 0,seasonal = FALSE)
  
}