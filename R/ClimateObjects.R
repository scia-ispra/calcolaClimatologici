#15 Febbraio 2017
library("xts")
library("lubridate")

#############################################################
#funzione ClimateData:crea un oggetto ClimateData
#
#Un oggetto ClimateData è un oggetto xts con classe ClimateData e una delle tre seguenti classi:
#daily/monthly/yearly. L'oggetto così creato ha un attributo "parametro" ricavabile mediante
#la funzione xtsAttribute. L'attributo "parametro" permette di stabilire che tipo di aggregazione
#utilizzare (somma o media)
#
#Input: Un dataframe di serie giornaliere, mensili o annuali in cui compare almeno la colonna yy (anni)
#Altre colonne opzionali sono le colonne mm (mesi) e le colonne dd (giorni). 

#Output: un oggetto di tipo ClimateData
#Se i dati sono giornalieri, viene associato al dato anno-mese-giorno
#Se i dati sono mensili, viene associato al dato anno-mese-01
#Se i dati sono annuali viene associato al dato anno-01-01

ClimateData<-function(x,...)
{
  UseMethod("ClimateData",x)  
}

ClimateData.data.frame<-function(x,param=NULL)
{

  #le colonne che identificano anno mese e giorno possono essere indicate come:
  #year, month, day
  #yy mm dd
  #le due forme possono essere mischiate
  
  stopifnot(is.data.frame(x) && ("yy" %in% names(x) || "year" %in% names(x)) )
  #se dd è presente deve essere presente anche mm
  if("dd" %in% names(x) || "day" %in% names(x) ) stopifnot("mm" %in% names(x) || "month" %in% names(x))

  if(!is.character(param) || !(param %in% c("pr","tas","tasmax","tasmin","prpc","tmax","tmin","tmean","prcp") )) 
    stop("Parametro param mancante o non corretto")

  #ricodifichiamo year, month e day se presenti come nomi colonne
  names(x)[names(x) %in% c("year")]<-"yy"
  names(x)[names(x) %in% c("month")]<-"mm"
  names(x)[names(x) %in% c("day")]<-"dd"
  
    
  #altra possibilità di nomi
  if((param=="prpc") || (param=="prcp")){
    param<-"pr"    
  }else if(param=="tmax"){
    param<-"tasmax"
  } else if(param=="tmin"){
    param<-"tasmin"
  }
  
  names(x)->nomi
  lStr<-list("daily"=c("yy","mm","dd"),"monthly"=c("yy","mm"),"yearly"=c("yy"))
  
  if("dd" %in% nomi){
    "daily"->tipoDati
    yy<-x$yy
    mm<-x$mm
    dd<-x$dd
  }else if("mm" %in% nomi){
    "monthly"->tipoDati
    yy<-x$yy
    mm<-x$mm
    dd<-rep("01",length(mm))    
  }else{
    "yearly"->tipoDati
    yy<-x$yy
    mm<-rep("01",length(yy))
    dd<-rep("01",length(yy))        
  }
  
  #switch(tipoDati,"daily"=4,"monthly"=3,"yearly"=2)->colonna
  lStr[[tipoDati]]->stringa
  x[,!names(x) %in% stringa,drop=FALSE]->soloDati
  
  #caso assurdo che passo solo le colonne dei tempi
  if(!ncol(soloDati)) return(NULL)
  
  #variabile da parrase a xts per ordinare i dati
  tempo<-as.Date(paste(yy,mm,dd,sep="-"),format="%Y-%m-%d")

  #converte tutti i dati in double. Necessario perchè readr potrebbe leggere alcune colonne come character. 
  #Questo si verifica quando i primi dati su cui readr verifica il tipo dato sono tutti NA.
  #Avere dati double e non character è necessario per poi creare serie numeriche con xts.
  soloDati %>% purrr::map_if(is.character,as.double)->soloDati_dbl
  rm(soloDati)
  xts(x=as.data.frame(do.call("cbind",soloDati_dbl),optional=FALSE),order.by=tempo)->xDati
  class(xDati)<-c("ClimateData",tipoDati,class(xDati))
  xtsAttributes(xDati)<-list("parametro"=param)

  xDati

}#fine funzione ClimateDataset
#######################################################


########################################################
#
# una classe ClimateData la possiamo creare direttamente da un oggetto xts 
ClimateData.xts<-function(x,param=NULL)
{
  if(!is.character(param) || !(param %in% c("pr","tas","tasmax","tasmin","prpc","tmax","tmin") )) 
    stop("Parametro param mancante o non corretto")
  
  #altra possibilità di nomi
  if(param=="prpc"){
    param<-"pr"    
  }else if(param=="tmax"){
    param<-"tasmax"
  } else if(param=="tmin"){
    param<-"tasmin"
  }  
  

  tryCatch({
    periodicity(x)$scale
  },error=function(e){
    
    #nel caso dei climatologici annuali, periodicity fallisce
    if(nrow(x)==1){
      "yearly" 
    }else{
      ""
    }
    
  })->periodicita
      
  if(!periodicita %in% c("daily","monthly","yearly")) stop("Periodicità non riconosciuta")
  
  class(x)<-c("ClimateData",periodicita,class(x))
  xtsAttributes(x)<-list("parametro"=param)
  
  x
  
}  

#########################################################


#verifica se si tratta di un oggetto ClimateData
is.ClimateData<-function(x)
{
  ifelse(is.xts(x) && "ClimateData" %in% class(x),TRUE,FALSE)
  
}#fine is.ClimateData

#########################################################


as.data.frame.ClimateData<-function(x)
{
  as.data.frame(coredata(x),optional=TRUE)->y
  periodicity(x)$scale->periodicita

  data.frame(yymmdd=as.character(index(x))) %>% tidyr::separate(.,col=yymmdd,into=c("yy","mm","dd"),sep="-")->yymmdd
  
  
  if(periodicita=="yearly"){
    
    y %>% mutate(yy=yymmdd$yy) %>% dplyr::select(yy,everything())
    
  }else if(periodicita=="monthly"){
    
    y %>% mutate(yy=yymmdd$yy,mm=yymmdd$mm) %>% dplyr::select(yy,mm,everything())
    
  }else if(periodicita=="daily"){
    
    y %>% mutate(yy=yymmdd$yy,mm=yymmdd$mm,dd=yymmdd$dd) %>% dplyr::select(yy,mm,dd,everything())
    
  }else{
    
    stop("Periodicità non riconosciuta")
    
  }
  
}#fine as.data.frame


#Verifica se un oggetto climate data sia composto di soli NA
#Ad esempio se calcolo i climatologici mensili, il calcolo delle anomalie
#deve everificare che il climatologico non sia tutto NA.
allNA<-function(x){
  all(purrr::map(is.na(x),all) %>% flatten_lgl())
}


