#11 Novembre 2016: la funzione ClimateData per creare un oggetto ClimateData ora accetta anche 
#come colonne valide per identificare anni mesi giorni e i relativi parametri i nomi utilizzati
#nei file prodotti dai controlli di qualità sulle serie giornaliere: year month day prcp tmax tmin
#
#ClimateDataset object per gestire le serie annuali e mensili di precipitazione e temperatura
library("tidyverse")
library("lazyeval")
library("xts")
library("stringr")
library("lubridate")
source("ClimateObjects.R")

stagioni<-list("primavera"=c(3,4,5),"estate"=c(6,7,8),"autunno"=c(9,10,11),"inverno"=c(1,2,12))

#######################################################

    #verificaBlocchiNA:
    #- Input: un vettore numerico e max.size.block.na (lunghezza massima accettabile di un blocco di NA)
    #- Output: lista di due elementi:
    #  - un elemento "ris" boolean pari a TRUE se il vettore vserie contiene almeno un blocco di NA di lunghezza superiore a max.size.block.na
    #  - un elemento "nb" (number of blocks) riporta il numero di blocchi che superano max.size.block.na
  
    # Per ris==FALSE, nb deve essere pari a zero. 
    verificaBlocchiNA<-function(vserie,max.size.block.na=3)
    {
      which(is.na(vserie))->index.na
      #per identificare blocchi di NA contigui prima di tutto sostituiamo NA con un valore
      #impossibile: -3333. Necessaria questa trasformazione in quanto rle riporta ogni NA come
      #un blovcco di dimensione 1, anche quando contigui
      vserie[index.na]<- -3333
      rle(vserie)->ris.rle
      ris.rle$lengths[ris.rle$values==-3333]->repeatedNA
      #ho trovato un blocco di NA più lungo di quanto permesso
      if(any(repeatedNA>max.size.block.na)) return(list("ris"=TRUE,"nb"=length(repeatedNA[repeatedNA>max.size.block.na])))
      
      return(list("ris"=FALSE,"nb"=length(repeatedNA[repeatedNA>max.size.block.na])))
      
    }#fine verificaBlocchiNa
    
#######################################################

#aggregaCD: aggrega i dati climatici
#Input: oggetto ClimateData con dati giornalieri
#Output: oggetto ClimateData con dati mensili. La temperatura viene aggregata mediante media, la precipitazione mediante somma

#Input: oggetto ClimateData con dati mensili
#Output: oggetto ClimateData con dati annuali. La temperatura viene aggregata mediante media, la precipitazione mediante somma

#Input: oggetto ClimateData con dati annuali
#Output: oggetto ClimateData con dati annuali (stesso oggetto)
aggregaCD<-function(x,...)
{
  stopifnot(is.ClimateData(x))
  UseMethod("aggregaCD",x)
  
}  

#default
aggregaCD.default<-function(x) { stop("Non so aggregare un oggetto non Climate Data") }

#La funzione aggrega viene utilizzata sia da apply.monthly che da apply.yearly per 
#aggregare (somma o media) i dati di precipitazione/temperatura. Nel caso dell'aggregazione
#annuale è necessario passare seasonal=TRUE in modo che venga verificata la validità di tutte le
#stagioni in un anno. Una stagione si considera valida quando sono validi almeno due mesi su tre
aggrega<-function(serie,aggFun=myFun,max.na,rle.check,max.size.block.na,seasonalCheck)
{
  #aggrega lavora su oggetti xts. Ad esempio, per passare dal mensile all'annuale devo
  #fare un controllo sulle stagioni quindi ho bisogno dell'informazione temporale
  stopifnot(is.xts(serie))
  stopifnot(is.numeric(max.na) || is.logical(rle.check) || is.numeric(max.size.block.na) || is.logical(seasonalCheck))
  which(is.na(serie))->index.na
  length(index.na)->numeroNA
  #fallito già il primo controllo restituisco NA
  if(numeroNA > max.na) return(NA)
  
  #se ci sono NA e voglio fare il controllo sulla presenza di blocchi di NA contigui
  #eseguo il codice in if
  if(rle.check && numeroNA){ 
    risBlocchiNA<-verificaBlocchiNA(c(coredata(serie)),max.size.block.na=max.size.block.na)$ris
    if(risBlocchiNA) return(NA)
  }  
  
  #Per aggregare i mensili a livello annuale dobbiamo sincerarci che tutte le stagioni siano valide
  #Una stagione è valida quando almeno due mesi su tre non sono NA
  if(seasonalCheck){
    if(numeroNA>=3){
      .indexmon(serie[index.na])+1->mesi
      #verifichiamo che ci siano almeno due stagioni valide
      unlist(lapply(stagioni,FUN=function(ss){length(mesi[mesi %in% ss])}))->len.out
      #len.out mi dice, per ogni stagione, quanti NA ho trovato. Se trovo un valore maggiore di 1
      #significa che almeno un stagione ha due NA ovvero non è valida, quindi non valido l'anno
      if(any(len.out>1)) return(NA)
    }
  }

  aggFun(serie,na.rm=TRUE)
  
}  


#max.na: massimo numero di na all'interno di un mese o all'interno di NA che sono disposto a tollerare
#rle.check: se TRUE verifica che all'interno del mese/anno non vi siano blocchi di NA di lunghezza
#           superiore a "max.size.block.na". Il controllo si rivela utile per la temperatura
#           dove il numero di NA che siamo disposti ad accettare è maggiore di quello della precipitazione.
#           Nel caso della precipitazione ponendo max.na=3 rende superfluo il controllo "rle.check"

aggregaCD.daily<-function(x,max.na,rle.check,max.size.block.na)
{

  xtsAttributes(x)$parametro->parametro

  #i giornalieri si aggregano diversamente a seconda del parametro pr/tas..
  ifelse(parametro=="pr",sum,mean)->funzione
  
  #in base al tipo di dati di x determino il nuovo tipo di dati che verrà prodotto da
  #aggregaCD. Se si tratta di dati daily otterrò dati monthly. Nel caso di dati monthly
  #otterrò dati yearly. In base al tipo di dati di output determino anche la funzione
  #da applicare sulle serie ovvero: apply.monthly (per i dati daily) o apply.yearly (per i dati monthly)

  purrr::map(1:ncol(x),.f=~(apply.monthly(x[,.],FUN=aggrega,aggFun=funzione,max.na=max.na,
                  rle.check=rle.check,max.size.block.na=max.size.block.na,seasonal=FALSE)

  ))->xx
  
  #cbind aggiunge "X" davanti ai codici numerici

  reduce(xx,cbind.xts) %>% ClimateData(x=.,param=parametro)->xx2
  rm(xx)
  
  names(xx2)<-stringr::str_replace_all(names(xx2),"^X","")
  
  xx2      
  
  
}#fine aggregaCD  


aggregaCD.monthly<-function(x,ignore.par,max.na,rle.check,max.size.block.na,seasonal)
{
  
  if(missing(ignore.par) || missing(max.na) || missing(rle.check) || missing(max.size.block.na) || missing(seasonal)) stop("Parametri mancanti")    
  stopifnot(is.logical(ignore.par) && is.numeric(max.na) || is.logical(rle.check) || is.numeric(max.size.block.na) || is.logical(seasonal))  
  xtsAttributes(x)$parametro->parametro
  
  #Per passare dai valori mensili agli annuali medio (tas) o sommo (pr). Questo
  #non vale quando calcolo l'anomalia media annuale per la
  #precipitazione: media appunto, non sommo le anomalie dei mesi
  ifelse(ignore.par,mean,ifelse(parametro=="pr",sum,mean))->funzione      

  #in base al tipo di dati di x determino il nuovo tipo di dati che verrà prodotto da
  #aggregaCD. Se si tratta di dati daily otterrò dati monthly. Nel caso di dati monthly
  #otterrò dati yearly. In base al tipo di dati di output determino anche la funzione
  #da applicare sulle serie ovvero: apply.monthly (per i dati daily) o apply.yearly (per i dati monthly)


  purrr::map(1:ncol(x),.f=~(apply.yearly(x[,.],FUN=aggrega,aggFun=funzione,max.na=max.na,
                 rle.check=rle.check,max.size.block.na=max.size.block.na,seasonal=seasonal)    
  )) %>% reduce(cbind)->ris

  names(ris)<-names(x)  
  ClimateData(x=ris,param=parametro)

}#fine aggregaCD  


#aggregaCD.yearly: a differenza delle corrispettive stazioni daily e monthly, aggregaCD.yearly
#prende dati annuali e restituisce ancora un oggetto con classe "yearly". Questa funzione viene
#utilizzata per creare un valore singolo (somma o media) di valori annuali
#
#Nel caso "yearly" dobbiamo prevedere che per la precipitazione sono possibili due operazioni:
# - somma (se calcolo l'aggregato annuale)
# - media: l'operazione media viene utilizzata per costruire il climatologico mensile (ad esempio il climatologico del mese di gennaio)


aggregaCD.yearly<-function(x,ignore.par,max.na,rle.check,max.size.block.na,seasonal)
{

  if(missing(ignore.par) || missing(max.na) || missing(rle.check) || missing(max.size.block.na) || missing(seasonal)) stop("Parametri mancanti")  

  xtsAttributes(x)$parametro->parametro

  #se ignore.par==TRUE, dobbiamo fare la media. Questo serve per calcolare il
  #climatologico mensile della precipitazione
  #In ogni altro caso saraà il parametro a determinare la funzione di aggregazione
  if(ignore.par){
    mean->funzione
  }else{  
    if(parametro=="pr"){
      sum->funzione
    }else{  
      mean->funzione
    }  
  }  

  #in base al tipo di dati di x determino il nuovo tipo di dati che verrà prodotto da
  #aggregaCD. Se si tratta di dati daily otterrò dati monthly. Nel caso di dati monthly
  #otterrò dati yearly. In base al tipo di dati di output determino anche la funzione
  #da applicare sulle serie ovvero: apply.monthly (per i dati daily) o apply.yearly (per i dati monthly)
  purrr::map(1:ncol(x),.f=~(aggrega(x[,.],aggFun=funzione,max.na=max.na,
            rle.check=rle.check,max.size.block.na=max.size.block.na,seasonal=seasonal))) %>%
            flatten_dbl()->ris
  
  names(ris)<-names(x)  
  #l'anno di riferimento è l'anno a metà tra l'inizio e la fine del periodo coperto
  #con mese "01" e giorno "01".
  yearS<-lubridate::year(index(x)[1])
  yearE<-lubridate::year(index(x)[nrow(x)])  
  mese<-lubridate::month(index(x)[1])
  xts(x=as.data.frame(t(ris),optional=TRUE),order.by =as.Date(paste(floor((yearS+yearE)/2),mese,"01",sep="-"),format="%Y-%m-%d"))->xris
  rm(ris)
  
  class(xris)<-c("ClimateData","yearly",class(xris))
  xtsAttributes(xris)<-list("parametro"=parametro)
  
  #ClimateData(x=xris,param=parametro) <--- non posso utilizzare ClimateData.xts perchè ho un singolo valore
  #e il calcolo della periodicita fallisce
  
  xris
    
}#fine aggregaCD 


##############################################
#checkSeriesValidity: verifica la validità di una serie annuale secondo i controlli applicati
#per la selezione delle serie SCIA utilizzate per la stima dei trend.
# 
# - L'86% dei dati deve essere valido
# - Un blocco di NA può al massimo essere lungo 4 anni
# - Gli ultimi 3 anni devono essere disponibili
#
# Input: un oggetto ClimateData con attributo "tipoDati"=="yearly"
# Output: lo stesso oggetto con le sole serie valide

checkSeriesValidity<-function(x,percentualeAnniPresenti=86,max.size.block.na=4)
{
  stopifnot(is.ClimateData(x) && "yearly" %in% class(x))
  stopifnot(is.numeric(percentualeAnniPresenti) && is.numeric(max.size.block.na))
  
  controllo<-function(serie){
    
    length(serie)->lenSerie
    
    #la serie deve avere almeno il "percentualeDati"% dei dati
    floor(lenSerie*(percentualeAnniPresenti/100) )->lunghezzaMinimaSerie
    which(!is.na(serie))->index.not.na
    if(length(index.not.na)<lunghezzaMinimaSerie) return(FALSE)
    
    #i blocchi di dati mancanti possono essere al più lunghi 4 anni
    verificaBlocchiNA(serie,max.size.block.na=max.size.block.na)->risBlocchiNA
    if(risBlocchiNA$ris) return(FALSE)
    
    #le serie devono avere gli ultimi tre anni non tutti NA
    if(all(is.na(serie[(lenSerie-2):lenSerie]))) return(FALSE)
    
    return(TRUE)
    
  }  

  #risControllo: vettore di FALSE (serie non valida) o TRUE (serie valida) con cui filtriamo
  #l'oggetto xts
  as.data.frame(coredata(x),optional=TRUE) %>% purrr::map_lgl(controllo)->risControllo
  
  if(all(risControllo==FALSE)) return(NULL)
  #restituisco x filtrato con le sole serie valide
  x[,names(x)[risControllo] ]

}#fine checkSeries



checkSeriesValidity2<-function(x,percentualeAnniPresenti,max.size.block.na,minLen=10)
{
  stopifnot(is.ClimateData(x) && "yearly" %in% class(x))
  stopifnot(is.numeric(minLen) || minLen>=0)
  stopifnot(is.numeric(percentualeAnniPresenti) && is.numeric(max.size.block.na))
  
  controllo<-function(serie){
      
    imin<-1
    imax<-length(serie)  
    
    ii<-1

    while(TRUE){

      #questo if serve a gestire i casi in cui i controlli di qualita
      #falliscono e imax e imin puntano a valori non NA. Se non ci fosse
      #questo codice il programma girerebbe senza fermarsi. ii>1 perchè
      #nel primo ciclo questo controllo nn va fatto
      if(ii>1 && forzaIf){
        forzaIf<-FALSE
        which(is.na(serie))->index.na
        stopifnot(length(index.na)!=0)
        #degli NA ci debbono essere per forza altrimenti i controlli precedenti
        #avrebbero restituito serie di NA (serie non abbstanza lunga,minore di minLen)
        #oppure la serie vera e propria (serie valida)
        
        #imin potrebbe essersi spostato. lo stesso imax.
        #devo escludere gli NA sui cui già si sono mossi imin e/o imax
        #altrimenti loop infinito
        index.na[index.na>imin & index.na<imax]->index.na
        #anche in questo caso, index.na deve essere non nulla altrimenti i controlli
        #avrebbero fatto terminare la funzione
        
        min(index.na)->minNA
        max(index.na)->maxNA

        #a partire da imin cerco un NA, lo stesso a partire da imax
        #vedo la distanza tra l'NA e imin e imax. Mi sposto sull'Na
        #che dista meno da imin o imax. Voglio infatti mantenere la 
        #parte della serie che ha più dati non NA
        if((minNA-imin)<=(imax-maxNA)){
          imin<-(minNA)
        }else{
          imax<-(maxNA)        
        }      
      }#su ii>1      

      (imax-imin+1)->lenSubSerie
      #se imin/imax definiscono una serie di partenza inferiore a minLen inutile
      #fare altri controlli, la serie va esclusa senza spendere cicli. Questo controllo non basta:
      #mi dice solo che la serie tra imin e imax è lunga almeno minLen. Questo non esclude che sia una serie
      #con molti NA per cui il numero di dati non NA non raggiunge minLen. Questo primo controllo ha solo lo scopo 
      #di scartare serie corte che non meritano cicli di calcolo
      if(lenSubSerie<minLen){serie[1:length(serie)]<-NA;return(serie)}
      
      #imin e imax puntano a due valori NA, restringo da entrambi i latila serie
      if(is.na(serie[imin]) && is.na(serie[imax])){
        imin<imin+1
        imax<-imax-1
        next
      }else if(is.na(serie[imin])){
        imin<-imin+1
        next
      }else if(is.na(serie[imax])){ 
        imax<-imax-1
        next
      } 
    
      #subSerie: gli estremi ora sono due non NA
      serie[imin:imax]->subSerie
      
      #ok ora subSerie identifica una porzione di serie tra due valori non NA e lunga almeno
      #minLen. 
      #Due domande: 
      #1) I valori non NA in subSerie sono almeno pari a minLen?
      #2) Ok ho almeno minLen dati validi. La serie comunque potrebbe essere formata da tantissimi NA
      #per cui la percentuale di dati validi non rispetti il vincolo dato da "percentualeDati".
      #Se ci fossero troppi NA devo continuare a cercare la sottoserie che contiene almeno minLen anni validi
      #ma senza avere al suo interno troppi NA  
      floor(lenSubSerie*(percentualeAnniPresenti/100) )->lunghezzaMinimaSerie
      length(subSerie[!is.na(subSerie)])->lenSenzaNA
      #se togliendo gli NA la serie è minore di lunghezzaMinimaSerie inutile andare
      #avanti, a ogni ciclo la serie si restringe, quindi a un giro successivo nn
      #soddisferei questa condizione
      if(lenSenzaNA<lunghezzaMinimaSerie || lenSenzaNA<minLen){serie[1:length(serie)]<-NA;return(serie)}

      #i blocchi di dati mancanti possono essere al più lunghi 4 anni
      verificaBlocchiNA(subSerie,max.size.block.na=max.size.block.na)->risBlocchiNA
      #anche questo controllo fallito..ritorno su.
      #Devo incrementare ii perchèse  ora imin e imax puntano a due NA, da li non si schiodano
      #ma siccome i controlli sono falliti la funzione non termina. Devo forzare i min e imax
      #a muoversi verso un nuono NA. L'idea è eliminare parte dei dati validi
      #per trovare una sottoserie su cui rifare i controlli
      
      if(risBlocchiNA$ris){ii<-ii+1;forzaIf<-TRUE; next}

      #searrivo qui tutti i controlli sono passati!!! serie buona
      break
                    
    }#fine su while

    # ifelse(imin-1>0,imin-1,1)->imin
    # ifelse(imax+1>length(serie),imax,imax+1)->imax 
    
    #la parte fuori da imin:imax è la parte che nn permette di soddisfare i criteri di qualita
    #la devo annullare
    serie[-c(imin:imax)]<-NA
    return(serie)
    
  }  
  
  #risControllo: vettore di FALSE (serie non valida) o TRUE (serie valida) con cui filtriamo
  #l'oggetto xts

  as.data.frame(coredata(x),optional=TRUE) %>% purrrlyr::dmap(controllo)->risControllo

  risControllo %>% purrr::map_lgl(.f=function(x){any(!is.na(x))})->vvv
  if(all(vvv==FALSE)) return(NULL)

  #restituisco x filtrato con le sole serie valide
  x[,names(x)[vvv] ]
  
}#fine checkSeries


################################################
#######
#funzione climatologiciMensili: utilizza la funzione "calcolaClimatologico" per il calcolo
#dei valori climatologici mensili e annuali secondo i criteri del WMO (vedi seguito).
#
#Calcolo dei valori climatologici secondo i criteri del WMO.
# 1) Dai dati giornalieri (oggetto ClimateData daily) calcolare l'aggregato mensile con aggregateCD 
# passando l'opzione strict.wmo=TRUE che fissa max.na e max.block.size.na pari a quanto richiesto dal WMO
# 2) Utilizzare la funzione climatologiciMensili per ottenere il climatologico di ogni mese. I dati
# così ottenuti sono oggetti ClimateData monthly.
# 3) Utilizzare aggregaeCD (per l'esattezza: aggregaCD.monthly) con max.na=0 per ottenere il climatologico annuale
#
# Il climatologico mensile secondo il WMO:
# -Temperatura: l'80% dei dati debbono essere disponibili (ovvero 24 dati su trentennio)
#               al massimo tre anni consecutivi mancanti
# Precipitazione: l'80% dei dati debbono essere disponibili (ovvero 24 dati su trentennio)
#               al massimo tre anni consecutivi mancanti
#
#I climatologico mensile sia per la precipitazione che per la temperatura si ottiene per media
#
#A livello annuale invece i 12 valori climatologici debbono essere tutti presenti. Il climatologico
#annuale si ottiene come media per latmperatura e come somma per la precipitazione
climatologiciMensili<-function(x,yearS,yearE,max.na,rle.check,max.size.block.na)
{
  stopifnot(is.ClimateData(x) && "monthly" %in% class(x))
  
  if(missing(yearS) || missing(yearE)) stop("Specificare anno inizio/fine calcolo climatologico")
  stopifnot(is.numeric(yearS) && is.numeric(yearE))
  

  if(missing(max.na) || missing(rle.check) || missing(max.size.block.na))
    stop("Specificare max.na, rle.check e max.size.block.na")

  stopifnot(is.numeric(max.na) && is.logical(rle.check) && is.numeric(max.size.block.na))      

  numeroMinimoAnni<-(yearE-yearS+1)-max.na
  stopifnot(numeroMinimoAnni>0)
  
  x[paste0(yearS,"/",yearE)]->subx
  
  #subDati è vuoto
  stopifnot(nrow(subx)!=0)
  
  #Altro controllo da fare: se i dati partono dal 1961 e voglio calcolare
  #il climatologico 1951/1980 ammettendo al massimo 6 NA devo prima di tutto
  #verificare che subx contenga almeno 24 anni distinti. Se non faccio questo controllo
  #ottengo, ad esempio, che subx ha i valori dal 1961 al 1980 tutti completi per cui
  #aggregaCD fa le mie aggregazioni ...ma su 20 anni di dati!! Quando io ne devo avere
  #almeno 24!! Quindi: verifichiamo che in subx ci siano 24 anni distinti altrimenti
  #restituisco NA
  unique(year(index(subx)))->anniDistinti
  
  if(length(anniDistinti)<numeroMinimoAnni){
    ncol(x)->colonne
    floor((yearE+yearS)/2)->annoOut
    seq.Date(from=as.Date(paste0(annoOut,"-01-01")),to=as.Date(paste0(annoOut,"-12-01")),by="month")->indice
    matrix(rep(NA,colonne*12),nrow = 12,ncol=colonne)->matriceNA
    xts(matriceNA,order.by = indice)->climatolNA
    names(climatolNA)<-names(x)
    return(ClimateData(climatolNA,param=attributes(x)$parametro))
  }
    
  
  purrr::map(0:11,.f=function(mese){
    
    subx[.indexmon(subx)==mese,]->serieMeseAnni
    class(serieMeseAnni)[class(serieMeseAnni)=="monthly"]<-"yearly"
    
    #ignore.par serve per utilizzare l'operazione "media" anche se si tratta di precipitazione
    #aggregaCD(serieMeseAnni,ignore.par=TRUE)
    aggregaCD(x=serieMeseAnni,ignore.par=TRUE,max.na=max.na,rle.check=rle.check,max.size.block.na=max.size.block.na,seasonal=FALSE)
    
  }) %>% reduce(rbind)->ris
  
  names(ris)<-names(x)
  ris
  
}#fine climatologiciMensili 


#checkClimatol
checkClimatol<-function(climatol){
  #verifichiamo l'indice di climatol
  stopifnot(is.ClimateData(climatol) && "yearly" %in% class(climatol))
  stopifnot(length(unique(years(index(climatol))))==1) #un unico anno
  stopifnot(length(unique(months(index(climatol))))==12) #12 mesi tutti distinti
}#fine checkClimatol




#####################################
#ANOMALIE:
#
# Calcolo anomalie mensili. Parto dalle serie mensili dei dati e sottraggo la media
# climatologica. La media climatologica può essere calcolata mediante i criteri del WMO
# o seguendo dei criteri più blandi.


#differenza tra ciascuna serie mensile e i valori climatologici mensili
OpWithClimatol<-function(fun){

  function(x,y){
    
    purrr::map2(.x=as.data.frame(coredata(x),optional=TRUE),.y=as.data.frame(coredata(y),optional=TRUE),.f=function(x,y){

      as.numeric(fun(x,y))
      
    })
    
  }#fine diffWithClimatol    
  
}#fine 

diffWithClimatol<-OpWithClimatol(fun=`-`)
sumWithClimatol<-OpWithClimatol(fun=`+`)

anomalie<-function(x,...)
{
  UseMethod("anomalie",x) 
}

anomalie.monthly<-function(x,climatol)
{
  stopifnot(is.ClimateData(x) && "monthly" %in% class(x))
  stopifnot((xtsAttributes(x)$parametro)==(xtsAttributes(climatol)$parametro))
  stopifnot(names(climatol)==names(x))  
  checkClimatol(climatol)
    
  xtsAttributes(x)$parametro->parametro

  #differenza tra ciascuna serie mensile e i valori climatologici mensili
  diffWithClimatol(x,climatol)->ris

  xts(as.data.frame(ris,optional=TRUE),order.by = index(x))->xris
  rm(ris)
  class(xris)<-c("ClimateData","monthly",class(xris))
  xtsAttributes(xris)<-list("parametro"=xtsAttributes(x)$parametro)

  xris
  
}#fine anomalie.monthly


#anomalie serie annuali, differenza tra una serie annuale e il suo climatologico
anomalie.yearly<-function(x,climatol)
{
  stopifnot(is.ClimateData(x) && "yearly" %in% class(x))
  stopifnot(is.ClimateData(climatol) && "yearly" %in% class(climatol))
  stopifnot((xtsAttributes(x)$parametro)==(xtsAttributes(climatol)$parametro))
  stopifnot(names(climatol)==names(x))
    
  xtsAttributes(x)$parametro->parametro
  #differenza tra ciascuna serie annuale e i valori climatologici mensili
  diffWithClimatol(x,climatol)->ris  

  xts(as.data.frame(ris,optional=TRUE),order.by = index(x))->xris
  rm(ris)
  class(xris)<-c("ClimateData","yearly",class(xris))
  xtsAttributes(xris)<-list("parametro"=xtsAttributes(x)$parametro)
  
  xris

}#fine anomalie.yearly



#RICOSTRUISCI SERIE: operazione inversa a quella effettuata da anomalie
#Prende un oggetto prodotto da anomalie e ricostruisce le serie
ricostruisciSerie<-function(x,...)
{

  UseMethod("ricostruisciSerie",x) 
}


ricostruisciSerie.monthly<-function(x,climatol)
{
  
  stopifnot(is.ClimateData(x) && "monthly" %in% class(x))
  stopifnot(names(climatol)==names(x))
  stopifnot((xtsAttributes(x)$parametro)==(xtsAttributes(climatol)$parametro))
  checkClimatol(climatol)
  
  #differenza tra ciascuna serie mensile e i valori climatologici mensili
  sumWithClimatol(x,climatol)->ris  

  xts(as.data.frame(ris,optional=TRUE),order.by = index(x))->xris
  rm(ris)

  ClimateData(x=xris,param=xtsAttributes(x)$parametro)
  
}#fine anomalie.monthly


ricostruisciSerie.yearly<-function(x,climatol)
{
  
  stopifnot(is.ClimateData(x) && "yearly" %in% class(x))
  stopifnot(is.ClimateData(climatol) && "yearly" %in% class(x))  
  stopifnot(names(climatol)==names(x))
  stopifnot((xtsAttributes(x)$parametro)==(xtsAttributes(climatol)$parametro))

  
  #differenza tra ciascuna serie mensile e i valori climatologici mensili
  sumWithClimatol(x,climatol)->ris
  
  xts(as.data.frame(ris,optional=TRUE),order.by = index(x))->xris
  rm(ris)
  
  ClimateData(x=xris,param=xtsAttributes(x)$parametro)
  
}#fine anomalie.yearly







