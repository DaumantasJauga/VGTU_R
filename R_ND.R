###############################################################1.1 Nepadaryta

###############################################################1.2  JSON failu nuskaitymas
library(rjson)
library(jsonlite)
library(curl)
library(rlist)

#-----------------------1--------------------------
getwd()
setwd("C://Users/dauma/OneDrive/Documents")

#A
url <- "http://jeroen.github.io/data/nycflights13.json.gz"
con <- gzcon(curl(url))
rez <- stream_in(con)
toJSON(rez)

#B
i <- 1 
i
fun <- function(i){
  
  visas <- gzfile("nycflights13.json.gz", "r")
  js <- readLines(visas, i)
  list(js)
}

fun(1)

#C
vard <- function(vard){
  
  duom <- fromJSON(sprintf("[%s]", paste(readLines(con),collapse=",")))[vard]
  
  return (unlist(duom))
}

(readNames("year"))
  
#d
con <- gzfile("nycflights13.json.gz", "r")
lines <- readLines(con)
close(con)
df1 <- jsonlite::stream_in(textConnection(lines), verbose = FALSE)
df <- data.frame(df1)
splitData <- split(df, df$month) 
month1 <- toJSON(splitData[1]) 
month2 <- toJSON(splitData[2]) 
  
write_json(splitData)

#e
df <- data.frame(df1)
splitData <- split(df, df$carrier)
splitData
write.csv(splitData$`9E`, "C:/Users/dauma/Downloads/myfile.csv")
lapply(names(splitData), function(df) write.csv(splitData[[df]], file=paste0("C:/Users/dauma/Downloads/"), df, ".csv"))
  
#f
con <- gzfile("nycflights13.json.gz", "r")
lines <- readLines(con)
close(con)
df1 <- jsonlite::stream_in(textConnection(lines), verbose = FALSE)
df <- data.frame(df1)
df
  
#----------------Kita dalis
  
#A
unq <- unique(df$carrier)
unq
length(unq)
  
#B
tbl <- table(df$origin,df$dest)
tb <- data.frame(tbl)
subset(tb,tb$Freq == max(tb$Freq))
  
#c
tbl <- table(df$tailnum)
tb <- data.frame(tbl)
tb1 <- subset(tb, tb$Var1 != '')
tb2 <- subset(tb1,tb1$Freq == max(tb1$Freq))
tb2
  
#d
df1 <- sqldf('with main as(select flight, max(distance) as distance from df group by flight)
             select flight, max(distance) from main')
df1
max(tb$tp)
  
  
#e
df
tp <- data.frame(with(df, tapply(tailnum,month, length)))
tp
tp1 <-subset(tp,tp$with.df..tapply.tailnum..month..length.. == max(tp$with.df..tapply.tailnum..month..length..))
tp1
names(tp1) <- "len"
tp1
df
  
#f
i <- unique(df$tailnum)
i
#
  
df1 <- sqldf('with main as(select month, count(distinct(tailnum)) as kiekis from df group by month)
             select * from main')
  
df1

#g
df
df1 <-  sqldf('with main as (select count(distinct(tailnum)) as kiekis from df
               where origin = "JFK" and dest = "LAX" 
               group by month)
               
              select * from main')
  
df1
  
#h
tb <- tapply(df$distance ,df$carrier, mean)
tb
  
#i
df1 <-  sqldf('with main as (select flight, count(flight) as kiekis, arr_delay from df
               group by flight
               order by count(flight) desc)
               
              select flight,kiekis,avg(arr_delay) from main 
              group by flight
              order by kiekis desc limit 3')
df1
  
#j
df
df1 <-  sqldf('with main as (select flight, hour*60 + minute as laikas from df
               group by flight
               )
              
              select flight, avg(laikas) from main
              group by flight')
  

###############################################################1.3 OpenWeatherMap meteorologiniu duomenu failu nuskaitymas

###################################1
# a
library(jsonlite)
library(curl)
library("sqldf")
con <- gzcon(gzfile("daily_14.json.gz"))
rez <- fromJSON(paste("[", paste(readLines(con), collapse = ","), "]"))

df1 <- sqldf("select *  from rez where city.id = 3632308")


df1
#con <- gzcon(curl(url))
#rez <- stream_in(con)
#head(rez)

  
#b
con <- gzcon(gzfile("daily_14.json.gz"))
rez <- fromJSON(paste("[", paste(readLines(con), collapse = ","), "]"))
attributes(rez)
class(rez)
df1 <- sqldf('with main as(select flight, max(distance) as distance from rez group by flight)
             select flight, max(distance) from main')

#c
df1 <- sqldf('with main as(select max(city.coord.lat), city.name, city.country  from rez),
             
             main2 as(select min(city.coord.lat), city.name, city.country  from rez)
             
            select * from main
            
            union all
            
            select * from main2
         
             ') 

#c
df1 <- sqldf("select city.name, city.coord.lat, city.coord.lon from rez where city.country = 'LT'")

#################################2
#a
con <- gzcon(gzfile("daily_16.json.gz"))
rez <- fromJSON(paste("[", paste(readLines(con), collapse = ","), "]"))

fromJSON(con)

df1 <- sqldf("select count(distinct(city.id))  from rez where city.country = 'LT'")


#B

df1 <- sqldf("select city.name, city.coord.lon, city.coord.lat  from rez where city.country = 'LT'")

#C

df1 <- sqldf("with main as (select city.country, count(city.id) as kiek  from rez group by city.country)
             
             select * from main order by kiek desc")



###############################################################2.1 Tekstiniu duomenu irašymas i HDF5 formato faila

#a
library(rhdf5)
library(microseq)


duom <- readFasta("AE016828.ffn")
duom

i <- 1;
while(i <= nrow(duom)) {
  
  fileName <- paste0("AE016828_", i, ".txt")
  file.create(fileName)
  write(as.vector(unlist(duom[i,])), fileName)
  i <- i + 1

}

#b

h5createFile("AE016828.hdf5")
duom <- readFasta("AE016828.ffn")
duom

h5write(x[,2], "AE016828.hdf5", name = "nukleotidai")
h5read("AE016828.hdf5", "nukleotidai")



############################################################### 3.1 Duomenu lenteles iš XML failo sudarymas

library(XML)
#--------------------------------1
getwd()
doc <- htmlParse(file = 'bond.xml')
doc


doc_name <- getNodeSet(doc, path = '//name')
name <- xmlSApply(doc_name,xmlValue)
name

doc_year <- getNodeSet(doc, path = '//year')
year <- xmlSApply(doc_year,xmlValue)
year

doc_bud <- getNodeSet(doc, path = '//budget')
budget <- xmlSApply(doc_bud,xmlValue)

actors <- c('','','')

doc_bo <- getNodeSet(doc, path = '//boxoffice')
boxoffice <- xmlSApply(doc_bo,xmlValue)


data.frame(name,year,actors,budget,boxoffice)
xmlToDataFrame('bond.xml')

################################################################ 3.2 Elementu iš XML failo išrinkimas

#1
download_xml("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=20074371&retmode=xml", "NDxml.xml")
doc <- xmlParse(file = 'NDxml.xml')
#A
data <- getNodeSet(doc, path = '//ArticleTitle')
xmlSApply(data,xmlValue)


#B
data <- getNodeSet(doc, path = "//ArticleIdList/ArticleId[@IdType = 'doi']")
xmlSApply(data,xmlValue)

#C
data <- getNodeSet(doc, path = "//History/PubMedPubDate[@PubStatus = 'pubmed']/child::node()")
temp <- xmlSApply(data,xmlValue)
paste(temp,collapse = ' ')

#D
data <- getNodeSet(doc, path = "//Abstract/AbstractText [@Label = 'RESULTS']")
temp <- xmlSApply(data,xmlValue)
summary(temp)

#E
data <- getNodeSet(doc, path = "//MeshHeadingList/MeshHeading/DescriptorName/text()")
data

#F
data <- getNodeSet(doc, path = "//MeshHeadingList/MeshHeading/DescriptorName[@UI = 'D016014']/text()")
data

#G
data <- getNodeSet(doc, path = "//ReferenceList/Reference/Citation/ArticleIdListtext()")
data


#H
data <- getNodeSet(doc, path = "//ReferenceList/Reference/Citation/text()")
ln <- length(data)
ln


################################################################ 4.1 Duomenu lenteliu nuskaitymas
  
library(sqldf)
library(data.table)
library(plyr)
library(rvest)
library(xml2)
library(dplyr)

################1
#A
url <- "https://www.scimagojr.com/countryrank.php"
doc <- read_html(url)

lentele <- html_table(doc, fill = TRUE)
lentele

#B
url <- "https://www.scimagojr.com/journalrank.php"
doc <- read_html(url)
?read_html
lentele <- html_table(doc, fill = TRUE)
DF <- data.frame(lentele)
ptn = 'Math'
DF[grep(ptn, DF$Title, perl=T),]


################2
#A
url <- "https://en.wikipedia.org/wiki/List_of_World_Heritage_in_Danger"
doc <- read_html(url)
lentele <- html_table(doc, fill = TRUE)
lentele[[2]]

#2B
ptn <- "Syria|Libya"
ptn
lent <- lentele[[2]][grep(ptn, lent$Location, perl=T),]
lent[,c("Name","Location")]

################################################################ 4.2 Duomenu iš HTML puslapio nuskaitymas
#a
fun <- function(i){
  ur <- "http://www.piliakalniai.lt/piliakalnis.php?piliakalnis_id="
  url <- paste(ur,i,collapse = NULL,sep = '')
  doc <- read_html(url)
  lentele <- html_table(doc, fill = TRUE)
  return(lentele)
}
fun(2)

#b
df <- data.frame(matrix(ncol = 12, nrow = 0))
i<-1
i
while(TRUE){
  ur <- "http://www.piliakalniai.lt/piliakalnis.php?piliakalnis_id="
  url <- paste(ur,i,collapse = NULL,sep = '')
  doc <- read_html(url)
  lentele <- data.frame(html_table(doc, fill = TRUE))
  df <- rbind(df, t(lentele[,2]))
  i = i+1
  #Padaryta iki 20, kad nereiktų ilgai laukti dėl visų 869
  if(i > 20 ) break()
}
x <- c('Pavadinimas','Kiti pavadinimai','Apskritis','Rajonas','Seniunija','LKS','WGS','Registro nr.','MC','Kiti MC','Apsaugos tikslas','Rušis')
colnames(df) <- x 
df

#c
df_frq <- data.frame(tapply(df$Pavadinimas,df$Rušis,length))
colnames(df_frq) <- "frequency"
df_frq

#d
sav <- tapply(df$Pavadinimas,df$Rajonas,length)
sav <- data.frame(sav)
colnames(sav) <- "skaičius"
sav

aps <- tapply(df$Pavadinimas,df$Apskritis,length)
aps <- data.frame(aps)
colnames(aps) <- "skaičius"
aps


################################################################ 5.1 Objektu su nurodytomis koordinatemis atvaizdavimas žemelapyje
  
library(maps)
library(sqldf)
library(data.table)
library(plyr)
library(rvest)
library(xml2)
library(ggmap)
library(dplyr)



#a
url <- "http://www.qrz.lt/lhfa/index.php/en/lhf-list-s-3/97-visi-piliakalniai"
doc <- read_html(url)

lentele <- html_table(doc, fill = TRUE)
lentele <- lentele[[2]]
lentele
lentele <- sqldf('select x4 from lentele where length(x4) = 29')
lentele
lentele <- lentele$X4


#b
pir = substr(lentele, 2, 3)
pir
antr = substr(lentele, 6, 7)
antr
trec = substr(lentele, 10, 11)
trec

long = paste0(pir,".",antr,trec)
long

latpir = substr(lentele, 17, 18)
latpir
latantr = substr(lentele, 21, 22)
latantr
lattrec = substr(lentele, 25, 26)
lattrec

lat = paste0(latpir,".",latantr,lattrec)
mode(lat)


df <-data.frame(long,lat)
df

map("world", regions = "Lithuania")
plot(žemėlapis, type = "l", ann = FALSE, frame = FALSE)
grid()
points(lat,long, pch = 1, cex = 1)
?points

#c
visas <- sqldf('with main as (select  lat, long, row_number() over (order by long desc) as row_num from df),
                     main2 as (select  lat, long, row_number() over (order by lat desc) as row_num from df)
                     
                     select * from main where row_num in ( 1, 908)
                     
               union 
              
                     select * from main2 where row_num in ( 1, 908) ')
visas
visas_long <- as.vector(visas$long)
visas_long
visas_lat <- as.vector(visas$lat)
class(visas)
plot(žemėlapis, type = "l", ann = FALSE, frame = FALSE)
points(visas_lat,visas_long, pch = 1, cex = 1)







