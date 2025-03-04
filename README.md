# extractOBISdata
*Author: Narimane Dorey : narimane.dorey_at_gmail.com from an initial idea from C. Clerc*

## R tool to extract data from the OBIS database (OBIS.org) and plot data on a map

This script will retrieve data for 5 species of cold-water corals from the OBIS database (*Lophelia pertusa*, *Madrepora oculata*, *Solenosmilia variabilis*, *Bathelia candida*, *Enallopsammia profunda*) and plot their occurence and average depth of discovery on a world map. 

It uses the packages [robis](https://github.com/iobis/robis), ggplot2 and leaflet.

> [OBIS](https://obis.org) **: Ocean Biodiversity Information System** is a global open-access data and information clearing-house on marine biodiversity for science, conservation and sustainable development.



### Map of the occurence of Desmophyllum (_Lophelia pertusa_) registered in OBIS, with ggplot :

![Desmophyllum map](https://github.com/DrUrchin/extractOBISdata/blob/main/Desmophyllumpoints_map.png)

### Map of the occurence of Madrepora (_M. oculata_),registered in OBIS, with leaflet :
Screenshot of the interactive html map produced.

![Madrepora map](https://github.com/DrUrchin/extractOBISdata/blob/main/Madrepora_leaflet_interactivemap.png)
