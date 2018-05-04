---
title: "Shiny tests"
date: "May 3, 2018"
output: html_document
---

View(linnarsson)

library(png)
library(rasterImage)

#edited dendrogram from http://mousebrain.org/dokuwiki/doku.php?id=clusters
png <- readPNG('./dendrogram-01.png')

#everything is 9px, starting at 0
