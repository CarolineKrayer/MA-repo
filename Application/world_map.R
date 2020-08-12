########################################################
# Create worldmap with group membership of each country.
########################################################
library(sf)
library(tmap)
library(spData)

rm(list=ls())
out_path = "./Application/Results/"

# Load list with group membership.
groups = list.load(file=paste0(out_path, "Group_membership_all_countries", all_countries, ".rdata"))
group_1 = groups[[1]]
group_2 = groups[[2]]

# Alter countries whose name did not coincide.
#`%!in%` = Negate(`%in%`)
#world_group1 = world[world$name_long %in% group_1, ]
#world_group2 = world[world$name_long %in% group_2, ]
#group_1[group_1 %!in% world_group1$name_long]
#group_2[group_2 %!in% world_group2$name_long]
group_1 = replace(group_1, group_1=="Central African Rep.", "Central African Republic")
group_1 = replace(group_1, group_1=="Congo, Republic of", "Republic of the Congo")
group_1 = replace(group_1, group_1=="Congo, Dem. Rep. of", "Democratic Republic of the Congo")
group_1 = replace(group_1, group_1=="Korea", "Republic of Korea")
group_1 = replace(group_1, group_1=="Gambia, The", "The Gambia")
group_1 = replace(group_1, group_1=="Lao People's Dem.Rep", "Lao PDR")
group_1 = replace(group_1, group_1=="Syrian Arab Republic", "Syria")
group_1 = replace(group_1, group_1=="Swaziland", "eSwatini")

group_2 = replace(group_2, group_2=="C<f4>te d'Ivoire", "C?te d'Ivoire")
group_2 = replace(group_2, group_2=="Venezuela, Rep. Bol.", "Venezuela")
group_2 = replace(group_2, group_2=="Slovak Republic", "Slovakia")
group_2 = replace(group_2, group_2=="Iran, I.R. of", "Iran")
group_2 = replace(group_2, group_2=="Kyrgyz Republic", "Kyrgyzstan")

# Add column for group membership to world data frame of SpData package.
world["group"] = 0
world$group[world$name_long %in% group_1] = 1
world$group[world$name_long %in% group_2] = 2

# Map countries with groups coloured.
worldmap = tm_shape(world) + 
           tm_fill(style="cat", 
                   col="group", 
                   title="Classification", 
                   labels=c("Not included", "Group 1", "Group 2"),
                   palette = c("#eaeab4", "#f98c1e", "#CC4D00")) +
           tm_legend(position=c(0.01, 0.2),
                     frame=TRUE,
                     title.size=1.4,
                     text.size=1.1) +
           tm_borders()

# Save worldmap.
tmap_save(tm=worldmap, filename=paste0(out_path, "Worldmap.png"))

