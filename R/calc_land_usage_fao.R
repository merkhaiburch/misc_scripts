


# Calculate the total amount of crop land used worldwide for
# corn, soy, wheat, rice, and sorghum
# data from FAOSTAT's data page: https://www.fao.org/faostat/en/#data/QCL

# Load in packages
library(dplyr)
library(ggplot2)

# Load in data
fao <- read.csv("~/Box Sync/Cornell_PhD/publications/Khaipho-Burch 2022 Yield Letter/FAOSTAT_data_en_8-25-2022.csv")

# subset out kust yield data
yield <- fao %>% filter(Element == "Yield")

# Calculate sum of land usage
land <- fao %>% 
  filter(Element == "Area harvested") 
sum(land$Value, na.rm = TRUE)

# Plot total yield across countries by crop
ggplot(yield, aes(x = Year, y = Value, group = Item))+
  geom_line(aes(color = Item)) 


p<-ggplot(df2, aes(x=dose, y=len, group=supp)) +
  geom_line(aes(color=supp))+
  geom_point(aes(color=supp))
p