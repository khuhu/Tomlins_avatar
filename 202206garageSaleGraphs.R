library(ggplot2)

garageSale <- readxl::read_xlsx("/mnt/DATA5/tmp/kev/misc/20220604garageSale.xlsx",
                                col_names = FALSE)
colnames(garageSale) <- c("Item", "Catergory", "Person", "Amount")
garageSale <- garageSale[-which(garageSale$Person == "Russell"), ]
garageSale <- garageSale[-which(garageSale$Amount < 1),]
garageSale$Amount <- round(garageSale$Amount)

garageSale$Catergory <- str_replace_all(garageSale$Catergory, "Books", "Book")

ggplot(garageSale, aes(x = Person, fill = factor(Catergory))) +
  geom_bar(stat = "count")

### reformat for dollar in category count

finalCountTable <- data.frame()
for (i in 1:nrow(garageSale)) {
  tmp2 <- NULL
  tmp <- c(garageSale$Catergory[i], garageSale$Person[i])
  for (j in 1:(garageSale$Amount[i])) {
   tmp2 <- rbind(tmp2, tmp) 
  }
  
  finalCountTable <- rbind(finalCountTable, tmp2)
}

colnames(finalCountTable) <- c("categories", "person")
my_colors <- c("blue4", "darkgoldenrod4", "forestgreen",
               "indianred3", "mediumslateblue")

names(my_colors) <- levels(finalCountTable$categories)
colScale <- scale_colour_manual(name = "grp",values = my_colors)

ggplot(finalCountTable, aes(x = person, fill = categories)) +
  geom_bar(stat = "count") + 
  theme_bw() + colScale + xlab("Person") + ylab("Total dollars") +
  ggtitle("Amount everyone made by item category") + 
  theme(plot.title = element_text(hjust = 0.5))



