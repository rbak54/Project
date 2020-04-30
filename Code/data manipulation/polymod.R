library(socialmixr)
data("polymod")
data <- polymod$contacts
contact_table <- table(data$part_id, data$duration_multi)
del <- as.numeric(names(which(contact_table[, 1] > 0)))
data <- data[which(data$part_id != del), ]
contact_table <- table(data$part_id, data$duration_multi)
con <- as.matrix(contact_table)
midpoints <- c(2.5, 10, 15 + (60 - 45) / 2, 150, 4 * 60)

time <- midpoints * con / 60
timebyrow <- apply(time, 1, sum)
#timebyrow
mean(timebyrow)
sqrt(var(timebyrow))
max(timebyrow)
min(timebyrow)

data <- polymod$contacts
contact_table <- table(data$part_id, data$duration_multi)
del <- as.numeric(names(which(contact_table[, 1] > 0)))
data <- data[which(data$part_id != del), ]
over18s <-
  polymod$participants$part_id[which(polymod$participants$part_age > 18)]
data <- data[data$part_id %in% over18s,]
#data
contact_table <- table(data$part_id, data$duration_multi)
con <- as.matrix(contact_table)
midpoints <- c(2.5, 10, 15 + (60 - 45) / 2, 150, 4 * 60)
time
time <- midpoints * con / 60
timebyrow <- apply(time, 1, sum)
#timebyrow
mean(timebyrow)
sqrt(var(timebyrow))
min(timebyrow)
max(timebyrow)

