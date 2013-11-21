is.filtered <- function(df = filter.stats) {
  filter.scores <- df[ , grep("^filter", colnames(df))]
  df$is.filtered <-
    if(is.null(dim(filter.scores))) filter.scores else apply(filter.scores, 1, sum)
  df
}

filter1 <- function(df = filter.stats, min.ratio = 2) {
  filter <-
    ((df$ng.lt.ratio >= min.ratio) & (df$g.lt.ratio <= min.ratio)) |
    ((df$ng.rt.ratio >= min.ratio) & (df$g.rt.ratio <= min.ratio))
  df$filter1 <- 0
  df$filter1[ filter ] <- 1
  is.filtered(df)
}

filter2 <- function(df = filter.stats, min.ratio = 1.5) {
  filter <-
    (df$lt.ratio >= min.ratio) |
    (df$rt.ratio >= min.ratio)
  df$filter2 <- 0
  df$filter2[ filter ] <- 1
  is.filtered(df)
}

filter3 <- function(df = filter.stats, min.ratio = 1.5) {
  filter <-
    ((df$lt.ratio >= min.ratio) & (df$rt.ratio <= min.ratio)) |
    ((df$rt.ratio >= min.ratio) & (df$lt.ratio <= min.ratio))
  df$filter3 <- 0
  df$filter3[ filter ] <- 1
  is.filtered(df)
}

filter4 <- function(df = filter.stats, ratio.1 = 1.5, ratio.2 = 0.75) {
  filter <-
    ((df$lt.ratio >= ratio.1) & (df$rt.ratio <= ratio.2)) |
    ((df$rt.ratio >= ratio.1) & (df$lt.ratio <= ratio.2))
  df$filter4 <- 0
  df$filter4[ filter ] <- 1
  is.filtered(df)
}

prep.data <- function() {
  header <- c("chr", "pos", "ref", "a.ct", "c.ct", "g.ct", "t.ct", "del.ct",
              "alt", "ng.lt", "ng.snp", "ng.rt", "g.lt", "g.snp", "g.rt",
              "ng.lt.ratio", "g.lt.ratio", "ng.rt.ratio", "g.rt.ratio",
              "lt.ratio", "rt.ratio")
  yes.filtered <- read.table("yes-filtered", col.names = header, sep = "\t", as.is = T)
  not.filtered <- read.table("not-filtered", col.names = header, sep = "\t", as.is = T)

  yes.filtered$should.filter <- "YES"
  not.filtered$should.filter <- "NO"

  filter.stats <- rbind( yes.filtered, not.filtered )
  filter.stats <- filter.stats[, c(1:3, 9, 16:22)]
  filter.stats$is.filtered <- 0

  filter.stats
}

filter.stats <- prep.data()
filter.stats <- filter1()
# filter.stats <- filter2()
# filter.stats <- filter3()
filter.stats <- filter4(ratio.1 = 1.2, ratio.2 = 0.8)
filter.stats



> filter.stats$rt.lt <- apply(cbind(filter.stats$rt.ratio, filter.stats$lt.inv), 1, mean)
> filter.stats$lt.rt <- apply(cbind(filter.stats$lt.ratio, filter.stats$rt.inv), 1, mean)
> filter.stats$max <- apply(cbind(filter.stats$rt.lt, filter.stats$lt.rt), 1, max)


filter.stats$max.rat <- apply(cbind(filter.stats$lt.rat, filter.stats$rt.rat), 1, max)
filter.stats$min.rat <- apply(cbind(filter.stats$lt.rat, filter.stats$rt.rat), 1, min)




###########
# TESTING #
###########

header <- c("chr", "pos", "ref", "a.ct", "c.ct", "g.ct", "t.ct", "del.ct",
            "alt", "ng.lt", "ng.snp", "ng.rt", "g.lt", "g.snp", "g.rt")

df <- read.table("snps/IMB211.A01.snps.nogap.gap.csv", col.names = header,
                sep = ",", as.is = T, header = T)

df$ng.lt.ratio <- df$ng.lt / df$ng.snp
df$g.lt.ratio <- df$g.lt / df$g.snp
df$ng.rt.ratio <- df$ng.rt / df$ng.snp
df$g.rt.ratio <- df$g.rt / df$g.snp
df$lt.ratio <- df$ng.lt.ratio / df$g.lt.ratio
df$rt.ratio <- df$ng.rt.ratio / df$g.rt.ratio

head(df[ ((df$lt.ratio >= ratio.1) & (df$rt.ratio <= ratio.2)) | ((df$rt.ratio >= ratio.1) & (df$lt.ratio <= ratio.2)), ], 200)

head(df[ !( ((df$lt.ratio >= ratio.1) & (df$rt.ratio <= ratio.2)) | ((df$rt.ratio >= ratio.1) & (df$lt.ratio <= ratio.2)) ), ], 200)

df.NOTfiltered <- df[ !( ((df$lt.ratio >= ratio.1) & (df$rt.ratio <= ratio.2)) | ((df$rt.ratio >= ratio.1) & (df$lt.ratio <= ratio.2)) ), c(1:2, 20:21)]

head(df.NOTfiltered[with(df.NOTfiltered, order(lt.ratio, -rt.ratio)), ], 20)


df.filtered <- df[ ((df$lt.ratio >= ratio.1) & (df$rt.ratio <= ratio.2)) | ((df$rt.ratio >= ratio.1) & (df$lt.ratio <= ratio.2)), c(1:2, 20:21)]

head(df.filtered[with(df.filtered, order(lt.ratio, -rt.ratio)), ], 20)


########
# JUNK #
########

filter.stats <- filter1()

filter.stats$filter1 <- 0
filter.stats$filter1[ filter1() ] <- 1

filter.stats$is.filtered <- is.filtered()


length(filter1())




filter.stats[filter.stats$alt == "del", ]
