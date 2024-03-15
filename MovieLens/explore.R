ratings <- read_csv("MovieLens/ml-latest/ratings.csv",show_col_types = FALSE)
# links <- read_csv("MovieLens/ml-latest/links.csv",show_col_types = FALSE)
movies <- read_csv("MovieLens/ml-latest/movies.csv",show_col_types = FALSE)
genome_scores <- read_csv("MovieLens/ml-latest/genome-scores.csv",show_col_types = FALSE)
genome_tags <- read_csv("MovieLens/ml-latest/genome-tags.csv",show_col_types = FALSE)
tags <- read_csv("MovieLens/ml-latest/tags.csv",show_col_types = FALSE)


head(tags)
head(genome_tags)
head(genome_scores)
head(movies)
head(ratings)
