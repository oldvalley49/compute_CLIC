library(ggplot2)

mouse <- read.csv("output/scores/mouse_scores.csv")
human <- read.csv("output/scores/human_scores.csv")

mouse$source <- "mouse"
human$source <- "human"

combined <- rbind(human, mouse)
p <- ggplot(combined, aes(x=source, y = pearson_correlation, fill=source)) +
        geom_violin() +
        labs(title = "Scores Across Species", x="", y="Pearson's Correlation Coefficient")
ggsave("/dcs07/hongkai/data/tomo/ENCODE_score/plots/score_vln.jpeg")