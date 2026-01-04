# Load required libraries
library(dplyr)
library(ggplot2)
setwd("~/Project_data/mittag_group/temperature_Transcriptome/Cell_tracking/2025_raw_tracks/")
spots <- read.csv("result_Files/MotGraph_Spot/20250519_Clamy_28_P1_spots.csv")
colnames(spots)


spots <- read.csv("result_Files/MotGraph_Spot/20250610_R1_Clamy_18_spots_P1.csv")
colnames(spots)

###IDS for 28
##3106, 2324, 2325

##IDS for 18
#12111, 9382
# If your file has multiple header rows, skip them by adjusting 'skip' parameter:
# spots <- read.table("trial.spots.txt", header = TRUE, sep = "\t", comment.char = "#", skip = 3, stringsAsFactors = FALSE)

# Ensure the key columns are present and named correctly
# If not, check the column names:
# print(colnames(spots))

# Convert columns to numeric if needed
spots$TRACK_ID <- as.numeric(spots$TRACK_ID)
spots$POSITION_X <- as.numeric(spots$POSITION_X)
spots$POSITION_Y <- as.numeric(spots$POSITION_Y)

# Remove any rows with missing values in key columns
spots <- spots %>%
  filter(!is.na(TRACK_ID), !is.na(POSITION_X), !is.na(POSITION_Y))


table(spots$TRACK_ID)
desired_id <- c(12111)  # Change this to the TRACK_ID you want

single_track <- spots %>% 
  filter(TRACK_ID %in% desired_id) %>%
  arrange(TRACK_ID, POSITION_T) %>%
  group_by(TRACK_ID) %>%
  mutate(start_time = first(POSITION_T)) %>% 
  filter(POSITION_T <= start_time + 60) %>% #Sets time
  # Normalize positions
  mutate(
    X0 = POSITION_X - first(POSITION_X),
    Y0 = POSITION_Y - first(POSITION_Y)
  ) %>%
  ungroup()

# Normalize each track to start at (0,0)
spots_norm <- single_track %>%
  group_by(TRACK_ID) %>%
  mutate(
    X0 = POSITION_X - dplyr::first(POSITION_X),
    Y0 = POSITION_Y - dplyr::first(POSITION_Y)
  )

#write.table(spots_norm, "28_normalizedSpots.txt",
#            sep="\t", quote = F, row.names = F)

# For each TRACK_ID, get the last point (by time)
last_points <- spots_norm %>%
  arrange(TRACK_ID, POSITION_T) %>%
  group_by(TRACK_ID) %>%
  slice_tail(n = 1)

##18
safe_colorblind_palette <- c(
  "#258CCC")

##28
safe_colorblind_palette <- c(
  "#428033"
)

# Use in ggplot2:
scale_color_manual(values = safe_colorblind_palette)


ggplot(spots_norm, aes(x = X0, y = Y0, group = TRACK_ID, color = as.factor(TRACK_ID))) +
  geom_path(size = 1) +
  scale_color_manual(values = safe_colorblind_palette) +
  #geom_point(data = last_points, aes(x = X0, y = Y0, fill = as.factor(TRACK_ID)),
  #           size = 3, shape =21 , stroke = 2) +
  annotate("point", x = 0, y = 0, color = "black", size = 9, shape = 16) +  # Mark the origin
  labs(
    title = "Selected Tracks (Normalized to Origin)",
    x = "X (relative to start)",
    y = "Y (relative to start)",
    color = "Track ID",
    fill = "Track ID"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 25),
    axis.title.y = element_text(face = "bold", size = 25),
    axis.text.x = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 22),
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    panel.grid.major = element_line(color = "grey", size = 0.1), 
    axis.line.x = element_line(size = 1.5, color = "black"),
    axis.line.y = element_line(size = 1.5, color = "black"),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 25),
    strip.text.x = element_text(size=26, color="black"),
    strip.background = element_rect(fill="grey")
  )


####Unused, axis through the origin.
ggplot(spots_norm, aes(x = X0, y = Y0, group = TRACK_ID, color = as.factor(TRACK_ID))) +
  geom_path(size = 1) +
  scale_color_manual(values = safe_colorblind_palette) +
  geom_hline(yintercept = 0, color = "black", size = 1.2) +   # Add horizontal axis at y=0
  geom_vline(xintercept = 0, color = "black", size = 1.2) +   # Add vertical axis at x=0
  #geom_point(data = last_points, aes(x = X0, y = Y0, fill = as.factor(TRACK_ID)),
  #           size = 3, shape =21 , stroke = 2) +
  annotate("point", x = 0, y = 0, color = "black", size = 6, shape = 16) +  # Mark the origin
  labs(
    title = "Selected Tracks (Normalized to Origin)",
    x = "X (relative to start)",
    y = "Y (relative to start)",
    color = "Track ID",
    fill = "Track ID"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 25),
    axis.text.x = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 22),
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
    panel.grid.major = element_line(color = "grey", size = 0.1), 
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 25),
    strip.text.x = element_text(size=26, color="black"),
    strip.background = element_rect(fill="grey")
  )


