library(dplyr)
library(ggplot2)

# Slack trends by person
messages_person <- read.csv("~/Downloads/Buckler Lab Member Analytics All time - Aug 16, 2023.csv") %>% 
  select("Name", "Account.created..UTC.", "Messages.posted", "Position")

colnames(messages_person) <- c("Name", "Date", "Count", "Role")

# Remove extra formatting around year
messages_person$Year <- gsub("-.*", "", messages_person$Date)
messages_person$Year <- as.numeric(messages_person$Year)

# Get days since joining
messages_person$daysSince <-  as.Date(messages_person$Date)
messages_person$daysSince <- as.numeric(Sys.Date() - messages_person$daysSince)

# Filter out low posting individials 
messages_person <- messages_person %>% 
  filter(Count > 250 & daysSince > 50 & Role != "Visiting Scientisist")

# Remove people with no assigned role
messages_person <- messages_person[-which(messages_person$Role == ""), ]

# Run model
summary(lm(Count ~ daysSince, data = messages_person))
summary(lm(Count ~ daysSince + Role, data = messages_person))

# Make a plot
daysSince <- ggplot(messages_person, aes(x = daysSince, y = Count, color = Role)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm")+
  facet_wrap(.~Role, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Days Since Joining Slack") +
  ylab("Total Number of Slack Messages Sent")
ggsave("~/Downloads/number_of_messages_days_since_joining.png", daysSince)


# ------------------------------------------------------------------------------
# Slack trends by person in the past year with reactions
# ------------------------------------------------------------------------------

# Slack trends by person
messages_2022 <- read.csv("~/Downloads/Buckler Lab Member Analytics Aug 1, 2022 - Aug 1, 2023.csv") %>% 
  select("Name", "Account.created..UTC.", "Days.active", "Messages.posted", "Reactions.added")
colnames(messages_2022) <- c("Name", "Date", "Days_active", "Count", "Reaction_Count")

# Subset last data for merging
messages_person_sub <- messages_person %>% select("Name", "Role", "daysSince")

# Merge with messages to get role information
messages_2022 <- merge(messages_2022, messages_person_sub, by = "Name")

# Filter out low posting individials 
messages_2022 <- messages_2022 %>% 
  filter(Count > 250)

# Run models
summary(lm(Count ~ Reaction_Count, data = messages_2022))
summary(lm(Count ~ Reaction_Count + Role, data = messages_2022))

# the number of reactions vs total messages
ggplot(messages_2022, aes(x = Count, y = Reaction_Count, color = Role)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm")+
  facet_wrap(.~Role) +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Total Count of Messages") +
  ylab("Total Number of Message Reactions (emojis)")
ggsave("~/Downloads/number_of_reactions_days_since_joining.png", )

# Days since joining slack and the number of reactions
ggplot(messages_2022, aes(x = daysSince, y = Reaction_Count, color = Role)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm")+
  facet_wrap(.~Role, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Days Since Joining Slack") +
  ylab("Total Number of Message Reactions (emojis)")
ggsave("~/Downloads/number_of_reactions_days_since_joining.png", )

# Total active days and number of messages
ggplot(messages_2022, aes(x = Days_active, y = Count, color = Role)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm")+
  facet_wrap(.~Role, scales = "free") +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Days Active On Slack") +
  ylab("Total Number of Messages")
ggsave("~/Downloads/number_of_reactions_days_since_joining.png", )


# histograms -------------------------------------------------------------------

# Get top reactors
messages_2022 %>% 
  select(Name, Reaction_Count) %>% 
  arrange(desc(Reaction_Count)) %>% 
  slice_max(Reaction_Count, n = 5)

# Get least reactors
messages_2022 %>% 
  select(Name, Reaction_Count) %>% 
  arrange(desc(Reaction_Count)) %>% 
  slice_min(Reaction_Count, n = 15)

library(ggridges)
messages_2022_sub <- messages_2022 %>% filter(Role != "Undergrad")

v10 <- ggplot(messages_2022_sub, aes(x = Reaction_Count, y = Role, fill = Role, color = Role)) + 
  geom_density_ridges(scale = 1) +
  theme(legend.position = "none")
ggsave("~/Downloads/number_of_reactions_aug2022_aug2023.png", v10)
