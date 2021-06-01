
# Steps to perform all cuts and export final data. "Stnd" is a placeholder for the name of the standard being called.

# Find initial cutpoints
data_cuts <- findCutpoints(data$D48CDES_Final)

# Perform first cut
data_firstcut <- data[data$D48CDES_Final >= data_cuts[1] & data$D48CDES_Final <= data_cuts[2] ]

# Get summary stats following first cut
mean(data_firstcut)
sd(data_firstcut)
length(data_firstcut)

tiff(file = "Data_firstcut.tiff", width = 5, height = 5, units = "in", res = 800, compression = "lzw")
plot(density(data_firstcut), main = expression(paste("Standard name (",Delta[48]*" \u2030), Config number")))
abline(v=mean(data_firstcut)+sd(data_firstcut), col = "#FCA50AFF")
abline(v=mean(data_firstcut)-sd(data_firstcut), col = "#FCA50AFF")
abline(v=mean(data_firstcut)+(2*sd(data_firstcut)), col = "#DD513AFF")
abline(v=mean(data_firstcut)-(2*sd(data_firstcut)), col = "#DD513AFF")
abline(v=mean(data_firstcut)+(3*sd(data_firstcut)), col = "#6B186EFF")
abline(v=mean(data_firstcut)-(3*sd(data_firstcut)), col = "#6B186EFF")
abline(v=mean(data_firstcut), lty=2, col = "#170C3AFF")
dev.off()

# Perform 3sigma exclusion as the second cut. If the Shapiro-Wilk test (line 31) indicates non-normality, use 2sigma (2*sd) or 1sigma (1*sd) as needed.
Stnd_secondcut <- Stnd_firstcut[Stnd_firstcut >= (mean(Stnd_firstcut)-(3*sd(Stnd_firstcut))) & 
                                Stnd_firstcut <= (mean(Stnd_firstcut)+(3*sd(Stnd_firstcut)))]

# Perform a Shapiro-Wilk test for normality following the 3 sigma exclusion
shapiro.test(Stnd_secondcut)

# Get summary stats following second cut
mean(Stnd_secondcut)
sd(Stnd_secondcut)
length(Stnd_secondcut)

# Quick plot to visualize the final data
plot(density(Stnd_secondcut), main = "Stnd Final Data")

# Apply the calculated cuts to the full dataset
Stnd_final <- data[data$Standard == "Stnd" & data$D48CDES_Final >= range(Stnd_secondcut)[1] & data$D48CDES_Final <= range(Stnd_secondcut)[2],]

# Create an Excel spreadsheet for the data
wb <- createWorkbook("FinalData") # Create an empty workbook in the Global Environment
addWorksheet(wb, "Stnd") # Add a blank sheet for the standard being cleaned
writeData(wb, sheet = "Stnd", Stnd_final) # Write the final, cleaned data to the blank sheet

# Repeat above for each standard, then:

# Final step is to save the complete workbook of data
saveWorkbook(wb, "Final Data.xlsx")

# We suggest tracking basic summary stats, cuts used, final sigma, etc. for each standard following the example provided in All_standards_exclusions_config1.xlsx