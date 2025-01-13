# The Demographic Transition in the United States ----
# Chiara Cattani


# Prepare environment
rm(list = ls())
graphics.off()
setwd(file_path <- dirname(rstudioapi::getSourceEditorContext()$path))

# Packages
library(tidyverse)
library(ggplot2)
library(demography)
library(pyramid)
library(dplyr)
library(mFilter)
library(StMoMo)
library(RColorBrewer)


# Data ----

# http://www.mortality.org/
username <- "chiara.cattani9@studio.unibo.it"
password <- "mi7i!sX69MDSDzC"


## Mortality rates and population ----

us <- read.demogdata("data/us.mx.txt", "data/us.exp", type = "mortality", label = "USA", skip = 2)
#us <- hmd.mx("USA", username, password)


## Life expectancy ----

us0 <- as.data.frame(hmd.e0("USA", username, password))
year <- seq(1933, 2021, 1)
us_ok <- cbind(us0, year)

us_female <- us_ok %>%
  select(Female, year) %>%
  mutate(sex = "Female", Value = Female) %>%
  select(-Female)

us_male <- us_ok %>%
  select(Male, year) %>%
  mutate(sex = "Male", Value = Male) %>%
  select(-Male)

us_ok_long <- rbind(us_male, us_female)


## Total fertility rate ----

us.tfr <- read.table("data/us.tfr.txt", header = TRUE, skip = 2)


## GDP ----

us.gdp <- read.csv("data/us.gdp.csv")

us.gdp <- us.gdp %>%
  filter(country == "United States") %>%
  mutate(across(starts_with("X"), as.character)) %>%
  pivot_longer(cols = starts_with("X"), 
               names_to = "Year", 
               values_to = "GDP")

us.gdp$Year <- as.numeric(gsub("X", "", us.gdp$Year))

us.gdp <- us.gdp %>%
  filter(Year >= 1933 & Year <= 2021) %>%
  mutate(GDP = as.numeric(sub("k", "", GDP)))

head(us.gdp)


## Employment Rate ----

us.emp <- read.csv("data/us.emp.csv")

us.emp <- us.emp %>%
  filter(country == "United States") %>%
  mutate(across(starts_with("X"), as.character)) %>%
  pivot_longer(cols = starts_with("X"), 
               names_to = "Year", 
               values_to = "Employment")

us.emp$Year <- as.numeric(gsub("X", "", us.emp$Year))

us.emp <- us.emp %>%
  filter(Year >= 1933 & Year <= 2021)

us.emp$Employment <- as.numeric(us.emp$Employment)

head(us.emp)


## Education ----

us.edu <- read.csv("data/us.edu.csv")

us.edu <- us.edu %>%
  filter(country == "United States") %>%
  mutate(across(starts_with("X"), as.character)) %>%
  pivot_longer(cols = starts_with("X"), 
               names_to = "year", 
               values_to = "Education")

us.edu$year <- as.numeric(gsub("X", "", us.edu$year))

us.edu <- us.edu %>%
  filter(year >= 1933 & year <= 2021)

us.edu$Education <- as.numeric(us.edu$Education)

head(us.edu)

edu.e0 <- merge(us.edu, us_ok, by = "year")



# Life table ----

compute_lifetable <- function(data, series, actual_year = NULL) {
  
  max.age = min(110, max(data$age))
  
  if (is.null(actual_year)) {
    
    lt <- lifetable(data, series = names(data$rate)[series], years = data$year,
                    ages = data$age, max.age = max.age, type = c("period"))
    print(lt)
    return(lt)
    
  } else {
    
    year_index <- actual_year - min(data$year) + 1
    
    lt <- lifetable(data, series = names(data$rate)[series], years = data$year[year_index],
                    ages = data$age, max.age = max.age, type = c("period"))
    print(lt)
    return(lt)
    
  }
  
}

lt_f <- compute_lifetable(us, 1) # females
lt_m <- compute_lifetable(us, 2) # males
lt_t <- compute_lifetable(us, 3) # total

write.table(lt_f, "output/lt_f.txt")
write.table(lt_m, "output/lt_m.txt")
write.table(lt_t, "output/lt_t.txt")


# Decomposing differences in life expectancy ----

e0_1960 <- life.expectancy(us, series = names(us$rate)[3],
                           years = us$year[28], type = c("period"))
e0_1980 <- life.expectancy(us, series = names(us$rate)[3],
                           years = us$year[48], type = c("period"))
delta_e0_1980_1960 <- e0_1980[1] - e0_1960[1]
delta_e0_1980_1960

ust_1960.lt <- lifetable(us, series = names(us$rate)[3],
                         years = us$year[28], ages = us$age,
                         max.age = min(110, max(us$age)), type = c("period"))
lx_1960 <- ust_1960.lt$lx
Lx_1960 <- ust_1960.lt$Lx


ust_1980.lt <- lifetable(us, series = names(us$rate)[3],
                         years = us$year[48], ages = us$age,
                         max.age = min(110, max(us$age)), type = c("period"))
lx_1980 <- ust_1980.lt$lx
Lx_1980 <- ust_1980.lt$Lx

delta00_04 <- (lx_1960[1]/lx_1960[1])  * (sum(Lx_1980[1:5])  /lx_1980[1]-sum(Lx_1960[1:5])/lx_1960[1])     + sum(Lx_1980[5:110]) /lx_1960[1]*(lx_1960[1]/lx_1980[1]-lx_1960[5]/lx_1980[5])
delta05_09 <- (lx_1960[6]/lx_1960[1])  * (sum(Lx_1980[6:10]) /lx_1980[6]-sum(Lx_1960[6:10])/lx_1960[6])    + sum(Lx_1980[10:110])/lx_1960[1]*(lx_1960[6]/lx_1980[6]-lx_1960[10]/lx_1980[10])

delta10_14 <- (lx_1960[11]/lx_1960[1]) * (sum(Lx_1980[11:15])/lx_1980[11]-sum(Lx_1960[11:15])/lx_1960[11]) + sum(Lx_1980[15:110])/lx_1960[1]*(lx_1960[11]/lx_1980[11]-lx_1960[25]/lx_1980[15])
delta15_19 <- (lx_1960[16]/lx_1960[1]) * (sum(Lx_1980[16:20])/lx_1980[16]-sum(Lx_1960[16:20])/lx_1960[16]) + sum(Lx_1980[20:110])/lx_1960[1]*(lx_1960[16]/lx_1980[16]-lx_1960[20]/lx_1980[20])

delta20_24 <- (lx_1960[21]/lx_1960[1]) * (sum(Lx_1980[21:25])/lx_1980[21]-sum(Lx_1960[21:25])/lx_1960[21]) + sum(Lx_1980[25:110])/lx_1960[1]*(lx_1960[21]/lx_1980[21]-lx_1960[25]/lx_1980[25])
delta25_29 <- (lx_1960[26]/lx_1960[1]) * (sum(Lx_1980[26:30])/lx_1980[26]-sum(Lx_1960[26:30])/lx_1960[26]) + sum(Lx_1980[30:110])/lx_1960[1]*(lx_1960[26]/lx_1980[26]-lx_1960[30]/lx_1980[30])

delta30_34 <- (lx_1960[31]/lx_1960[1]) * (sum(Lx_1980[31:35])/lx_1980[31]-sum(Lx_1960[31:35])/lx_1960[31]) + sum(Lx_1980[35:110])/lx_1960[1]*(lx_1960[31]/lx_1980[31]-lx_1960[35]/lx_1980[35])
delta35_39 <- (lx_1960[36]/lx_1960[1]) * (sum(Lx_1980[36:40])/lx_1980[36]-sum(Lx_1960[36:40])/lx_1960[36]) + sum(Lx_1980[40:110])/lx_1960[1]*(lx_1960[36]/lx_1980[36]-lx_1960[40]/lx_1980[40])

delta40_44 <- (lx_1960[41]/lx_1960[1]) * (sum(Lx_1980[41:45])/lx_1980[41]-sum(Lx_1960[41:45])/lx_1960[41]) + sum(Lx_1980[45:110])/lx_1960[1]*(lx_1960[41]/lx_1980[41]-lx_1960[45]/lx_1980[45])
delta45_49 <- (lx_1960[46]/lx_1960[1]) * (sum(Lx_1980[46:50])/lx_1980[46]-sum(Lx_1960[46:50])/lx_1960[46]) + sum(Lx_1980[50:110])/lx_1960[1]*(lx_1960[46]/lx_1980[46]-lx_1960[50]/lx_1980[50])

delta50_54 <- (lx_1960[51]/lx_1960[1]) * (sum(Lx_1980[51:55])/lx_1980[51]-sum(Lx_1960[51:55])/lx_1960[51]) + sum(Lx_1980[55:110])/lx_1960[1]*(lx_1960[51]/lx_1980[51]-lx_1960[55]/lx_1980[55])
delta55_59 <- (lx_1960[56]/lx_1960[1]) * (sum(Lx_1980[56:60])/lx_1980[56]-sum(Lx_1960[56:60])/lx_1960[56]) + sum(Lx_1980[60:110])/lx_1960[1]*(lx_1960[56]/lx_1980[56]-lx_1960[60]/lx_1980[60])

delta60_64 <- (lx_1960[61]/lx_1960[1]) * (sum(Lx_1980[61:65])/lx_1980[61]-sum(Lx_1960[61:65])/lx_1960[61]) + sum(Lx_1980[65:110])/lx_1960[1]*(lx_1960[61]/lx_1980[61]-lx_1960[65]/lx_1980[65])
delta65_69 <- (lx_1960[66]/lx_1960[1]) * (sum(Lx_1980[66:70])/lx_1980[66]-sum(Lx_1960[66:70])/lx_1960[66]) + sum(Lx_1980[70:110])/lx_1960[1]*(lx_1960[66]/lx_1980[66]-lx_1960[70]/lx_1980[70])

delta70_74 <- (lx_1960[71]/lx_1960[1]) * (sum(Lx_1980[71:75])/lx_1980[71]-sum(Lx_1960[71:75])/lx_1960[71]) + sum(Lx_1980[75:110])/lx_1960[1]*(lx_1960[71]/lx_1980[71]-lx_1960[75]/lx_1980[75])
delta75_79 <- (lx_1960[66]/lx_1960[1]) * (sum(Lx_1980[76:80])/lx_1980[76]-sum(Lx_1960[76:80])/lx_1960[76]) + sum(Lx_1980[80:110])/lx_1960[1]*(lx_1960[76]/lx_1980[76]-lx_1960[80]/lx_1980[80])

delta80_84 <- (lx_1960[81]/lx_1960[1]) * (sum(Lx_1980[81:85])/lx_1980[81]-sum(Lx_1960[81:85])/lx_1960[81]) + sum(Lx_1980[85:110])/lx_1960[1]*(lx_1960[81]/lx_1980[81]-lx_1960[85]/lx_1980[85])
delta85piu <- (lx_1960[86]/lx_1960[1]) * (sum(Lx_1980[86:110])/lx_1980[86]-sum(Lx_1960[86:106])/lx_1960[86])

delta <- c(delta00_04, delta05_09, delta10_14, delta15_19, delta20_24, delta25_29, delta30_34, delta35_39, delta40_44, delta45_49,
           delta50_54, delta55_59, delta60_64, delta65_69, delta70_74, delta75_79, delta80_84, delta85piu)

classi <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49",
            "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85+")

delta_tot <- sum(delta)
delta_perc <- delta/delta_tot*100

tab_delta <- matrix(rbind(delta, delta_perc), nrow = 18, ncol = 2, byrow = TRUE, dimnames = list(c(classi), c("5 Delta x", "5 Delta x %")))
tab_delta

ggplot(data = as.data.frame(delta), aes(x = classi, y = delta, fill = ifelse(delta < 0, "#FF9F7F", "#1f77b4"))) +
  geom_bar(stat = "identity", col = "black") +
  labs(x = "Age Class", y = "Difference", title = "Decomposition of differences in life expectancy between 1960 and 1980 by age class") +
  scale_fill_manual(values = c("#FF9F7F" = "#FF9F7F", "#1f77b4" = "#1f77b4")) +
  theme_minimal() +
  guides(fill = FALSE)



# GDP ----

ggplot(data = us.gdp, aes(x = Year, y = GDP)) +
  geom_line(color = "#FF9F7F", linewidth = 1) +
  geom_vline(xintercept = c(1945, 1960, 1980, 1990, 2019), linetype = "dashed", color = "darkgreen") +
  xlim(1933, 2021) +
  labs(x = "Year", y = "GDP in thousand $", title = "GDP per capita (price and inflation adjusted), 1933-2021") +
  theme_minimal()


## Employment rate ----

ggplot(data = us.emp, aes(x = Year, y = Employment)) +
  geom_line(col = "#FF9F7F", linewidth = 1) +
  labs(x = "Year", y = "Employment Rate", title = "Trend of employment rate over years - US, aged 15+") +
  theme_minimal()


## Education and e0 ----

ggplot(data = edu.e0, aes(x = Education)) +
  geom_point(aes(y = Female, color = "Female")) +
  geom_point(aes(y = Male, color = "Male")) +
  geom_smooth(aes(y = Female, color = "Female"), method = "lm", se = FALSE) +
  geom_smooth(aes(y = Male, color = "Male"), method = "lm", se = FALSE) +
  labs(x = "OWID Education Index", y = "Life Expectancy", color = "Sex",  title = "Relationship between life expectancy and educational attainment") +
  scale_color_manual(values = c("#FF9F7F", "#1f77b4"), labels = c("Female", "Male")) +
  theme_minimal()


# Age pyramid ----

plot_pyramid <- function(data, actual_year) {
  
  year_index <- actual_year - min(data$year) + 1
  
  females <- data$pop$female[, year_index]
  males <- data$pop$male[, year_index]
  ages <- data$age
  pyr_data <- data.frame(males, females, ages)
  
  pyramid(pyr_data, Llab = "Males", Rlab = "Females", Clab = "Age",
          Laxis = seq(0, 3000000, by = 1000000),
          AxisFM = "d", AxisBM = ",",
          Csize = 0.8, Cstep = 10,
          main = paste("Age Pyramid - USA, ", actual_year))
}

par(mfrow = c(1,3))
plot_pyramid(us, 1933)
plot_pyramid(us, 1977)
plot_pyramid(us, 2021)


# e0 trend

ggplot(data = us_ok_long, aes(x = year, y = Value, group = sex, colour = sex)) +
  geom_line(size = 1) +
  geom_vline(xintercept = c(1945, 1960, 1980, 1990, 2019), linetype = "dashed", color = "darkgreen") +
  scale_colour_manual(values = c("#FF9F7F", "#1f77b4"), name = "Sex",
                      breaks = c("Female", "Male"), labels = c("Female", "Male")) +
  labs(x = "Year", y = "Life expectancy at birth", title = "Life expectancy at birth, 1933-2021") +
  theme_minimal()


# TFR trend

ggplot(data = us.tfr, aes(x = Year, y = TFR)) +
  geom_line(color = "#FF9F7F", size = 1) +
  geom_vline(xintercept = c(1945, 1960, 1980, 1990, 2019), linetype = "dashed", color = "darkgreen") +
  ylim(0, 4) +
  xlim(1933, 2021) +
  labs(x = "Year", y = "TFR", title = "Total Fertility Rate, 1933-2021") +
  theme_minimal()



# Deaths dx ----

plot_deaths <- function(data, years, series, main) {
  
  par(mfrow = c(1,2))
  lt <- lifetable(data, series = names(data$rate)[series])
  age <- lt$age
  year <- lt$year
  
  persp(age, year, lt$dx, theta = 35, zlab = "dx",
        main = paste("Deaths dx,", min(data$year), "-", max(data$year)))
  
  lt_list <- list()
  colors <- brewer.pal(length(years), "Set1")  
  
  for (i in 1:length(years)) {
    
    year_index <- years[i] - min(data$year) + 1
    lt <- lifetable(data, series = names(data$rate)[series],
                    years = data$year[year_index], ages = data$age,
                    max.age = min(110, max(data$age)),
                    type = c("period"))
    lt_list[[i]] <- lt
    
  }
  
  plot(NULL, xlim = c(0, 119), ylim = c(0, 0.05), lwd = 2, xlab = "Age", ylab = "dx", main = main, bty = "n")
  
  for (i in 1:length(years)) {
    
    dx <- lt_list[[i]]$dx
    lines(dx, col = colors[i], lty = 1, lwd = 2)
    
  }
  
  legend("topright", legend = as.character(years), col = colors[1:length(years)], lty = 1, lwd = 2, cex = 0.9, bty = "n")
  grid(col = "#F5F5F5", lty = 1)
  
}

plot_deaths(us, c(1945, 1960, 1985, 2019, 2021), 3, main = "Deaths dx")



# Survivors lx ----

plot_survivors <- function(data, years, series, main) {
  
  par(mfrow = c(1,2))
  lt <- lifetable(data, series = names(data$rate)[series])
  age <- lt$age
  year <- lt$year
  
  persp(age, year, lt$lx, theta = 35, zlab = "lx",
        main = paste("Survivors lx,", min(data$year), "-", max(data$year)))
  
  
  lt_list <- list()
  colors <- brewer.pal(length(years), "Set1")  
  
  for (i in 1:length(years)) {
    
    year_index <- years[i] - min(data$year) + 1
    lt <- lifetable(data, series = names(data$rate)[series],
                    years = data$year[year_index], ages = data$age,
                    max.age = min(110, max(data$age)),
                    type = c("period"))
    lt_list[[i]] <- lt
    
  }
  
  plot(NULL, xlim = c(0, 119), ylim = c(0, 1), lwd = 2, xlab = "Age", ylab = "lx", main = main, bty = "n")
  
  for (i in 1:length(years)) {
    
    lx <- lt_list[[i]]$lx
    lines(lx, col = colors[i], lty = 1, lwd = 2)
    
  }
  
  legend("topright", legend = as.character(years), col = colors[1:length(years)], lty = 1, lwd = 2, cex = 0.9, bty = "n")
  grid(col = "#F5F5F5", lty = 1)
  
}

plot_survivors(us, c(1945, 1960, 1985, 2019, 2021), 3, main = "Survivors lx")




## Mortality rates mx ----

par(mfrow = c(1,3))

years = c(1933, 1945, 1960, 1985, 2019, 2021)

plot(us, series = names(us$rate)[1], years = years,
     main = "Mortality rates - USA, Females")
legend("bottomright", legend = years,
       col = rainbow(length(years)), lty = 1, cex = 0.9, bty ="n") 

plot(us, series = names(us$rate)[2], years = years,
     main = "Mortality rates - USA, Males")
legend("bottomright", legend = years,
       col = rainbow(length(years)), lty = 1, cex = 0.9, bty ="n") 

plot(us, series = names(us$rate)[3], years = years,
     main = "Mortality rates - USA, Total")
legend("bottomright", legend = years,
       col = rainbow(length(years)), lty = 1, cex = 0.9, bty ="n") 



## Life expectancy by age ----

years = c(1945, 1960, 1985, 2019, 2021)

par(mfrow = c(1,3))
plot(lt_f, years = years, xlim = c(60,110), ylim = c(0, 22), main = "Life Expectancy - USA, 1946-2020, Females")
legend("topright", legend = years,
       col = rainbow(length(years)), lty = 1, cex = 0.9, bty ="n") 

plot(lt_m, years = years, xlim = c(60,110), ylim = c(0, 22), main = "Life Expectancy - USA, 1946-2020, Males")
legend("topright", legend = years,
       col = rainbow(length(years)), lty = 1, cex = 0.9, bty ="n") 

plot(lt_t, years = years, xlim = c(60,110), ylim = c(0, 22), main = "Life Expectancy - USA, 1946-2020, Total")
legend("topright", legend = years,
       col = rainbow(length(years)), lty = 1, cex = 0.9, bty ="n") 




# Lee-Carter model and forecasting of e0 ----

run_country_analysis <- function(data, series, years, h) {
  
  ## Lee-Carter model ----
  LC_model <- lca(data, series = names(data$rate)[series], years = years)
  LC_model
  plot(LC_model)
  
  ## Residual analysis ----
  par(mfrow = c(1,3))
  residuals <- residuals(LC_model, "residuals")
  plot(rep(residuals$y, length(residuals$x)), residuals$z, xlab = "age", ylab = "residuals",
       main = "Residuals by ages from LC model")
  plot(rep(residuals$x, length(residuals$y)), residuals$z, xlab = "years", ylab = "residuals",
       main = "Residuals by years from LC model")
  plot(residuals, type = "image",
       main = "Residuals heatmap")
  
  ## Forecast Kt as a random walk with drift ARIMA(0,1,0) ----
  par(mfrow = c(1,1))
  forecast <- forecast(LC_model)
  plot(forecast$kt, bty = "n", main = "Forecasted Kt component")
  
  ## Compare observed, (LC fitted,) forecasted death rates ----
  par(mfrow = c(1,2))
  #plot(fitted(LC_model), main = paste("USA: LC total death rates (", years[1], "-", years[length(years)], ")"))
  plot(data, series = names(data$rate)[series], years = data$year, ylim = c(-10, 0), bty = "n",
       main = paste("Death rates,", min(years), "-", max(years)))
  grid(col = "#F5F5F5", lty = 1)
  plot(forecast,  ylim = c(-10, 0), bty = "n", 
       main = paste("Projected death rates,", max(years)+1, "-", max(years)+h))
  grid(col = "#F5F5F5", lty = 1)
  
  ## Projected lifetable
  lifetable_forecast <- lifetable(forecast) # ,years = 2050

  ## Life expectancy forecast with prediction interval ----
  e0_forecast <- e0(forecast, PI = TRUE, nsim = 200)
  par(mfrow = c(1,1))
  plot(e0_forecast, bty = "n", main = "Past and Forecasted Life Expectancy")
  
  invisible(NULL)
}

run_country_analysis(us, 3, us$year, 50)



# LC and RH models ----
Ext <- us$pop$male
Dxt <- us$rate$male * Ext
ages <- us$age
years = us$year

LC <- lc()
RH <- rh()

LCfit <- fit(LC, Dxt = Dxt, Ext = Ext, ages = ages, years = years, ages.fit = 60:89)
RHfit <- fit(RH, Dxt = Dxt, Ext = Ext, ages = ages, years = years, ages.fit = 60:89)

plot(LCfit)
plot(RHfit)

## Residual analysis ----
LCres <- residuals(LCfit)
plot(LCres, type = 'scatter')
plot(LCres, type = 'colourmap')

RHres <- residuals(RHfit)
plot(RHres, type = 'scatter')
plot(RHres, type = 'colourmap')


## Comparison ----
BIC <- matrix (nrow = 1, ncol = 2)
rownames(BIC) <- "BIC"
colnames(BIC) <- c("LC", "RH")
BIC[1,1] <- BIC(LCfit)
BIC[1,2] <- BIC(RHfit)
BIC

## Forecasting ----
LCfor <- forecast(LCfit, h = 50)
plot(LCfor, only.kt = TRUE)

RHfor <- forecast(RHfit, h = 50)
plot(RHfor, only.kt = TRUE)
plot(RHfor, only.gc = TRUE)


# Detrend Life expectancy ----

par(mfrow = c(1,1))
life_expect <- life.expectancy(us, series = names(us$rate)[3],  years = us$year, type = c("period"))
plot(life_expect)

par(mfrow = c(1, 2), mar = c(2.2, 2.2, 1, 1), cex = 0.8)

# Detrend data with a linear filter ----
lin.mod <- lm(life_expect ~ us$year)
lin.trend <- lin.mod$fitted.values
linear <- ts(lin.trend, start = c(1933, 1), frequency = 1)
lin.residual <- life_expect - linear

plot.ts(life_expect, ylab = "", main = "Life expectacny")
lines(linear, col = "red")
legend("bottomright", legend = c("data", "trend"), lty = 1, col = c("black", "red"), bty = "n")
plot.ts(lin.residual, ylab = "", main = "Residual trend")
legend("bottomright", legend = c("residuals"), lty = 1, col = c("black"), bty = "n")


# Detrend data with the Hodrick-Prescott filter ----

hp.decom <- hpfilter(life_expect, freq = 6.25, type = "lambda")

plot.ts(life_expect, ylab = "", main = "Life expectancy")
lines(hp.decom$trend, col = "red")
legend("bottomright", legend = c("data", "HPtrend"), lty = 1, col = c("black", "red"), bty = "n")
plot.ts(hp.decom$cycle, ylab = "", main = "Residual trend")
legend("bottomright", legend = c("HPresidual"), lty = 1, col = c("black"), bty = "n")


# Detrend data with the Baxter-King filter ----

bk.decom <- bkfilter(life_expect)

plot.ts(life_expect, ylab = "", main = "Life expectancy")
lines(bk.decom$trend, col = "red")
legend("bottomright", legend = c("data", "BKtrend"), lty = 1, col = c("black", "red"), bty = "n")
plot.ts(bk.decom$cycle, ylab = "", main = "Residual trend")
legend("bottomright", legend = c("BKresidual"), lty = 1, col = c("black"), bty = "n")





