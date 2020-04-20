library(tidyverse)
library(KFAS)
library(zoo)
library(snakecase)
library(tictoc)
library(lubridate)
library(rstan)

for(f in list.files("R", full.names = T)) source(f)
