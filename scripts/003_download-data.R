dat = read_csv(url)
linelist = NCoVUtils::get_linelist()

state_abbrev = 
  read_csv("https://raw.githubusercontent.com/jasonong/List-of-US-States/master/states.csv") %>% 
  mutate_at(vars(State), to_snake_case) %>% 
  setNames(tolower(colnames(.)))