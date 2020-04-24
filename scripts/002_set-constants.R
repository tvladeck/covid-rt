url = 'https://covidtracking.com/api/v1/states/daily.csv'

WINDOW = 20
SERIAL_INTERVAL = 4
GAMMA = 1 / SERIAL_INTERVAL

shutdown_dates = 
  list(
    alabama = ymd(20200404),
    alaska = ymd(20200328),
    california = ymd(20200319),
    colorado = ymd(20200326),
    connecticut = ymd(20200323),
    delaware = ymd(20200324),
    district_of_columbia = ymd(20200401),
    florida = ymd(20200403),
    georgia = ymd(20200403),
    hawaii = ymd(20200325),
    idaho = ymd(20200325),
    illinois = ymd(20200321),
    indiana = ymd(20200324),
    kansas = ymd(20200330),
    kentucky = ymd(20200326),
    louisiana = ymd(20200323),
    maine = ymd(20200402),
    maryland = ymd(20200330),
    massachusetts = ymd(20200324),
    michigan = ymd(20200324),
    minnesota =ymd(20200331),
    mississippi = ymd(20200403),
    missouri = ymd(20200406),
    montana = ymd(20200328),
    nevada = ymd(20200401),
    new_hampshire = ymd(20200327),
    new_jersey = ymd(20200321),
    new_mexico = ymd(20200322)
  )