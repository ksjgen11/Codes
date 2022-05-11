num <- round(runif(1)*100, digits = 0)
guess <- -1
cat('Guess a number between 0 and 100.\n')

while (guess != num){
  guess <- readline(prompt = 'Guess number: ')
  guess <- as.integer(guess)
  if (guess == num){
    cat("Congraturations! ", num, " is right\n")
  } else if (guess < num){
    cat("It's smaller!\n")
  } else if (guess >num){
    cat("It's bigger!\n")
  }
}

x <- 1
y <- 1

while (x != 0 | y != 0){
  x <- readline(prompt = '1: ')
  x <- as.integer(x)
  y <- readline(prompt = '2: ')
  y <- as.integer(y)
  cat (x*y)
}