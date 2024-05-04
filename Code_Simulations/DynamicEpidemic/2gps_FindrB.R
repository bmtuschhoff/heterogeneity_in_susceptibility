## read in values of cv, f_inf, f_a and from that find r_b
## then output r_b so can be piped into C++ code for dynamic contact tracing data in discrete case

arg <- commandArgs(trailingOnly=T)

cv <- as.numeric(arg[1])
f_inf <- as.numeric(arg[2])
f_a_old <- as.numeric(arg[3])

# function to find r_b with root finding
find_r_b <- function(r_b, f_a, cv, f_inf)
{
  abs(((-log(1 - (1/f_a)*(f_inf - (1 - exp(-r_b))*(1 - f_a))) - r_b)*sqrt(f_a*(1 - f_a))) / ((-log(1 - (1/f_a)*(f_inf - (1 - exp(-r_b))*(1 - f_a))))*f_a + r_b*(1 - f_a)) - cv)
}


if(cv == 0)
{
  rb <- -log(1 - f_inf)
} else
{
  # find risk values numerically solving for root
  if(f_inf < f_a_old)
  {
    if(f_inf < 1-f_a_old)
    {
      rb <- optimize(find_r_b, interval=c(0,-log(1 - (f_inf / (1 - f_a_old)))), f_a=f_a_old, cv=cv, f_inf=f_inf)$minimum
    } else
    {
      rb <- optimize(find_r_b, interval=c(0,50), f_a=f_a_old, cv=cv, f_inf=f_inf)$minimum
    }
  } else
  {
    if(f_inf < 1-f_a_old)
    {
      rb <- optimize(find_r_b, interval=c(-log(1 - ((f_inf - f_a_old) / (1 - f_a_old))),-log(1 - (f_inf / (1 - f_a_old)))), f_a=f_a_old, cv=cv, f_inf=f_inf)$minimum
    } else
    {
      rb <- optimize(find_r_b, interval=c(-log(1 - ((f_inf - f_a_old) / (1 - f_a_old))),50), f_a=f_a_old, cv=cv, f_inf=f_inf)$minimum
    }
  }
}

cat(rb)