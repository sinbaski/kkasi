function val = objective(args, w)
val = loglikelihood(args, w);
val = -val;

