function val = objective(args, model, MALags, w)
val = loglikelihood(args, model, MALags, w);
val = -val;

