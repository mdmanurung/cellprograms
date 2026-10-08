## summarize_sim_adj.R <dir of rep*.csv>  -> per scenario x arm: rejection rate (p_pair<0.05), batch corr, recovery
d <- do.call(rbind, lapply(list.files(commandArgs(TRUE)[1], "^rep.*csv$", full.names = TRUE), read.csv))
d$rej <- d$p_pair < 0.05
s <- aggregate(cbind(rej, max_cor_batch, best_cor_z, k_a) ~ shared + rho + arm, d, mean)
s$n <- aggregate(rep ~ shared + rho + arm, d, function(v) length(unique(v)))$rep
s[order(s$shared, s$rho, s$arm), ] |> format(digits = 2) |> print(row.names = FALSE)
