args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE))))
root <- normalizePath(file.path(here, "../.."))
for (f in list.files(file.path(root, "R"), full.names = TRUE)) source(f)
for (f in c("sim_coord.R", "evaluate.R")) source(file.path(here, f))
w <- matrix(c(.9, .85, .8, 0), 2, 2)
stopifnot(identical(coord_assign(w), c(2L, 1L)))
stopifnot(all(is.na(coord_assign(matrix(.2, 3, 2)))))
d <- gen_coord("act_only", c(R100=1,R50=.8,R20=.6), 41, G=2000, r20_size="small")
f <- coord_truth_fit(d)
map <- coord_match_nodes(f, d)
stopifnot(all(map$program_id == map$truth), all(d$truth$programs$R20$members %in% d$truth$small))
g <- do.call(rbind, lapply(seq_along(d$truth$programs), function(i) {
 p <- d$truth$programs[[i]]
 data.frame(cluster=i, cell_type=p$members, program_id=p$id)
}))
m <- coord_evaluate(list(clusters=g), f, d, "ideal", oracle=TRUE)
stopifnot(all(m$found), all(m$recall==1), all(m$precision==1), !any(m$spurious))
gene <- gen_coord("genes_only", 1, 42, G=2000)
gf <- coord_truth_fit(gene)
m <- coord_evaluate(list(clusters=g), gf, gene, "gene_match", oracle=TRUE)
stopifnot(!any(m$present), all(m$spurious))
d07 <- gen_coord("act_only", 1, 43, G=2000, act_cor=.7)
r <- cor(do.call(cbind, d07$truth$programs$R100$z))
stopifnot(abs(mean(r[upper.tri(r)])-.7)<.12)
cat("Evaluator independent-oracle checks passed\n")
# Duplicate half-clusters must not both earn a recovery; they are fragments.
g2 <- g
ix <- g2$program_id == "R100"
g2$cluster[which(ix)[1:5]] <- 4L
m <- coord_evaluate(list(clusters=g2), f, d, "ideal", oracle=TRUE)
stopifnot(all(m$spurious), all(m$fragments==1L), m$recall[m$program=="R100"]==.5)
# A nuisance-like loading cannot steal a perfectly matched recurring activity.
f2 <- f; ct <- d$truth$programs$R100$members[1]
f2$loadings[[ct]][,"R100"] <- d$truth$loadings[[ct]][,"N"]
mm <- coord_match_nodes(f2,d)
stopifnot(mm$truth[mm$cell_type==ct & mm$program_id=="R100"]=="R100")
