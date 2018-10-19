options(stringsAsFactors = FALSE)

dirw = file.path(Sys.getenv("misc2"), "grn23", "51.camoco")

### convert camoco coex table to R symmetric matrix
fc = file.path(dirw, "coex.csv")
vec = scan(fc, what = numeric())

roots = polyroot(c(-length(vec), -0.5, 0.5))
ng = Re(roots)[Re(roots) > 0]
stopifnot(ng*(ng-1)/2 == length(vec))

coex <- matrix(rep(0, ng*ng), nrow=ng)
coex[lower.tri(coex)] = vec
coex = t(coex)
coex[lower.tri(coex)] = vec
save(coex, file = file.path(dirw, "camoco.rda"))


