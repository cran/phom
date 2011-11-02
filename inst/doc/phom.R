### R code from vignette source 'phom.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ex1
###################################################
library(phom)

N <- 50

x1 <- rnorm(N) * 0.1
y1 <- rnorm(N) * 0.1
X1 <- t(as.matrix(rbind(x1, y1)))

x2 <- rnorm(N) * 0.1 + 0.5
y2 <- rnorm(N) * 0.1 + 0.5
X2 <- t(as.matrix(rbind(x2, y2)))

x <- cbind(x1, x2)
y <- cbind(y1, y2)

X <- as.matrix(rbind(X1, X2))

max_dim <- 0
max_f <- 0.8

intervals <- pHom(X, max_dim, max_f, metric="manhattan")



###################################################
### code chunk number 2: fig1plot
###################################################
plot(X)


###################################################
### code chunk number 3: fig1
###################################################
plot(X)


###################################################
### code chunk number 4: pd1plot
###################################################
plotBarcodeDiagram(intervals, max_dim, max_f, title="")


###################################################
### code chunk number 5: pd1
###################################################
plotBarcodeDiagram(intervals, max_dim, max_f, title="")


###################################################
### code chunk number 6: ex2
###################################################
library(phom)

t <- 2 * pi * runif(100)
x <- cos(t); y <- sin(t)
X <- t(as.matrix(rbind(x, y)))

max_dim <- 1
max_f <- 0.6

intervals <- pHom(X, max_dim, max_f)


###################################################
### code chunk number 7: fig2plot
###################################################
plot(X)


###################################################
### code chunk number 8: fig2
###################################################
plot(X)


###################################################
### code chunk number 9: pd2plot
###################################################
plotPersistenceDiagram(intervals, max_dim, max_f, 
	title="Random Points on S^1 with l_2 Norm")


###################################################
### code chunk number 10: pd2
###################################################
plotPersistenceDiagram(intervals, max_dim, max_f, 
	title="Random Points on S^1 with l_2 Norm")


###################################################
### code chunk number 11: ex3
###################################################
library(phom)

sphere_points <- function(n, d) {
	points <- matrix(rnorm(n*d), nrow = n, ncol = d)
	L <- apply(points, MARGIN = 1,
		   FUN = function(x){sqrt(sum(x*x))})
	D <- diag(1 / L)
	U <- D %*% points
	U
}

n <- 5000
d <- 4
X <- as.matrix(sphere_points(n, d))

max_dim <- d - 1
max_f <- 0.9
landmark_set_size <- 40
maxmin_samples <- 1000

intervals <- pHom(X, max_dim, max_f, mode="lw", metric="euclidean", 
	landmark_set_size=landmark_set_size, maxmin_samples=maxmin_samples)


###################################################
### code chunk number 12: pd3plot
###################################################
plotPersistenceDiagram(intervals, max_dim, max_f, 
	title="Random Points on S^3")


###################################################
### code chunk number 13: pd3
###################################################
plotPersistenceDiagram(intervals, max_dim, max_f, 
	title="Random Points on S^3")


